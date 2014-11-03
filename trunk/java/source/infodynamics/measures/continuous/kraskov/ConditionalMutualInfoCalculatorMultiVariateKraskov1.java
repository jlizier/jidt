/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.measures.continuous.kraskov;

import java.util.Calendar;
import java.util.PriorityQueue;

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.NeighbourNodeData;

/**
 * <p>Computes the differential conditional mutual information of two multivariate
 *  <code>double[][]</code> sets of observations, conditioned on another
 *  (implementing {@link ConditionalMutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see references below)
 *  <b>algorithm 1</b>.
 *  Most of the functionality is defined by the parent class 
 *  {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}.</p>
 *
 * <p>Crucially, the calculation is performed by examining
 * neighbours in the full joint space (as specified by Frenzel and Pompe)
 * rather than two MI calculators.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link ConditionalMutualInfoCalculatorMultiVariate},
 * and expanded on in {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}.
 * </p>
 *  
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Frenzel and Pompe, <a href="http://dx.doi.org/10.1103/physrevlett.99.204101">
 * 	"Partial Mutual Information for Coupling Analysis of Multivariate Time Series"</a>,
 * 	Physical Review Letters, <b>99</b>, p. 204101+ (2007).</li>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @author Ipek Özdemir
 */
public class ConditionalMutualInfoCalculatorMultiVariateKraskov1
	extends ConditionalMutualInfoCalculatorMultiVariateKraskov {
	
	public ConditionalMutualInfoCalculatorMultiVariateKraskov1() {
		super();
		isAlgorithm1 = true;
	}

	@Override
	protected double[] partialComputeFromObservations(int startTimePoint,
			int numTimePoints, boolean returnLocals) throws Exception {
		
		double startTime = Calendar.getInstance().getTimeInMillis();

		double[] localCondMi = null;
		if (returnLocals) {
			localCondMi = new double[numTimePoints];
		}
		
		// Count the average number of points within eps_xz and eps_yz and eps_z of each point
		double sumDiGammas = 0;
		double sumNxz = 0;
		double sumNyz = 0;
		double sumNz = 0;
		
		// Arrays used for fast searching on conditionals with a marginal:
		boolean[] isWithinRForConditionals = new boolean[totalObservations];
		int[] indicesWithinRForConditionals = new int[totalObservations+1];

		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps for this time step by
			//  finding the kth closest neighbour for point t:
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreeJoint.findKNearestNeighbours(k, t);
			// First element in the PQ is the kth NN,
			//  and epsilon = kthNnData.distance
			NeighbourNodeData kthNnData = nnPQ.poll();
			
			// Now count the points in the conditional space, and 
			//  the var1-conditional and var2-conditional spaces.
			// We have 3 coded options for how to do this:
			
			/* Option A -- straightforward way using each k-d tree separately:
			int n_xz = kdTreeVar1Conditional.countPointsStrictlyWithinR(
					t, kthNnData.distance);
			int n_yz = kdTreeVar2Conditional.countPointsStrictlyWithinR(
					t, kthNnData.distance);
			int n_z = nnSearcherConditional.countPointsStrictlyWithinR(
			t, kthNnData.distance);
			*/ // end option A
			
			/* Option B -- 
			 * Select all points within conditional z, then check x and y norms for
			 *  these points only. Works better than A if few points qualify for the conditionals
			 *  (e.g. large multivariate conditional) but not so well for 
			 *  many qualifying points (e.g. low dimensional conditional).
			 *   
			Collection<NeighbourNodeData> z_pointsWithinR = 
					nnSearcherConditional.findPointsStrictlyWithinR(
							t, kthNnData.distance);
			int n_z = z_pointsWithinR.size();
			
			int n_xz = 0, n_yz = 0;
			for(NeighbourNodeData zNeighbour :  z_pointsWithinR) {
				if (KdTree.normWithAbort(
						var1Observations[t],
						var1Observations[zNeighbour.sampleIndex],
						kthNnData.distance, normType) < kthNnData.distance) {
					n_xz++;
				}
				if (KdTree.normWithAbort(
						var2Observations[t],
						var2Observations[zNeighbour.sampleIndex],
						kthNnData.distance, normType) < kthNnData.distance) {
					n_yz++;
				}
			}
			 */ // end option B
			
			// Option C --
			// Identify the points satisfying the conditional criteria, then use
			//  the knowledge of which points made this cut to speed up the searching
			//  in the conditional-marginal spaces:
			// 1. Identify the n_z points within the conditional boundaries:
			nnSearcherConditional.findPointsWithinR(t, kthNnData.distance,
					false, isWithinRForConditionals, indicesWithinRForConditionals);
			// 2. Then compute n_xz and n_yz harnessing our knowledge of
			//  which points qualified for the conditional already:
			int n_xz = kdTreeVar1Conditional.countPointsWithinR(t, kthNnData.distance,
					false, 1, isWithinRForConditionals);
			int n_yz = kdTreeVar2Conditional.countPointsWithinR(t, kthNnData.distance,
					false, 1, isWithinRForConditionals);
			// 3. Finally, reset our boolean array for its next use while we count n_z:
			int n_z;
			for (n_z = 0; indicesWithinRForConditionals[n_z] != -1; n_z++) {
				isWithinRForConditionals[indicesWithinRForConditionals[n_z]] = false;
			}
			// end option C

			sumNxz += n_xz;
			sumNyz += n_yz;
			sumNz += n_z;
			// And take the digammas:
			double digammaNxzPlusOne = MathsUtils.digamma(n_xz+1);
			double digammaNyzPlusOne = MathsUtils.digamma(n_yz+1);
			double digammaNzPlusOne = MathsUtils.digamma(n_z+1);
			sumDiGammas += digammaNzPlusOne - digammaNxzPlusOne - digammaNyzPlusOne;
			
			if (returnLocals) {
				localCondMi[t-startTimePoint] = digammaK - digammaNxzPlusOne - digammaNyzPlusOne + digammaNzPlusOne;
				if (debug) {
					System.out.printf("t=%d, n_xz=%d, n_yz=%d, n_z=%d, local=%.4f\n",
							t, n_xz, n_yz, n_z, localCondMi[t-startTimePoint]);
				}
			}
		}
		
		if (debug) {
			Calendar rightNow2 = Calendar.getInstance();
			long endTime = rightNow2.getTimeInMillis();
			System.out.println("Subset " + startTimePoint + ":" +
					(startTimePoint + numTimePoints) + " Calculation time: " +
					((endTime - startTime)/1000.0) + " sec" );
		}

		// Select what to return:
		if (returnLocals) {
			return localCondMi;
		} else {
			// Pad return array with two values, to allow compatibility in 
			//  return length with algorithm 2
			return new double[] {sumDiGammas, sumNxz, sumNyz, sumNz, 0, 0};
		}		
	}
}
