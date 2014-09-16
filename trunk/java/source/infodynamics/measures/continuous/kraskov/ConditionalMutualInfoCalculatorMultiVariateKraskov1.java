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

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

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
	
	/**
	 * Multiplier used as hueristic for determining whether to use a linear search
	 *  for kth nearest neighbour or a binary search.
	 */
	protected static final double CUTOFF_MULTIPLIER = 1.5;

	public ConditionalMutualInfoCalculatorMultiVariateKraskov1() {
		super();
		isAlgorithm1 = true;
	}

	@Override
	protected double[] partialComputeFromObservations(int startTimePoint,
			int numTimePoints, boolean returnLocals) throws Exception {
		
		double startTime = Calendar.getInstance().getTimeInMillis();

		int N = var1Observations.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		
		double[] localCondMi = null;
		if (returnLocals) {
			localCondMi = new double[numTimePoints];
		}
		
		// Constants:
		double digammaK = MathsUtils.digamma(k);
		
		// Count the average number of points within eps_xz and eps_yz and eps_z of each point
		double sumDiGammas = 0;
		double sumNxz = 0;
		double sumNyz = 0;
		double sumNz = 0;
		
		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps for this time step:
			//  First get x and y and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity.
			double[][] xyzNorms = normCalculator.computeNorms(
					var1Observations, var2Observations, condObservations, t);
			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2] = Math.max(xyzNorms[t2][0], Math.max(xyzNorms[t2][1], xyzNorms[t2][2]));
			}
			// Then find the kth closest neighbour, using a heuristic to 
			// select whether to keep the k mins only or to do a sort.
			double epsilon = 0.0;
			if (k <= cutoffForKthMinLinear) {
				// just do a linear search for the minimum
				epsilon = MatrixUtils.kthMin(jointNorm, k);
			} else {
				// Sort the array of joint norms first
				java.util.Arrays.sort(jointNorm);
				// And find the distance to it's kth closest neighbour
				// (we subtract one since the array is indexed from zero)
				epsilon = jointNorm[k-1];
			}
			
			// Count the number of points whose z distance is less
			//  than eps, x and z distance is less than eps, and
			//  y and z distance is less than eps
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyzNorms[t2][2] < epsilon) {
					n_z++;
					if (xyzNorms[t2][0] < epsilon) {
						n_xz++;
					}
					if (xyzNorms[t2][1] < epsilon) {
						n_yz++;
					}
				}
			}
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

	
	@Override
	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e",
				k, MathsUtils.digamma(k));
		return constants;
	}
}
