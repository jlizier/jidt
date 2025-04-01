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

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.NeighbourNodeData;

/**
 * <p>Computes the differential mutual information of two given multivariate sets of
 *  observations (implementing {@link MutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below),
 *  <b>algorithm 1</b>.
 *  Most of the functionality is defined by the parent class 
 *  {@link MutualInfoCalculatorMultiVariateKraskov}.</p>
 *
 * <p>Usage is as per the paradigm outlined for {@link MutualInfoCalculatorMultiVariate},
 * and expanded on in {@link MutualInfoCalculatorMultiVariateKraskov}.
 * </p>
 *  
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @author Ipek Özdemir
 */
public class MutualInfoCalculatorMultiVariateKraskov1
	extends MutualInfoCalculatorMultiVariateKraskov {
	
	public MutualInfoCalculatorMultiVariateKraskov1() {
		super();
		isAlgorithm1 = true;
	}

	@Override
	protected double[] partialComputeFromObservations(
			int startTimePoint, int numTimePoints, boolean returnLocals) throws Exception {
		
		double startTime = Calendar.getInstance().getTimeInMillis();

		double[] localMi = null;
		if (returnLocals) {
			localMi = new double[numTimePoints];
		}
		
		// Count the average number of points within eps_x and eps_y of each point
		double sumDiGammas = 0;
		double sumNx = 0;
		double sumNy = 0;
		double sum2xkNNDist = 0;
		
		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps for this time step by
			//  finding the kth closest neighbour for point t:
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreeJoint.findKNearestNeighbours(k, t, dynCorrExclTime);
			// First element in the PQ is the kth NN,
			//  and epsilon = kthNnData.distance
			NeighbourNodeData kthNnData = nnPQ.poll();

			// Count the number of points whose x distance is less
			//  than eps, and whose y distance is less than
			//  epsilon = kthNnData.distance
			int n_x = nnSearcherSource.countPointsStrictlyWithinR(
							t, kthNnData.distance, dynCorrExclTime);
			int n_y = nnSearcherDest.countPointsStrictlyWithinR(
							t, kthNnData.distance, dynCorrExclTime);

			sum2xkNNDist += kthNnData.distance + kthNnData.distance;
			sumNx += n_x;
			sumNy += n_y;
			// And take the digammas:
			double digammaNxPlusOne = MathsUtils.digamma(n_x+1);
			double digammaNyPlusOne = MathsUtils.digamma(n_y+1);
			sumDiGammas += digammaNxPlusOne + digammaNyPlusOne;

			if (returnLocals) {
				// TODO should digammaN be adjusted if we are using 
				//  dynamic correlation exclusion, since we're not sampling from
				//  N after all? Think about this for all Kraskov estimators.
				localMi[t-startTimePoint] = digammaK - digammaNxPlusOne - digammaNyPlusOne + digammaN;
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
			return localMi;
		} else {
			return new double[] {sumDiGammas, sumNx, sumNy, sum2xkNNDist};
		}
	}

	/**
	 * Mirrors {@link #partialComputePredictionErrorFromObservations(int, int, int, boolean)}
	 * in implementing the guts of each Kraskov algorithm;
	 * however this is not intended to be used in a computation per se but for debugging
	 * purposes where the caller specifically wants to examine the neighbour counts
	 * for each data point (which is what is returned)
	 * 
	 * @param startTimePoint
	 * @param numTimePoints
	 * @return an array of arrays of neighbour counts for each sample (first index), where
	 *  the array of neighbour counts is for variable 1 or x (index 0) then variable 2 or y (index 1)
	 * @throws Exception
	 */
	public int[][] partialNeighbourCountFromObservations(
			int startTimePoint, int numTimePoints) throws Exception {
		
		double startTime = Calendar.getInstance().getTimeInMillis();
		
		ensureKdTreesConstructed();

		int[][] neighbourCounts = new int[numTimePoints][2];

		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps for this time step by
			//  finding the kth closest neighbour for point t:
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreeJoint.findKNearestNeighbours(k, t, dynCorrExclTime);
			// First element in the PQ is the kth NN,
			//  and epsilon = kthNnData.distance
			NeighbourNodeData kthNnData = nnPQ.poll();
			
			// Count the number of points whose x distance is less
			//  than eps, and whose y distance is less than
			//  epsilon = kthNnData.distance
			int n_x = nnSearcherSource.countPointsStrictlyWithinR(
							t, kthNnData.distance, dynCorrExclTime);
			int n_y = nnSearcherDest.countPointsStrictlyWithinR(
							t, kthNnData.distance, dynCorrExclTime);
			
			neighbourCounts[t - startTimePoint][0] = n_x;
			neighbourCounts[t - startTimePoint][1] = n_y;

		}
		
		if (debug) {
			Calendar rightNow2 = Calendar.getInstance();
			long endTime = rightNow2.getTimeInMillis();
			System.out.println("Subset " + startTimePoint + ":" +
					(startTimePoint + numTimePoints) + " Calculation time: " +
					((endTime - startTime)/1000.0) + " sec" );
		}

		return neighbourCounts;
		
	}

	@Override
	protected double[] partialComputeFromNewObservations(int startTimePoint, int numTimePoints,
			double[][] newVar1Observations, double[][] newVar2Observations, boolean returnLocals) throws Exception {

		double startTime = Calendar.getInstance().getTimeInMillis();

		double[] localMi = null;
		if (returnLocals) {
			localMi = new double[numTimePoints];
		}
		
		// Count the average number of points within eps_x and eps_y of each point
		double sumDiGammas = 0;
		double sumNx = 0;
		double sumNy = 0;
		
		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps for this time step by
			//  finding the kth closest neighbour for the new sample:
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreeJoint.findKNearestNeighbours(k,
							new double[][] {newVar1Observations[t], newVar2Observations[t]});
			// First element in the PQ is the kth NN,
			//  and epsilon = kthNnData.distance
			NeighbourNodeData kthNnData = nnPQ.poll();

			// Count the number of points whose x distance is less
			//  than eps, and whose y distance is less than
			//  epsilon = kthNnData.distance
			int n_x = nnSearcherSource.countPointsWithinR(
							new double[][] {newVar1Observations[t]},
							kthNnData.distance, false);
			int n_y = nnSearcherDest.countPointsWithinR(
						new double[][] {newVar2Observations[t]},
							kthNnData.distance, false);

			sumNx += n_x;
			sumNy += n_y;
			// And take the digammas:
			double digammaNxPlusOne = MathsUtils.digamma(n_x+1);
			double digammaNyPlusOne = MathsUtils.digamma(n_y+1);
			sumDiGammas += digammaNxPlusOne + digammaNyPlusOne;

			if (returnLocals) {
				// For new observations we're taking the probability counts over an extra point (no self-exclusion)
				//  so we don't use digamma(N) but digamma(N+1)
				localMi[t-startTimePoint] = digammaK - digammaNxPlusOne - digammaNyPlusOne + MathsUtils.digamma(totalObservations+1);
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
			return localMi;
		} else {
			return new double[] {sumDiGammas, sumNx, sumNy};
		}
	}	

}
