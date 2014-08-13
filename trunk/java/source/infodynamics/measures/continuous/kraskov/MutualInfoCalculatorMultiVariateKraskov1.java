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

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

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
 */
public class MutualInfoCalculatorMultiVariateKraskov1
	extends MutualInfoCalculatorMultiVariateKraskov {
	
	/**
	 * Multiplier used as hueristic for determining whether to use a linear search
	 *  for kth nearest neighbour or a binary search.
	 */
	protected static final double CUTOFF_MULTIPLIER = 1.5;

	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		return computeAverageLocalOfObservations(null);
	}

	@Override
	public double computeAverageLocalOfObservations(int[] reordering) throws Exception {
		if (!tryKeepAllPairsNorms || (sourceObservations.length > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			double[][] originalData2 = destObservations;
			if (reordering != null) {
				// Generate a new re-ordered data2
				destObservations = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData2, reordering);
			}
			// Compute the MI
			double newMI = computeAverageLocalOfObservationsWhileComputingDistances();
			// restore data2
			destObservations = originalData2;
			return newMI;
		}
		
		if (xNorms == null) {
			computeNorms();
		}
		int N = sourceObservations.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		
		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double avNx = 0;
		double avNy = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  using x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			
			int tForY = (reordering == null) ? t : reordering[t];

			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				int t2ForY = (reordering == null) ? t2 : reordering[t2];
				jointNorm[t2] = Math.max(xNorms[t][t2], yNorms[tForY][t2ForY]);
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
			
			// Count the number of points whose x distance is less
			//  than eps, and whose y distance is less than eps
			int n_x = 0;
			int n_y = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xNorms[t][t2] < epsilon) {
					n_x++;
				}
				int t2ForY = (reordering == null) ? t2 : reordering[t2];
				if (yNorms[tForY][t2ForY] < epsilon) {
					n_y++;
				}
			}
			avNx += n_x;
			avNy += n_y;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_x+1) + MathsUtils.digamma(n_y+1);
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
		}
		
		double average = MathsUtils.digamma(k) - averageDiGammas + MathsUtils.digamma(N);
		miComputed = true;
		if (reordering == null) {
			lastAverage = average;
		}
		return average;
	}
	
	/**
	 * This method correctly computes the average MI, but recomputes the x and y 
	 *  distances between all tuples in time.
	 * Kept here for cases where we have too many observations
	 *  to keep the norm between all pairs, and for testing purposes.
	 * 
	 * @see #computeAverageLocalOfObservations()
	 * @return average MI value in nats not bits
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservationsWhileComputingDistances() throws Exception {
		int N = sourceObservations.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		
		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double avNx = 0;
		double avNy = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  First get x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xyNorms = normCalculator.computeNorms(sourceObservations, destObservations, t);
			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2] = Math.max(xyNorms[t2][0], xyNorms[t2][1]);
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
			
			// Count the number of points whose x distance is less
			//  than eps, and whose y distance is less than eps
			int n_x = 0;
			int n_y = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyNorms[t2][0] < epsilon) {
					n_x++;
				}
				if (xyNorms[t2][1] < epsilon) {
					n_y++;
				}
			}
			avNx += n_x;
			avNy += n_y;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_x+1) + MathsUtils.digamma(n_y+1);
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
		}
		
		lastAverage = MathsUtils.digamma(k) - averageDiGammas + MathsUtils.digamma(N);
		miComputed = true;
		return lastAverage;
	}

	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		int N = sourceObservations.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		double[] localMi = new double[N];
		
		// Constants:
		double digammaK = MathsUtils.digamma(k);
		double digammaN = MathsUtils.digamma(N);
		
		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double avNx = 0;
		double avNy = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  First get x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity.
			double[][] xyNorms = normCalculator.computeNorms(sourceObservations, destObservations, t);
			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2] = Math.max(xyNorms[t2][0], xyNorms[t2][1]);
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
			
			// Count the number of points whose x distance is less
			//  than eps, and whose y distance is less than eps
			int n_x = 0;
			int n_y = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyNorms[t2][0] < epsilon) {
					n_x++;
				}
				if (xyNorms[t2][1] < epsilon) {
					n_y++;
				}
			}
			// And take the digamma:
			double digammaNxPlusOne = MathsUtils.digamma(n_x+1);
			double digammaNyPlusOne = MathsUtils.digamma(n_y+1);
			
			localMi[t] = digammaK - digammaNxPlusOne - digammaNyPlusOne + digammaN;
			
			avNx += n_x;
			avNy += n_y;
			// And keep track of the average
			averageDiGammas += digammaNxPlusOne + digammaNyPlusOne;
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
		}
		
		lastAverage = digammaK - averageDiGammas + digammaN;
		miComputed = true;
		return localMi;
	}

	@Override
	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e + digamma(N=%d)=%.3e => %.3e",
				k, MathsUtils.digamma(k), N, MathsUtils.digamma(N),
				(MathsUtils.digamma(k) + MathsUtils.digamma(N)));
		return constants;
	}
}
