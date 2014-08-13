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
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential mutual information of two given multivariate sets of
 *  observations (implementing {@link MutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below),
 *  <b>algorithm 2</b>.
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
public class MutualInfoCalculatorMultiVariateKraskov2
	extends MutualInfoCalculatorMultiVariateKraskov {

	protected static final int JOINT_NORM_VAL_COLUMN = 0;
	protected static final int JOINT_NORM_TIMESTEP_COLUMN = 1;

	/**
	 * Multiplier used as hueristic for determining whether to use a linear search
	 *  for kth nearest neighbour or a binary search.
	 */
	protected static final double CUTOFF_MULTIPLIER = 1.5;
	
	public double computeAverageLocalOfObservations(int[] reordering) throws Exception {
		if (!tryKeepAllPairsNorms || (sourceObservations.length > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			double[][] originalData2 = destObservations;
			// Generate a new re-ordered data2
			destObservations = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData2, reordering);
			// Compute the MI
			double newMI = computeAverageLocalOfObservationsWhileComputingDistances();
			// restore data2
			destObservations = originalData2;
			return newMI;
		}
		
		// Otherwise we will use the norms we've already computed, and use a "virtual"
		//  reordered data2.
		
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
			// Compute eps_x and eps_y for this time step:
			//  First get x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).

			int tForY = reordering[t];

			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				int t2ForY = reordering[t2];
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = Math.max(xNorms[t][t2], yNorms[tForY][t2ForY]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_y = 0.0;
			int[] timeStepsOfKthMins = null;
			if (k <= cutoffForKthMinLinear) {
				// just do a linear search for the minimum epsilon value
				timeStepsOfKthMins = MatrixUtils.kMinIndices(jointNorm, JOINT_NORM_VAL_COLUMN, k);
			} else {
				// Sort the array of joint norms
				java.util.Arrays.sort(jointNorm, FirstIndexComparatorDouble.getInstance());
				// and now we have the closest k points.
				timeStepsOfKthMins = new int[k];
				for (int j = 0; j < k; j++) {
					timeStepsOfKthMins[j] = (int) jointNorm[j][JOINT_NORM_TIMESTEP_COLUMN];
				}
			}
			// and now we have the closest k points.
			// Find eps_{x,y} as the maximum x and y norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xNorms[t][timeStepOfJthPoint] > eps_x) {
					eps_x = xNorms[t][timeStepOfJthPoint];
				}
				if (yNorms[tForY][reordering[timeStepOfJthPoint]] > eps_y) {
					eps_y = yNorms[tForY][reordering[timeStepOfJthPoint]];
				}
			}
			
			// Count the number of points whose x distance is less
			//  than or equal to eps_x, and whose y distance is less
			//  than or equal to eps_y
			int n_x = 0;
			int n_y = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xNorms[t][t2] <= eps_x) {
					n_x++;
				}
				if (yNorms[tForY][reordering[t2]] <= eps_y) {
					n_y++;
				}
			}
			avNx += n_x;
			avNy += n_y;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_x) + MathsUtils.digamma(n_y);
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
		}
		
		double average = MathsUtils.digamma(k) - 1.0/(double)k - averageDiGammas + MathsUtils.digamma(N);
		miComputed = true;
		if (reordering == null) {
			lastAverage = average;
		}
		return average;
	}

	public double computeAverageLocalOfObservations() throws Exception {
		if (!tryKeepAllPairsNorms || (sourceObservations.length > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			return computeAverageLocalOfObservationsWhileComputingDistances();
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
			// Compute eps_x and eps_y for this time step:
			//  using x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = Math.max(xNorms[t][t2], yNorms[t][t2]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_y = 0.0;
			int[] timeStepsOfKthMins = null;
			if (k <= cutoffForKthMinLinear) {
				// just do a linear search for the minimum epsilon value
				timeStepsOfKthMins = MatrixUtils.kMinIndices(jointNorm, JOINT_NORM_VAL_COLUMN, k);
			} else {
				// Sort the array of joint norms
				java.util.Arrays.sort(jointNorm, FirstIndexComparatorDouble.getInstance());
				// and now we have the closest k points.
				timeStepsOfKthMins = new int[k];
				for (int j = 0; j < k; j++) {
					timeStepsOfKthMins[j] = (int) jointNorm[j][JOINT_NORM_TIMESTEP_COLUMN];
				}
			}
			// and now we have the closest k points.
			// Find eps_{x,y} as the maximum x and y norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xNorms[t][timeStepOfJthPoint] > eps_x) {
					eps_x = xNorms[t][timeStepOfJthPoint];
				}
				if (yNorms[t][timeStepOfJthPoint] > eps_y) {
					eps_y = yNorms[t][timeStepOfJthPoint];
				}
			}
			
			// Count the number of points whose x distance is less
			//  than or equal to eps_x, and whose y distance is less
			//  than or equal to eps_y
			int n_x = 0;
			int n_y = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xNorms[t][t2] <= eps_x) {
					n_x++;
				}
				if (yNorms[t][t2] <= eps_y) {
					n_y++;
				}
			}
			avNx += n_x;
			avNy += n_y;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_x) + MathsUtils.digamma(n_y);
			// if (debug) {
			// 	System.out.printf("n=%d, \n");
			// }
		}
		averageDiGammas /= (double) N;
		lastAverage = MathsUtils.digamma(k) - 1.0/(double)k - averageDiGammas + MathsUtils.digamma(N);
		miComputed = true;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.printf("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy);
			System.out.printf("psi(k=%d)=%.4f - 1/k=%.4f - averageDiGammas=%.4f -psi(N)=%.4f => %.4f\n",
					k, MathsUtils.digamma(k), 1.0/(double)k, averageDiGammas, MathsUtils.digamma(N), lastAverage);
		}
		
		return lastAverage;
	}

	/**
	 * This method correctly computes the average local MI, but recomputes the x and y 
	 *  distances between all tuples in time.
	 * Kept here for cases where we have too many observations
	 *  to keep the norm between all pairs, and for testing purposes.
	 * 
	 * @see #computeAverageLocalOfObservations()
	 * @return
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
			// Compute eps_x and eps_y for this time step:
			//  First get x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xyNorms = normCalculator.computeNorms(sourceObservations, destObservations, t);
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = Math.max(xyNorms[t2][0], xyNorms[t2][1]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_y = 0.0;
			int[] timeStepsOfKthMins = null;
			if (k <= cutoffForKthMinLinear) {
				// just do a linear search for the minimum epsilon value
				timeStepsOfKthMins = MatrixUtils.kMinIndices(jointNorm, JOINT_NORM_VAL_COLUMN, k);
			} else {
				// Sort the array of joint norms
				java.util.Arrays.sort(jointNorm, FirstIndexComparatorDouble.getInstance());
				// and now we have the closest k points.
				timeStepsOfKthMins = new int[k];
				for (int j = 0; j < k; j++) {
					timeStepsOfKthMins[j] = (int) jointNorm[j][JOINT_NORM_TIMESTEP_COLUMN];
				}
			}
			// and now we have the closest k points.
			// Find eps_{x,y} as the maximum x and y norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xyNorms[timeStepOfJthPoint][0] > eps_x) {
					eps_x = xyNorms[timeStepOfJthPoint][0];
				}
				if (xyNorms[timeStepOfJthPoint][1] > eps_y) {
					eps_y = xyNorms[timeStepOfJthPoint][1];
				}
			}

			// Count the number of points whose x distance is less
			//  than or equal to eps_x, and whose y distance is less
			//  than or equal to eps_y
			int n_x = 0;
			int n_y = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyNorms[t2][0] <= eps_x) {
					n_x++;
				}
				if (xyNorms[t2][1] <= eps_y) {
					n_y++;
				}
			}
			avNx += n_x;
			avNy += n_y;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_x) + MathsUtils.digamma(n_y);
		}
		averageDiGammas /= (double) N;
		lastAverage = MathsUtils.digamma(k) - 1.0/(double)k - averageDiGammas + MathsUtils.digamma(N);
		miComputed = true;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.printf("Average n_x=%.3f, Average n_y=%.3f\n", avNx, avNy);
			System.out.printf("psi(k=%d)=%.4f - 1/k=%.4f - averageDiGammas=%.4f + psi(N)=%.4f => %.4f\n",
					k, MathsUtils.digamma(k), 1.0/(double)k, averageDiGammas, MathsUtils.digamma(N), lastAverage);
		}
		
		return lastAverage;
	}

	public double[] computeLocalOfPreviousObservations() throws Exception {
		int N = sourceObservations.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		double[] localMi = new double[N];
		
		// Constants:
		double digammaK = MathsUtils.digamma(k);
		double invK = 1.0 / (double)k;
		double digammaN = MathsUtils.digamma(N);
		
		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double avNx = 0;
		double avNy = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_y for this time step:
			//  First get x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xyNorms = normCalculator.computeNorms(sourceObservations, destObservations, t);
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = Math.max(xyNorms[t2][0], xyNorms[t2][1]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_y = 0.0;
			int[] timeStepsOfKthMins = null;
			if (k <= cutoffForKthMinLinear) {
				// just do a linear search for the minimum epsilon value
				timeStepsOfKthMins = MatrixUtils.kMinIndices(jointNorm, JOINT_NORM_VAL_COLUMN, k);
			} else {
				// Sort the array of joint norms
				java.util.Arrays.sort(jointNorm, FirstIndexComparatorDouble.getInstance());
				// and now we have the closest k points.
				timeStepsOfKthMins = new int[k];
				for (int j = 0; j < k; j++) {
					timeStepsOfKthMins[j] = (int) jointNorm[j][JOINT_NORM_TIMESTEP_COLUMN];
				}
			}
			// and now we have the closest k points.
			// Find eps_{x,y} as the maximum x and y norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xyNorms[timeStepOfJthPoint][0] > eps_x) {
					eps_x = xyNorms[timeStepOfJthPoint][0];
				}
				if (xyNorms[timeStepOfJthPoint][1] > eps_y) {
					eps_y = xyNorms[timeStepOfJthPoint][1];
				}
			}

			// Count the number of points whose x distance is less
			//  than or equal to eps_x, and whose y distance is less
			//  than or equal to eps_y
			int n_x = 0;
			int n_y = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyNorms[t2][0] <= eps_x) {
					n_x++;
				}
				if (xyNorms[t2][1] <= eps_y) {
					n_y++;
				}
			}
			avNx += n_x;
			avNy += n_y;
			// And take the digamma:
			double digammaNx = MathsUtils.digamma(n_x);
			double digammaNy = MathsUtils.digamma(n_y);
			
			localMi[t] = digammaK - invK - digammaNx - digammaNy + digammaN;

			averageDiGammas += digammaNx + digammaNy;
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
		}
		
		lastAverage = digammaK - invK - averageDiGammas + digammaN;
		miComputed = true;
		return localMi;
	}

	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e - 1/k=%.3e + digamma(N=%d)=%.3e => %.3e",
				k, MathsUtils.digamma(k), 1.0/(double)k, N, MathsUtils.digamma(N),
				(MathsUtils.digamma(k) - 1.0/(double)k + MathsUtils.digamma(N)));
		return constants;
	}
}
