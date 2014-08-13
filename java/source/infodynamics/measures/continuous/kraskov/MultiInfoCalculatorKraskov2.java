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

import infodynamics.measures.continuous.MultiInfoCalculator;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential multi-information of two given multivariate
 *  sets of
 *  observations (implementing {@link MultiInfoCalculator}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below),
 *  <b>algorithm 2</b>.
 *  Most of the functionality is defined by the parent class 
 *  {@link MultiInfoCalculatorKraskov}.</p>
 *
 * <p>Usage is as per the paradigm outlined for {@link MultiInfoCalculator},
 * and expanded on in {@link MultiInfoCalculatorKraskov}.
 * </p>
 *  
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * TODO should use the kthMinIndices code from MatrixUtils
 * as per MutualInfoCalculatorMultiVariateKraskov2
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MultiInfoCalculatorKraskov2
	extends MultiInfoCalculatorKraskov {

	protected static final int JOINT_NORM_VAL_COLUMN = 0;
	protected static final int JOINT_NORM_TIMESTEP_COLUMN = 1;

	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		if (miComputed) {
			return mi;
		}
		return computeAverageLocalOfObservations(null);
	}

	@Override
	public double computeAverageLocalOfObservations(int[][] reordering) throws Exception {
		if (V == 1) {
			miComputed = true;
			return 0.0;
		}
		if (!tryKeepAllPairsNorms || (data.length * V > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			double[][] originalData = data;
			// Generate a new re-ordered data
			if (reordering != null) {
				// Generate a new re-ordered data
				data = MatrixUtils.reorderDataForVariables(originalData, reordering);
			}
			// Compute the MI
			double newMI = computeAverageLocalOfObservationsWhileComputingDistances();
			// restore data
			data = originalData;
			return newMI;
		}
		
		// Otherwise we will use the norms we've already computed, and use a "virtual"
		//  reordered data2.
		
		if (norms == null) {
			computeNorms();
		}
		
		// Count the average number of points within eps_x[v]
		double averageDiGammas = 0;
		double[] avNx = new double[V];
		
		for (int t = 0; t < N; t++) {
			// Compute eps for each marginal for this time step:
			//  First grab marginal norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).

			// Get the reordered time steps for the variables
			int[] tForEachMarginal = reorderedTimeStepsForEachMarginal(reordering, t);

			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				int[] t2ForEachMarginal = reorderedTimeStepsForEachMarginal(reordering, t2);
				// Find the max marginal norm between the vector at t and the reordered vector
				//  at t2:
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = 0;
				for (int v = 0; v < V; v++) {
					double normForThisVar = norms[v][tForEachMarginal[v]][t2ForEachMarginal[v]];
					if (normForThisVar > jointNorm[t2][JOINT_NORM_VAL_COLUMN]) {
						jointNorm[t2][JOINT_NORM_VAL_COLUMN] = normForThisVar;
					}
				}
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours:
			double[] eps_x = new double[V];
			if (k == 1) {
				// just do a linear search for the minimum epsilon value
				int timeStepOfMin = MatrixUtils.minIndex(jointNorm, JOINT_NORM_VAL_COLUMN);
				int[] timeStepOfMinForEachMarginal =
					reorderedTimeStepsForEachMarginal(reordering, timeStepOfMin);
				for (int v = 0; v < V; v++) {
					eps_x[v] = norms[v][tForEachMarginal[v]][timeStepOfMinForEachMarginal[v]];
				}
			} else {
				// Sort the array of joint norms
				java.util.Arrays.sort(jointNorm, FirstIndexComparatorDouble.getInstance());
				// and now we have the closest k points.
				// Find eps_{x,y} as the maximum x and y norms amongst this set:
				for (int j = 0; j < k; j++) {
					int timeStepOfJthPoint = (int)jointNorm[j][JOINT_NORM_TIMESTEP_COLUMN];
					int[] timeStepOfJthPointForEachMarginal =
						reorderedTimeStepsForEachMarginal(reordering, timeStepOfJthPoint);
					for (int v = 0; v < V; v++) {
						if (norms[v][tForEachMarginal[v]][timeStepOfJthPointForEachMarginal[v]] > eps_x[v]) {
							eps_x[v] = norms[v][tForEachMarginal[v]][timeStepOfJthPointForEachMarginal[v]];
						}
					}
				}
			}
			
			// Count the number of points (in each marginal variable) 
			//  whose marginal distance is less than eps in that marginal dimension
			int[] n_x = new int[V];
			for (int t2 = 0; t2 < N; t2++) {
				int[] t2ForEachMarginal = reorderedTimeStepsForEachMarginal(reordering, t2);
				for (int v = 0; v < V; v++) {
					if (norms[v][tForEachMarginal[v]][t2ForEachMarginal[v]] <= eps_x[v]) {
						n_x[v]++;
					}
				}
			}
			// Track the averages, and take the digamma before adding into the 
			//  average:
			for (int v = 0; v < V; v++) {
				avNx[v] += n_x[v];
				averageDiGammas += MathsUtils.digamma(n_x[v]);
			}
		}
		averageDiGammas /= (double) N;
		if (debug) {
			for (int v = 0; v < V; v++) {
				avNx[v] /= (double)N;
				System.out.print(String.format("Average n_x[%d]=%.3f, ", v, avNx[v]));
			}
			System.out.println();
		}
		
		mi = MathsUtils.digamma(k) - (double) (V - 1) /(double)k - averageDiGammas +
				(double) (V - 1) * MathsUtils.digamma(N);
		miComputed = true;
		return mi;
	}

	/**
	 * This method correctly computes the average multi-info, but recomputes the
	 *  marginal distances between all tuples in time.
	 * Kept here for cases where we have too many observations
	 *  to keep the norm between all pairs, and for testing purposes.
	 * 
	 * @see #computeAverageLocalOfObservations()
	 * @return average multi-info
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservationsWhileComputingDistances() throws Exception {

		if (V == 1) {
			miComputed = true;
			return 0.0;
		}
		// Count the average number of points within eps for each marginal variable
		double averageDiGammas = 0;
		double[] avNx = new double[V];
	
		for (int t = 0; t < N; t++) {
			// Compute eps_x (for each marginal) for this time step:
			//  First get the marginal norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			
			double[][] normsForT = EuclideanUtils.computeNorms(data, t);
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = MatrixUtils.max(normsForT[t2]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			
			// Then find the k closest neighbours:
			double[] eps_x = new double[V];
			if (k == 1) {
				// just do a linear search for the minimum epsilon value
				int timeStepOfMin = MatrixUtils.minIndex(jointNorm, JOINT_NORM_VAL_COLUMN);
				for (int v = 0; v < V; v++) {
					eps_x[v] = normsForT[timeStepOfMin][v];
				}
			} else {
				// Sort the array of joint norms
				java.util.Arrays.sort(jointNorm, FirstIndexComparatorDouble.getInstance());
				// and now we have the closest k points.
				// Find eps_{x,y} as the maximum x and y norms amongst this set:
				for (int j = 0; j < k; j++) {
					int timeStepOfJthPoint = (int)jointNorm[j][JOINT_NORM_TIMESTEP_COLUMN];
					for (int v = 0; v < V; v++) {
						if (normsForT[timeStepOfJthPoint][v] > eps_x[v]) {
							eps_x[v] = normsForT[timeStepOfJthPoint][v];
						}
					}
				}
			}
			
			// Count the number of points (in each marginal variable) 
			//  whose marginal distance is less than eps in that marginal dimension
			int[] n_x = new int[V];
			for (int t2 = 0; t2 < N; t2++) {
				for (int v = 0; v < V; v++) {
					if (normsForT[t2][v] <= eps_x[v]) {
						n_x[v]++;
					}
				}
			}
			// Track the averages, and take the digamma before adding into the 
			//  average:
			for (int v = 0; v < V; v++) {
				avNx[v] += n_x[v];
				averageDiGammas += MathsUtils.digamma(n_x[v]);
			}
		}
		averageDiGammas /= (double) N;
		if (debug) {
			for (int v = 0; v < V; v++) {
				avNx[v] /= (double)N;
				System.out.print(String.format("Average n_x[%d]=%.3f, ", v, avNx[v]));
			}
			System.out.println();
		}
		
		mi = MathsUtils.digamma(k) - (double) (V - 1) /(double)k - averageDiGammas +
				(double) (V - 1) * MathsUtils.digamma(N);
		miComputed = true;
		return mi;
	}

	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] localMi = new double[N];
		if (V == 1) {
			miComputed = true;
			return localMi;
		}

		// Constants:
		double digammaK = MathsUtils.digamma(k);
		double Vminus1TimesDigammaN = (double) (V - 1) * MathsUtils.digamma(N);
		double Vminus1TimesInvK = (double) (V - 1) / (double)k;
		
		// Count the average number of points within eps_x[v] for each marginal v
		double averageDiGammas = 0;
		double[] avNx = new double[V];
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x (for each marginal) for this time step:
			//  First get the marginal norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			
			double[][] normsForT = EuclideanUtils.computeNorms(data, t);
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = MatrixUtils.max(normsForT[t2]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}

			// Then find the k closest neighbours:
			double[] eps_x = new double[V];
			if (k == 1) {
				// just do a linear search for the minimum epsilon value
				int timeStepOfMin = MatrixUtils.minIndex(jointNorm, JOINT_NORM_VAL_COLUMN);
				for (int v = 0; v < V; v++) {
					eps_x[v] = normsForT[timeStepOfMin][v];
				}
			} else {
				// Sort the array of joint norms
				java.util.Arrays.sort(jointNorm, FirstIndexComparatorDouble.getInstance());
				// and now we have the closest k points.
				// Find eps_{x,y} as the maximum x and y norms amongst this set:
				for (int j = 0; j < k; j++) {
					int timeStepOfJthPoint = (int)jointNorm[j][JOINT_NORM_TIMESTEP_COLUMN];
					for (int v = 0; v < V; v++) {
						if (normsForT[timeStepOfJthPoint][v] > eps_x[v]) {
							eps_x[v] = normsForT[timeStepOfJthPoint][v];
						}
					}
				}
			}
			
			// Count the number of points (in each marginal variable) 
			//  whose marginal distance is less than eps in that marginal dimension
			int[] n_x = new int[V];
			for (int t2 = 0; t2 < N; t2++) {
				for (int v = 0; v < V; v++) {
					if (normsForT[t2][v] <= eps_x[v]) {
						n_x[v]++;
					}
				}
			}

			// Track the averages, and take the digamma before adding into the 
			//  local:
			localMi[t] = digammaK - Vminus1TimesInvK + Vminus1TimesDigammaN;
			for (int v = 0; v < V; v++) {
				double digammaNx = MathsUtils.digamma(n_x[v]);
				localMi[t] -= digammaNx;
				// And keep track of the averages
				avNx[v] += n_x[v];
				averageDiGammas += digammaNx;
			}
		}
		if (debug) {
			for (int v = 0; v < V; v++) {
				avNx[v] /= (double)N;
				System.out.print(String.format("Average n_x[%d]=%.3f, ", v, avNx[v]));
			}
			System.out.println();
		}
		
		mi = digammaK - Vminus1TimesInvK - averageDiGammas + Vminus1TimesDigammaN;
		miComputed = true;
		return localMi;
	}

	@Override
	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e - 1/k=%.3e + digamma(N=%d)=%.3e => %.3e",
				k, MathsUtils.digamma(k), 1.0/(double)k, N, MathsUtils.digamma(N),
				(MathsUtils.digamma(k) - 1.0/(double)k + MathsUtils.digamma(N)));
		return constants;
	}
}
