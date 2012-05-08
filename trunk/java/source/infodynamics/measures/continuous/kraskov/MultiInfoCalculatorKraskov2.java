package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Compute the Multi Info using the Kraskov estimation method.
 * Uses the second algorithm (defined on p.5 of the paper)</p>
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * @see http://dx.doi.org/10.1103/PhysRevE.69.066138
 * 
 * TODO should use the kthMinIndices code from MatrixUtils as per MutualInfoCalculatorMultiVariateKraskov2
 * 
 * @author Joseph Lizier
 *
 */
public class MultiInfoCalculatorKraskov2
	extends MultiInfoCalculatorKraskov {

	protected static final int JOINT_NORM_VAL_COLUMN = 0;
	protected static final int JOINT_NORM_TIMESTEP_COLUMN = 1;

	/**
	 * Compute the average MI from the previously set observations
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		if (miComputed) {
			return mi;
		}
		return computeAverageLocalOfObservations(null);
	}

	/**
	 * Compute what the average MI would look like were the 2nd, 3rd, etc time series reordered
	 *  as per the array of time indices in reordering.
	 * The user should ensure that all values 0..N-1 are represented exactly once in the
	 *  array reordering and that no other values are included here. 
	 * 
	 * @param reordering the specific new orderings to use. First index is the variable number
	 *  (can be for all variables, or one less than all if the first is not to be reordered),
	 *  second index is the time step, the value is the reordered time step to use
	 *  for that variable at the given time step.
	 *  If null, no reordering is performed.
	 * @return
	 * @throws Exception
	 */
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
	 * This method correctly computes the average local MI, but recomputes the x and y 
	 *  distances between all tuples in time.
	 * Kept here for cases where we have too many observations
	 *  to keep the norm between all pairs, and for testing purposes.
	 * 
	 * @return
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

	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e - 1/k=%.3e + digamma(N=%d)=%.3e => %.3e",
				k, MathsUtils.digamma(k), 1.0/(double)k, N, MathsUtils.digamma(N),
				(MathsUtils.digamma(k) - 1.0/(double)k + MathsUtils.digamma(N)));
		return constants;
	}
}
