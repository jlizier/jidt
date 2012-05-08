package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Compute the Multi Info using the Kraskov estimation method.
 * Uses the first algorithm (defined on p.4/5 of the paper)</p>
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * @see http://dx.doi.org/10.1103/PhysRevE.69.066138
 * 
 * @author Joseph Lizier
 *
 */
public class MultiInfoCalculatorKraskov1
	extends MultiInfoCalculatorKraskov {
	
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
	 * If reordering is null, it is assumed there is no reordering of
	 *  the y variable.
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
		
		if (norms == null) {
			computeNorms();
		}
		
		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double[] avNx = new double[V];
		
		int cutoffForKthMinLinear = (int) (Math.log(N) / Math.log(2.0));

		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  First grab marginal norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).			
			
			// Create storage for the reordered time steps for the variables
			int[] tForEachMarginal = reorderedTimeStepsForEachMarginal(reordering, t);

			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				int[] t2ForEachMarginal = reorderedTimeStepsForEachMarginal(reordering, t2);
				// Find the max marginal norm between the vector at t and the reordered vector
				//  at t2:
				jointNorm[t2] = 0;
				for (int v = 0; v < V; v++) {
					double normForThisVar = norms[v][tForEachMarginal[v]][t2ForEachMarginal[v]];
					if (normForThisVar > jointNorm[t2]) {
						jointNorm[t2] = normForThisVar;
					}
				}
			}
			// Then find the kth closest neighbour:
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
			
			// Count the number of points (in each marginal variable) 
			//  whose marginal distance is less than eps
			int[] n_x = new int[V];
			for (int t2 = 0; t2 < N; t2++) {
				int[] t2ForEachMarginal = reorderedTimeStepsForEachMarginal(reordering, t2);
				for (int v = 0; v < V; v++) {
					if (norms[v][tForEachMarginal[v]][t2ForEachMarginal[v]] < epsilon) {
						n_x[v]++;
					}
				}
			}
			// Track the averages, and take the digamma before adding into the 
			//  average:
			for (int v = 0; v < V; v++) {
				avNx[v] += n_x[v];
				averageDiGammas += MathsUtils.digamma(n_x[v]+1);
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
		
		mi = MathsUtils.digamma(k) - averageDiGammas + (double) (V - 1) * MathsUtils.digamma(N);
		miComputed = true;
		return mi;
	}
	
	/**
	 * This method correctly computes the average local MI, but recomputes the 
	 *  distances for each marginal variable between all tuples in time.
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
		int cutoffForKthMinLinear = (int) (Math.log(N) / Math.log(2.0));

		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  First get marginal norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			
			double[][] normsForT = EuclideanUtils.computeNorms(data, t);
			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2] = MatrixUtils.max(normsForT[t2]);
			}
			// Then find the kth closest neighbour:
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
			
			// Count the number of points (in each marginal variable) 
			//  whose marginal distance is less than eps
			int[] n_x = new int[V];
			for (int t2 = 0; t2 < N; t2++) {
				for (int v = 0; v < V; v++) {
					if (normsForT[t2][v] < epsilon) {
						n_x[v]++;
					}
				}
			}
			for (int v = 0; v < V; v++) {
				avNx[v] += n_x[v];
			}
			// And take the digamma before adding into the 
			//  average:
			for (int v = 0; v < V; v++) {
				averageDiGammas += MathsUtils.digamma(n_x[v]+1);
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
		
		mi = MathsUtils.digamma(k) - averageDiGammas + (double) (V - 1) * MathsUtils.digamma(N);
		miComputed = true;
		return mi;
	}

	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] localMi = new double[N];
		int cutoffForKthMinLinear = (int) (Math.log(N) / Math.log(2.0));

		if (V == 1) {
			miComputed = true;
			return localMi;
		}

		// Constants:
		double digammaK = MathsUtils.digamma(k);
		double Vminus1TimesdigammaN = (double) (V - 1) * MathsUtils.digamma(N);
		
		// Count the average number of points within eps_x[v] for each marginal v
		double averageDiGammas = 0;
		double[] avNx = new double[V];
		
		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  First get marginal norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			
			double[][] norms = EuclideanUtils.computeNorms(data, t);
			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2] = MatrixUtils.max(norms[t2]);
			}
			
			// Then find the kth closest neighbour:
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
						
			// Count the number of points (in each marginal variable) 
			//  whose marginal distance is less than eps
			int[] n_x = new int[V];
			for (int t2 = 0; t2 < N; t2++) {
				for (int v = 0; v < V; v++) {
					if (norms[t2][v] < epsilon) {
						n_x[v]++;
					}
				}
			}

			// And take the digammas, and add into the local
			localMi[t] = digammaK + Vminus1TimesdigammaN;
			for (int v = 0; v < V; v++) {
				double digammaNxPlusOne = MathsUtils.digamma(n_x[v]+1);
				localMi[t] -= digammaNxPlusOne;
				// And keep track of the averages
				averageDiGammas += digammaNxPlusOne;
				avNx[v] += n_x[v];
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
		
		mi = digammaK - averageDiGammas + Vminus1TimesdigammaN;
		miComputed = true;
		return localMi;
	}

	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e + digamma(N=%d)=%.3e => %.3e",
				k, MathsUtils.digamma(k), N, MathsUtils.digamma(N),
				(MathsUtils.digamma(k) + MathsUtils.digamma(N)));
		return constants;
	}
}
