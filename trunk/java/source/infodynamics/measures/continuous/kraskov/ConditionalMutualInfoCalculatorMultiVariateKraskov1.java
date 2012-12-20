package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Compute the Conditional Mutual Info using the Kraskov estimation method,
 * as extended by Frenzel and Pompe.
 * Uses the first algorithm (defined at end of p.2 of the Kraskov paper)</p>
 * <p>Computes this directly looking at the marginal space for each variable, rather than
 * using the multi-info (or integration) in the marginal spaces.
 * </p>
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * http://dx.doi.org/10.1103/PhysRevE.69.066138
 * @see "Partial Mutual Information for Coupling Analysis of Multivariate Time Series", Frenzel and Pompe, 2007
 * 
 * @author Joseph Lizier
 *
 */
public class ConditionalMutualInfoCalculatorMultiVariateKraskov1
	extends ConditionalMutualInfoCalculatorMultiVariateKraskov {
	
	// Multiplier used in hueristic for determining whether to use a linear search
	//  for min kth element or a binary search.
	protected static final double CUTOFF_MULTIPLIER = 1.5;

	/**
	 * Compute the average conditional MI from the previously set observations
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		return computeAverageLocalOfObservations(null);
	}

	/**
	 * Compute what the average conditional MI would look like were the second time series reordered
	 *  as per the array of time indices in reordering.
	 * The user should ensure that all values 0..N-1 are represented exactly once in the
	 *  array reordering and that no other values are included here. 
	 * If reordering is null, it is assumed there is no reordering of
	 *  the y variable.
	 * 
	 * @param reordering the reordered time steps of the y variable
	 * @return
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[] reordering) throws Exception {
		if (!tryKeepAllPairsNorms || (data1.length > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			double[][] originalData2 = data2;
			if (reordering != null) {
				// Generate a new re-ordered data2
				data2 = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData2, reordering);
			}
			// Compute the MI
			double newMI = computeAverageLocalOfObservationsWhileComputingDistances();
			// restore data2
			data2 = originalData2;
			return newMI;
		}
		
		// Else we'll use the arrays of marginal distances
		
		if (xNorms == null) {
			computeNorms();
		}
		int N = data1.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		
		// Count the average number of points within eps_xz and eps_yz and eps_z
		double averageDiGammas = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  using x, y and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			
			int tForY = (reordering == null) ? t : reordering[t];

			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				int t2ForY = (reordering == null) ? t2 : reordering[t2];
				// Joint norm is the max of all three marginals
				jointNorm[t2] = Math.max(xNorms[t][t2], Math.max(yNorms[tForY][t2ForY], zNorms[t][t2]));
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
			
			// Count the number of points whose (x,z) distance is less
			//  than eps, and whose (y,z) distance is less than eps, and
			//  whose z distance is less than eps.
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (zNorms[t][t2] < epsilon) {
					n_z++;
					if (xNorms[t][t2] < epsilon) {
						n_xz++;
					}
					int t2ForY = (reordering == null) ? t2 : reordering[t2];
					if (yNorms[tForY][t2ForY] < epsilon) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			// Note: we're using digamma function which has opposite sign to the harmonic
			//  number used by Frenzel and Pompe, and is also offset by a constant (though
			//  this cancels out)
			averageDiGammas += MathsUtils.digamma(n_z+1) - MathsUtils.digamma(n_xz+1)
							- MathsUtils.digamma(n_yz+1);
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.println(String.format("Average n_xz=%.3f, Average n_yz=%.3f, Average n_z=%.3f", avNxz, avNyz, avNz));
		}
		
		condMi = MathsUtils.digamma(k) + averageDiGammas;
		condMiComputed = true;
		return condMi;
	}
	
	/**
	 * This method correctly computes the average local MI, but recomputes the x, y and z 
	 *  distances between all tuples in time.
	 * Kept here for cases where we have too many observations
	 *  to keep the norm between all pairs, and for testing purposes.
	 * 
	 * @return
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservationsWhileComputingDistances() throws Exception {
		int N = data1.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		
		// Count the average number of points within eps_xz, eps_yz and eps_z
		double averageDiGammas = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  First get x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xyzNorms = normCalculator.computeNorms(data1, data2, dataCond, t);
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
			
			// Count the number of points whose (x,z) distance is less
			//  than eps, whose (y,z) distance is less than eps, and whose
			//  z distance is less than eps
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
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_z+1) - MathsUtils.digamma(n_xz+1)
						- MathsUtils.digamma(n_yz+1);
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			System.out.println(String.format("Average n_xz=%.3f, Average n_yz=%.3f, Average n_z=%.3f",
					avNxz, avNyz, avNz));
		}
		
		condMi = MathsUtils.digamma(k) + averageDiGammas;
		condMiComputed = true;
		return condMi;
	}

	public double[] computeLocalOfPreviousObservations() throws Exception {
		int N = data1.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		double[] localCondMi = new double[N];
		
		// Constants:
		double digammaK = MathsUtils.digamma(k);
		
		// Count the average number of points within eps_xz and eps_yz and eps_z
		double averageDiGammas = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  First get x and y and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity.
			double[][] xyzNorms = normCalculator.computeNorms(data1, data2, dataCond, t);
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
			
			// Count the number of points whose x distance is less
			//  than eps, and whose y distance is less than eps
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
			// And take the digamma:
			double digammaNxzPlusOne = MathsUtils.digamma(n_xz+1);
			double digammaNyzPlusOne = MathsUtils.digamma(n_yz+1);
			double digammaNzPlusOne = MathsUtils.digamma(n_z+1);
			
			localCondMi[t] = digammaK - digammaNxzPlusOne - digammaNyzPlusOne + digammaNzPlusOne;
			
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And keep track of the average
			averageDiGammas += digammaNzPlusOne - digammaNxzPlusOne - digammaNyzPlusOne;
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.println(String.format("Average n_xz=%.3f, Average n_yz=%.3f, Average n_z=%.3f",
					avNxz, avNyz, avNz));
		}
		
		condMi = digammaK + averageDiGammas;
		condMiComputed = true;
		return localCondMi;
	}

	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e",
				k, MathsUtils.digamma(k));
		return constants;
	}
}
