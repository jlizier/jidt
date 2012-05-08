package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Compute the Conditional Mutual Info using the Kraskov estimation method,
 * as extended by Frenzel and Pompe.</p>
 * <p>
 * Uses the second algorithm (defined at start of p.3 of the Kraskov paper) - 
 *  note that Frenzel and Pompe only extended the technique directly for the
 *  first algorithm, though here we define it for the second.
 *  It is unclear exactly how to do this, since we need to account for
 *  the 1/k factor, which changes in the space of the second
 *  MI. I've taken a guess, though use of 
 *  {@link ConditionalMutualInfoCalculatorMultiVariateKraskov1}
 *  is perhaps recommended instead of this class.</p>
 *  
 * <p>Computes this directly looking at the marginal space for each variable, rather than
 * using the multi-info (or integration) in the marginal spaces.
 * </p>
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * 		http://dx.doi.org/10.1103/PhysRevE.69.066138
 * @see "Partial Mutual Information for Coupling Analysis of Multivariate Time Series", Frenzel and Pompe, 2007
 * 
 * @author Joseph Lizier
 *
 */
public class ConditionalMutualInfoCalculatorMultiVariateKraskov2
	extends ConditionalMutualInfoCalculatorMultiVariateKraskov {

	protected static final int JOINT_NORM_VAL_COLUMN = 0;
	protected static final int JOINT_NORM_TIMESTEP_COLUMN = 1;

	// Multiplier used in hueristic for determining whether to use a linear search
	//  for min kth element or a binary search.
	protected static final double CUTOFF_MULTIPLIER = 1.5;
	
	/**
	 * Compute what the average conditional MI would look like were the second time series reordered
	 *  as per the array of time indices in reordering.
	 * The user should ensure that all values 0..N-1 are represented exactly once in the
	 *  array reordering and that no other values are included here. 
	 * 
	 * @param reordering
	 * @return
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[] reordering) throws Exception {
		if (!tryKeepAllPairsNorms || (data1.length > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			double[][] originalData2 = data2;
			// Generate a new re-ordered data2
			data2 = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData2, reordering);
			// Compute the MI
			double newMI = computeAverageLocalOfObservationsWhileComputingDistances();
			// restore data2
			data2 = originalData2;
			return newMI;
		}
		
		// Otherwise we will use the norms we've already computed, and use a "virtual"
		//  reordered data2.
		
		if (xNorms == null) {
			computeNorms();
		}
		int N = data1.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));

		// Count the average number of points within eps_xz and eps_yz and eps_z
		double averageDiGammas = 0;
		double averageInverseCountInJointYZ = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_xz and eps_yz ad eps_z for this time step:
			//  First get x and y and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).

			int tForY = reordering[t];

			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				int t2ForY = reordering[t2];
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = Math.max(xNorms[t][t2],
							Math.max(yNorms[tForY][t2ForY], zNorms[t][t2]));
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_y = 0.0;
			double eps_z = 0.0;
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
			// Find eps_{x,y,z} as the maximum x and y and z norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xNorms[t][timeStepOfJthPoint] > eps_x) {
					eps_x = xNorms[t][timeStepOfJthPoint];
				}
				if (yNorms[tForY][reordering[timeStepOfJthPoint]] > eps_y) {
					eps_y = yNorms[tForY][reordering[timeStepOfJthPoint]];
				}
				if (zNorms[t][timeStepOfJthPoint] > eps_z) {
					eps_z = zNorms[t][timeStepOfJthPoint];
				}
			}
			
			// Count the number of points whose x distance is less
			//  than or equal to eps_x, and whose y distance is less
			//  than or equal to eps_y, and whose z distance is less
			//  than or equal to eps_z
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (zNorms[t][t2] <= eps_z) {
					n_z++;
					if (xNorms[t][t2] <= eps_x) {
						n_xz++;
					}
					if (yNorms[tForY][reordering[t2]] <= eps_y) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_z) - MathsUtils.digamma(n_xz)
								- MathsUtils.digamma(n_yz);
			averageInverseCountInJointYZ += 1.0 / (double) n_yz;
		}
		averageDiGammas /= (double) N;
		averageInverseCountInJointYZ /= (double) N;
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.println(String.format("Average n_xz=%.3f, Average n_yz=%.3f, Average n_z=%.3f",
						avNxz, avNyz, avNz));
		}
		
		condMi = MathsUtils.digamma(k) - 1.0 / (double) k + 
			+ averageInverseCountInJointYZ + averageDiGammas;
		condMiComputed = true;
		return condMi;
	}

	public double computeAverageLocalOfObservations() throws Exception {
		if (!tryKeepAllPairsNorms || (data1.length > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			return computeAverageLocalOfObservationsWhileComputingDistances();
		}
		
		if (xNorms == null) {
			computeNorms();
		}
		int N = data1.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double averageInverseCountInJointYZ = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_y and eps_z for this time step:
			//  using x and y and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = Math.max(xNorms[t][t2],
							Math.max(yNorms[t][t2], zNorms[t][t2]));
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_y = 0.0;
			double eps_z = 0.0;
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
			// Find eps_{x,y,z} as the maximum x and y and z norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xNorms[t][timeStepOfJthPoint] > eps_x) {
					eps_x = xNorms[t][timeStepOfJthPoint];
				}
				if (yNorms[t][timeStepOfJthPoint] > eps_y) {
					eps_y = yNorms[t][timeStepOfJthPoint];
				}
				if (zNorms[t][timeStepOfJthPoint] > eps_z) {
					eps_z = zNorms[t][timeStepOfJthPoint];
				}
			}
			
			// Count the number of points whose x distance is less
			//  than or equal to eps_x, and whose y distance is less
			//  than or equal to eps_y, and whose z distance is less
			//  than or equal to eps_z
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (zNorms[t][t2] <= eps_z) {
					n_z++;
					if (xNorms[t][t2] <= eps_x) {
						n_xz++;
					}
					if (yNorms[t][t2] <= eps_y) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_z) - MathsUtils.digamma(n_xz)
								- MathsUtils.digamma(n_yz);
			averageInverseCountInJointYZ += 1.0 / (double) n_yz;
		}
		averageDiGammas /= (double) N;
		averageInverseCountInJointYZ /= (double) N;
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f, Average n_y=%.3f",
					avNxz, avNyz, avNz));
		}
		
		condMi = MathsUtils.digamma(k) - 1.0 / (double) k +
			averageInverseCountInJointYZ + averageDiGammas;
		condMiComputed = true;
		return condMi;
	}

	/**
	 * This method correctly computes the average local MI, but recomputes the x and y and z
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

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double averageInverseCountInJointYZ = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_y and eps_z for this time step:
			//  First get x and y and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xyzNorms = EuclideanUtils.computeNorms(data1, data2, dataCond, t);
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = Math.max(xyzNorms[t2][0], 
						Math.max(xyzNorms[t2][1], xyzNorms[t2][2]));
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_y = 0.0;
			double eps_z = 0.0;
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
				if (xyzNorms[timeStepOfJthPoint][0] > eps_x) {
					eps_x = xyzNorms[timeStepOfJthPoint][0];
				}
				if (xyzNorms[timeStepOfJthPoint][1] > eps_y) {
					eps_y = xyzNorms[timeStepOfJthPoint][1];
				}
				if (xyzNorms[timeStepOfJthPoint][2] > eps_z) {
					eps_z = xyzNorms[timeStepOfJthPoint][2];
				}
			}

			// Count the number of points whose x distance is less
			//  than or equal to eps_x, and whose y distance is less
			//  than or equal to eps_y, and whose z distance is less
			//  than or equal to eps_z
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyzNorms[t2][2] <= eps_z) {
					n_z++;
					if (xyzNorms[t2][0] <= eps_x) {
						n_xz++;
					}
					if (xyzNorms[t2][1] <= eps_y) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_z) - MathsUtils.digamma(n_xz)
								- MathsUtils.digamma(n_yz);
			averageInverseCountInJointYZ += 1.0 / (double) n_yz;
		}
		averageDiGammas /= (double) N;
		averageInverseCountInJointYZ /= (double) N;
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.println(String.format("Average n_zx=%.3f, Average n_yz=%.3f, Average n_z=%.3f",
						avNxz, avNyz, avNz));
		}
		
		condMi = MathsUtils.digamma(k) - 1.0 / (double) k +
				averageInverseCountInJointYZ + averageDiGammas;
		condMiComputed = true;
		return condMi;
	}

	public double[] computeLocalOfPreviousObservations() throws Exception {
		int N = data1.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		double[] localCondMi = new double[N];
		
		// Constants:
		double digammaK = MathsUtils.digamma(k);
		
		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double averageInverseCountInJointYZ = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_y and eps_z for this time step:
			//  First get x and y and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xyzNorms = EuclideanUtils.computeNorms(data1, data2, dataCond, t);
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = Math.max(xyzNorms[t2][0],
							Math.max(xyzNorms[t2][1], xyzNorms[t2][2]));
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_y = 0.0;
			double eps_z = 0.0;
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
			// Find eps_{x,y,z} as the maximum x and y norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xyzNorms[timeStepOfJthPoint][0] > eps_x) {
					eps_x = xyzNorms[timeStepOfJthPoint][0];
				}
				if (xyzNorms[timeStepOfJthPoint][1] > eps_y) {
					eps_y = xyzNorms[timeStepOfJthPoint][1];
				}
				if (xyzNorms[timeStepOfJthPoint][0] > eps_z) {
					eps_z = xyzNorms[timeStepOfJthPoint][2];
				}
			}

			// Count the number of points whose x distance is less
			//  than or equal to eps_x, and whose y distance is less
			//  than or equal to eps_y, and whose z distance is less
			//  than or equal to eps_z
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyzNorms[t2][2] <= eps_z) {
					n_z++;
					if (xyzNorms[t2][0] <= eps_x) {
						n_xz++;
					}
					if (xyzNorms[t2][1] <= eps_y) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma:
			double digammaNxz = MathsUtils.digamma(n_xz);
			double digammaNyz = MathsUtils.digamma(n_yz);
			double digammaNz = MathsUtils.digamma(n_z);
			
			localCondMi[t] = digammaK - digammaNxz - digammaNyz + digammaNz
						- 1.0 / (double) k + 1.0/(double) n_yz ;

			averageDiGammas += digammaNz - digammaNxz - digammaNyz;
			averageInverseCountInJointYZ += 1.0 / (double) n_yz;
		}
		averageDiGammas /= (double) N;
		averageInverseCountInJointYZ /= (double) N;
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.println(String.format("Average n_xz=%.3f, Average n_yz=%.3f, Average n_z=%.3f",
					avNxz, avNyz, avNz));
		}
		
		condMi = digammaK + averageDiGammas - 1.0 / (double) k +
			averageInverseCountInJointYZ;
		condMiComputed = true;
		return localCondMi;
	}

	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e",
				k, MathsUtils.digamma(k));
		return constants;
	}
}
