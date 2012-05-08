package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Compute the Mutual Info using the Kraskov estimation method.
 * Uses the second algorithm (defined at start of p.3 of the paper)</p>
 * <p>Computes this directly looking at the marginal space for each variable, rather than
 * using the multi-info (or integration) in the marginal spaces.
 * </p>
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * @see http://dx.doi.org/10.1103/PhysRevE.69.066138
 * 
 * @author Joseph Lizier
 *
 */
public class MutualInfoCalculatorMultiVariateKraskov2
	extends MutualInfoCalculatorMultiVariateKraskov {

	protected static final int JOINT_NORM_VAL_COLUMN = 0;
	protected static final int JOINT_NORM_TIMESTEP_COLUMN = 1;

	// Multiplier used in hueristic for determining whether to use a linear search
	//  for min kth element or a binary search.
	protected static final double CUTOFF_MULTIPLIER = 1.5;
	
	/**
	 * Compute what the average MI would look like were the second time series reordered
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
		
		mi = MathsUtils.digamma(k) - 1.0/(double)k - averageDiGammas + MathsUtils.digamma(N);
		miComputed = true;
		return mi;
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
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
		}
		
		mi = MathsUtils.digamma(k) - 1.0/(double)k - averageDiGammas + MathsUtils.digamma(N);
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
		int N = data1.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double avNx = 0;
		double avNy = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_y for this time step:
			//  First get x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xyNorms = EuclideanUtils.computeNorms(data1, data2, t);
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
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
		}
		
		mi = MathsUtils.digamma(k) - 1.0/(double)k - averageDiGammas + MathsUtils.digamma(N);
		miComputed = true;
		return mi;
	}

	public double[] computeLocalOfPreviousObservations() throws Exception {
		int N = data1.length; // number of observations
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
			double[][] xyNorms = EuclideanUtils.computeNorms(data1, data2, t);
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
		
		mi = digammaK - invK - averageDiGammas + digammaN;
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
