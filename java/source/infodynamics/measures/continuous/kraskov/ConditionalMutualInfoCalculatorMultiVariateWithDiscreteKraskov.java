package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSource;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Compute the Conditional Mutual Information between a discrete variable and a 
 *  vector of continuous variables, conditioned on another vector of continuous variables
 *  using the Kraskov estimation method.</p>
 * <p>Uses Kraskov method type 2, since type 1 only looks at points with
 * distances strictly less than the kth variable, which won't work for one marginal
 * being discrete.</p>
 * 
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * @see http://dx.doi.org/10.1103/PhysRevE.69.066138
 * 
 * @author Joseph Lizier
 */
public class ConditionalMutualInfoCalculatorMultiVariateWithDiscreteKraskov implements ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSource {

	// Multiplier used in hueristic for determining whether to use a linear search
	//  for min kth element or a binary search.
	protected static final double CUTOFF_MULTIPLIER = 1.5;
	
	/**
	 * we compute distances to the kth neighbour
	 */
	protected int k;
	protected double[][] continuousDataX;
	protected double[][] conditionedDataZ;
	protected int[] discreteDataY;
	protected int[] counts;
	protected int base;
	protected boolean debug;
	protected double condMi;
	protected boolean miComputed;
	
	// Storage for the norms from each observation to each other one
	protected double[][] xNorms;
	protected double[][] zNorms;
	protected double[][] xzNorms;
	// Keep the norms each time (making reordering very quick)
	//  (Should only be set to false for testing)
	public static boolean tryKeepAllPairsNorms = true;
	public static int MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM = 2000;
	
	public final static String PROP_K = "k";
	public final static String PROP_NORM_TYPE = "NORM_TYPE";	
	public static final String PROP_NORMALISE = "NORMALISE";
	private boolean normalise = true;

	public ConditionalMutualInfoCalculatorMultiVariateWithDiscreteKraskov() {
		super();
		k = 1; // by default
	}

	/**
	 * Initialise the calculator.
	 * 
	 * @param dimensions number of joint continuous variables
	 * @param base number of discrete states
	 * @param dimensionsCond the number of joint continuous variables
	 *     to condition on
	 */
	public void initialise(int dimensions, int base, int dimensionsCond) {
		condMi = 0.0;
		miComputed = false;
		xNorms = null;
		continuousDataX = null;
		discreteDataY = null;
		// No need to keep the dimenions for the conditional variables here
		this.base = base;
	}

	/**
	 * 
	 * @param propertyName name of the property to set
	 * @param propertyValue value to set on that property
	 */
	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			EuclideanUtils.setNormToUse(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		}
	}

	public void addObservations(double[][] source, double[][] destination) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void addObservations(double[][] source, double[][] destination, int startTime, int numTimeSteps) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination, boolean[] sourceValid, boolean[] destValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination, boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void startAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void finaliseAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations, double[][] conditionedObservations)
			throws Exception {
		
		if ((continuousObservations.length != discreteObservations.length) ||
			(continuousObservations.length != conditionedObservations.length)) {
			throw new Exception("Time steps for observations2 " +
					discreteObservations.length + " does not match the length " +
					"of observations1 " + continuousObservations.length +
					" and of conditionedObservations " + conditionedObservations.length);
		}
		if (continuousObservations[0].length == 0) {
			throw new Exception("Computing MI with a null set of data");
		}
		if (conditionedObservations[0].length == 0) {
			throw new Exception("Computing MI with a null set of conditioned data");
		}
		continuousDataX = continuousObservations;
		discreteDataY = discreteObservations;
		conditionedDataZ = conditionedObservations;
		if (normalise) {
			// Take a copy since we're going to normalise it
			continuousDataX = MatrixUtils.normaliseIntoNewArray(continuousObservations);
			conditionedDataZ = MatrixUtils.normaliseIntoNewArray(conditionedObservations);
		}
		// count the discrete states:
		counts = new int[base];
		for (int t = 0; t < discreteDataY.length; t++) {
			counts[discreteDataY[t]]++;
		}
		for (int b = 0; b < counts.length; b++) {
			if (counts[b] < k) {
				throw new RuntimeException("This implementation assumes there are at least k items in each discrete bin");
			}
		}
	}

	/**
	 * Compute the norms for each marginal time series
	 *
	 */
	protected void computeNorms() {
		int N = continuousDataX.length; // number of observations
		
		xNorms = new double[N][N];
		zNorms = new double[N][N];
		xzNorms = new double[N][N];
		for (int t = 0; t < N; t++) {
			// Compute the norms from t to all other time points
			for (int t2 = 0; t2 < N; t2++) {
				if (t2 == t) {
					xNorms[t][t2] = Double.POSITIVE_INFINITY;
					zNorms[t][t2] = Double.POSITIVE_INFINITY;
					xzNorms[t][t2] = Double.POSITIVE_INFINITY;
					continue;
				}
				// Compute norm in the continuous space
				xNorms[t][t2] = EuclideanUtils.norm(continuousDataX[t], continuousDataX[t2]);
				zNorms[t][t2] = EuclideanUtils.norm(conditionedDataZ[t], conditionedDataZ[t2]);
				xzNorms[t][t2] = Math.max(xNorms[t][t2], zNorms[t][t2]);
			}
		}
	}
	
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
		int N = continuousDataX.length; // number of observations
		if (!tryKeepAllPairsNorms || (N > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			// Generate a new re-ordered set of discrete data
			int[] originalDiscreteData = discreteDataY;
			discreteDataY = MatrixUtils.extractSelectedTimePoints(discreteDataY, reordering);
			// Compute the MI
			double newMI = computeAverageLocalOfObservationsWhileComputingDistances();
			// restore data2
			discreteDataY = originalDiscreteData;
			return newMI;
		}
		
		// Otherwise we will use the norms we've already computed, and use a "virtual"
		//  reordered data2.
		int[] reorderedDiscreteData = MatrixUtils.extractSelectedTimePoints(discreteDataY, reordering);

		if (xNorms == null) {
			computeNorms();
		}

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double averageInverseCountInJointYZ = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_z for this time step:
			//  using max of x and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).

			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][0] = Math.max(xNorms[t][t2], zNorms[t][t2]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][1] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_z = 0.0;
			int[] timeStepsOfKthMins = null;
			// just do a linear search for the minimum epsilon value
			timeStepsOfKthMins = MatrixUtils.kMinIndicesSubjectTo(
					jointNorm, 0, k, reorderedDiscreteData, reorderedDiscreteData[t]);
			// and now we have the closest k points.
			// Find eps_{x,y,z} as the maximum x and y and z norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xNorms[t][timeStepOfJthPoint] > eps_x) {
					eps_x = xNorms[t][timeStepOfJthPoint];
				}
				if (zNorms[t][timeStepOfJthPoint] > eps_z) {
					eps_z = zNorms[t][timeStepOfJthPoint];
				}
			}

			// Count the number of points whose distances are less
			//  than or equal to eps in each required joint space
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (zNorms[t][t2] <= eps_z) {
					n_z++;
					if (xNorms[t][t2] <= eps_x) {
						n_xz++;
					}
					if (reorderedDiscreteData[t] == reorderedDiscreteData[t2]) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_xz) + MathsUtils.digamma(n_yz)
					- MathsUtils.digamma(n_z);
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
		
		condMi = MathsUtils.digamma(k) - 1.0/(double)k +
			averageInverseCountInJointYZ - averageDiGammas;
		miComputed = true;
		return condMi;
	}

	public double computeAverageLocalOfObservations() throws Exception {
		if (!tryKeepAllPairsNorms || (continuousDataX.length > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			return computeAverageLocalOfObservationsWhileComputingDistances();
		}
		
		if (xNorms == null) {
			computeNorms();
		}
		int N = continuousDataX.length; // number of observations

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double averageInverseCountInJointYZ = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_z for this time step:
			//  using x,z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).

			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][0] = Math.max(xNorms[t][t2], zNorms[t][t2]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][1] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_z = 0.0;
			int[] timeStepsOfKthMins = null;
			// just do a linear search for the minimum epsilon value
			timeStepsOfKthMins = MatrixUtils.kMinIndicesSubjectTo(
					jointNorm, 0, k, discreteDataY, discreteDataY[t]);
			// and now we have the closest k points.
			// Find eps_{x,y,z} as the maximum x and y and z norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xNorms[t][timeStepOfJthPoint] > eps_x) {
					eps_x = xNorms[t][timeStepOfJthPoint];
				}
				if (zNorms[t][timeStepOfJthPoint] > eps_z) {
					eps_z = zNorms[t][timeStepOfJthPoint];
				}
			}

			// Count the number of points whose x,y,z distances are less
			//  than or equal to eps (not including this point)
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (zNorms[t][t2] <= eps_z) {
					n_z++;
					if (xNorms[t][t2] <= eps_x) {
						n_xz++;
					}
					if (discreteDataY[t] == discreteDataY[t2]) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_xz) + MathsUtils.digamma(n_yz)
						- MathsUtils.digamma(n_z);
			averageInverseCountInJointYZ += 1.0 / (double) n_yz;
		}
		averageDiGammas /= (double) N;
		averageInverseCountInJointYZ /= (double) N;
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double) N;
			System.out.printf("Average n_xz=%.3f (-> digam=%.3f %.3f), Average n_yz=%.3f (-> digam=%.3f)",
					avNxz, MathsUtils.digamma((int) avNxz), MathsUtils.digamma((int) avNxz - 1), avNyz, MathsUtils.digamma((int) avNyz));
			System.out.printf(", Average n_z=%.3f (-> digam=%.3f)\n", avNz, MathsUtils.digamma((int) avNz));
			System.out.printf("Independent average num in joint box is %.3f\n", (avNxz * avNyz / (double) N));
			System.out.println(String.format("digamma(k)=%.3f - 1/k=%.3f + <1/n_yz>=%.3f - averageDiGammas=%.3f\n",
					MathsUtils.digamma(k), 1.0/(double)k, averageInverseCountInJointYZ,
					averageDiGammas));
		}
		
		condMi = MathsUtils.digamma(k) - 1.0/(double)k +
			averageInverseCountInJointYZ - averageDiGammas;
		miComputed = true;
		return condMi;
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
		int N = continuousDataX.length; // number of observations

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double averageInverseCountInJointYZ = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_* for this time step:
			//  First get xz norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xzNorms = EuclideanUtils.computeNorms(continuousDataX, conditionedDataZ, t);
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][0] = Math.max(xzNorms[t2][0], 
						xzNorms[t2][1]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][1] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_z = 0.0;
			int[] timeStepsOfKthMins = null;
			// just do a linear search for the minimum epsilon value
			//  subject to the discrete variable value
			timeStepsOfKthMins = MatrixUtils.kMinIndicesSubjectTo(
					jointNorm, 0, k, discreteDataY, discreteDataY[t]);
			// and now we have the closest k points.
			// Find eps_{x,y} as the maximum x and y norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j];
				if (xzNorms[timeStepOfJthPoint][0] > eps_x) {
					eps_x = xzNorms[timeStepOfJthPoint][0];
				}
				if (xzNorms[timeStepOfJthPoint][1] > eps_z) {
					eps_z = xzNorms[timeStepOfJthPoint][1];
				}
			}

			// Count the number of points whose distances is less
			//  than or equal to eps
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xzNorms[t2][1] <= eps_z) {
					n_z++;
					if (xzNorms[t2][0] <= eps_x) {
						n_xz++;
					}
					if (discreteDataY[t] == discreteDataY[t2]) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_xz) + MathsUtils.digamma(n_yz)
				- MathsUtils.digamma(n_z);
		}
		averageDiGammas /= (double) N;
		averageInverseCountInJointYZ /= (double) N;
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.printf("Average n_xz=%.3f, Average n_yz=%.3f, Average n_z=%.3f\n",
					avNxz, avNyz, avNz);
		}
		
		condMi = MathsUtils.digamma(k) - 1.0/(double)k +
			averageInverseCountInJointYZ - averageDiGammas;
		miComputed = true;
		return condMi;
	}

	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param numPermutationsToCheck
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public synchronized EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(continuousDataX.length, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param newOrderings the specific new orderings to use
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		// Store the real observations and their MI:
		double actualMI = condMi;
		
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		
		int countWhereMiIsMoreSignificantThanOriginal = 0;
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Compute the MI under this reordering
			double newMI = computeAverageLocalOfObservations(newOrderings[i]);
			measDistribution.distribution[i] = newMI;
			if (debug){
				System.out.println("New MI was " + newMI);
			}
			if (newMI >= actualMI) {
				countWhereMiIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the actual MI and the observations
		condMi = actualMI;

		// And return the significance
		measDistribution.pValue = (double) countWhereMiIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = condMi;
		return measDistribution;
	}

	public double[] computeLocalUsingPreviousObservations(double[][] continuousStates,
			int[] discreteStates) throws Exception {
		throw new Exception("Local method not implemented yet");
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return condMi;
	}
	
	public int getNumObservations() {
		return continuousDataX.length;
	}

	public double[] computeLocalUsingPreviousObservations(
			double[][] contStates, int[] discreteStates,
			double[][] conditionedStates) throws Exception {
		// TODO Auto-generated method stub
		throw new Exception("Not implemented yet");
	}
}
