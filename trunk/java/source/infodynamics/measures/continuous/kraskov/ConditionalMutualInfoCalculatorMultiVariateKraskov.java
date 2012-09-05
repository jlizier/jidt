package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Compute the Conditional Mutual Information between two vectors,
 * conditioned on a third, using the Kraskov estimation method,
 * as extended by Frenzel and Pompe.</p>
 * <p>Computes this directly looking at the marginal space for each variable, rather than
 * using the multi-info (or integration) in the marginal spaces.
 * Two child classes actually implement the two algorithms in the Kraskov paper.</p>
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * @see http://dx.doi.org/10.1103/PhysRevE.69.066138
 * @see "Partial Mutual Information for Coupling Analysis of Multivariate Time Series", Frenzel and Pompe, 2007
 * 
 * TODO Finish writing this class - changing it from original Kraskov one
 * 
 * @author Joseph Lizier
 */
public abstract class ConditionalMutualInfoCalculatorMultiVariateKraskov {

	/**
	 * we compute distances to the kth neighbour in the joint space
	 */
	protected int k;
	protected double[][] data1;
	protected double[][] data2;
	protected double[][] dataCond;
	protected boolean debug;
	protected double condMi;
	protected boolean condMiComputed;
	
	// Storage for the norms from each observation to each other one
	protected double[][] xNorms;
	protected double[][] yNorms;
	protected double[][] zNorms;
	// Keep the norms each time (making reordering very quick)
	//  (Should only be set to false for testing)
	public static boolean tryKeepAllPairsNorms = true;
	public static int MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM = 2000;
	
	public final static String PROP_K = "k";
	public final static String PROP_NORM_TYPE = "NORM_TYPE";
	public static final String PROP_NORMALISE = "NORMALISE";
	private boolean normalise = true;

	public ConditionalMutualInfoCalculatorMultiVariateKraskov() {
		super();
		k = 1; // by default
	}

	public void initialise(int dimensions1, int dimensions2, int dimensions3) {
		condMi = 0.0;
		condMiComputed = false;
		xNorms = null;
		yNorms = null;
		zNorms = null;
		data1 = null;
		data2 = null;
		// No need to keep the dimensions here
	}

	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			EuclideanUtils.setNormToUse(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		}
	}

	public void addObservations(double[][] var1, double[][] var2,
			double[][] conditionedVar) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void addObservations(double[][] var1, double[][] var2,
			double[][] conditionedVar, int startTime, int numTimeSteps) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] var1, double[][] var2,
			double[][] conditionedVar, boolean[] var1Valid,
			boolean[] var2Valid, boolean[] conditionedValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] var1, double[][] var2,
			double[][] conditionedVar, boolean[][] var1Valid,
			boolean[][] var2Valid, boolean[][] conditionedValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void startAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void finaliseAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] observations1,
			double[][] observations2, double[][] obsConditioned) throws Exception {
		if (observations1.length != observations2.length) {
			throw new Exception("Time steps for observations2 " +
					observations2.length + " does not match the length " +
					"of observations1 " + observations1.length);
		}
		if ((observations1[0].length == 0) || (observations2[0].length == 0)) {
			throw new Exception("Computing MI with a null set of data");
		}
		// Normalise it if required
		if (normalise) {
			// Take a copy since we're going to normalise it
			data1 = MatrixUtils.normaliseIntoNewArray(observations1);
			data2 = MatrixUtils.normaliseIntoNewArray(observations2);
			dataCond = MatrixUtils.normaliseIntoNewArray(obsConditioned);
		} else {
			data1 = observations1;
			data2 = observations2;
			dataCond = obsConditioned;
		}
	}

	/**
	 * Compute the norms for each time series
	 *
	 */
	protected void computeNorms() {
		int N = data1.length; // number of observations
		
		xNorms = new double[N][N];
		yNorms = new double[N][N];
		zNorms = new double[N][N];
		for (int t = 0; t < N; t++) {
			// Compute the norms from t to all other time points
			double[][] xyzNormsForT = EuclideanUtils.computeNorms(data1,
					data2, dataCond, t);
			for (int t2 = 0; t2 < N; t2++) {
				xNorms[t][t2] = xyzNormsForT[t2][0];
				yNorms[t][t2] = xyzNormsForT[t2][1];
				zNorms[t][t2] = xyzNormsForT[t2][2];
			}
		}
	}
	
	public abstract double computeAverageLocalOfObservations() throws Exception;

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
	public abstract double computeAverageLocalOfObservations(int[] reordering) throws Exception;

	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y,z) correlations, while retaining the p(x,z), p(y) marginals, to check how
	 *  significant this conditional mutual information actually was.
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
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(data1.length, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y,z) correlations, while retaining the p(x,z), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param newOrderings the specific new orderings to use
	 * @return the proportion of conditional MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception {
		int numPermutationsToCheck = newOrderings.length;
		if (!condMiComputed) {
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

	public abstract double[] computeLocalOfPreviousObservations() throws Exception;

	public double[] computeLocalUsingPreviousObservations(double[][] states1, double[][] states2) throws Exception {
		// If implemented, will need to incorporate any normalisation here
		//  (normalising the incoming data the same way the previously
		//   supplied observations were normalised).
		throw new Exception("Local method not implemented yet");
	}
	
	public abstract String printConstants(int N) throws Exception ;

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return condMi;
	}
	
	public int getNumObservations() {
		return data1.length;
	}
}
