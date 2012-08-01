package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Compute the Mutual Information between two vectors using the Kraskov estimation method.</p>
 * <p>Computes this directly looking at the marginal space for each variable, rather than
 * using the multi-info (or integration) in the marginal spaces.
 * Two child classes actually implement the two algorithms in the Kraskov paper.</p>
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * @see http://dx.doi.org/10.1103/PhysRevE.69.066138
 * 
 * @author Joseph Lizier
 */
public abstract class MutualInfoCalculatorMultiVariateKraskov implements
		MutualInfoCalculatorMultiVariate {

	/**
	 * we compute distances to the kth neighbour
	 */
	protected int k;
	protected double[][] data1;
	protected double[][] data2;
	protected boolean debug;
	protected double mi;
	protected boolean miComputed;
	
	// Storage for the norms from each observation to each other one
	protected double[][] xNorms;
	protected double[][] yNorms;
	// Keep the norms each time (making reordering very quick)
	//  (Should only be set to false for testing)
	public static boolean tryKeepAllPairsNorms = true;
	public static int MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM = 2000;
	
	public final static String PROP_K = "k";
	public final static String PROP_NORM_TYPE = "NORM_TYPE";
	public static final String PROP_NORMALISE = "NORMALISE";
	private boolean normalise = true;
	private int timeDiff = 0;

	public MutualInfoCalculatorMultiVariateKraskov() {
		super();
		k = 1; // by default
	}

	public void initialise(int dimensions1, int dimensions2) {
		mi = 0.0;
		miComputed = false;
		xNorms = null;
		yNorms = null;
		data1 = null;
		data2 = null;
		// No need to keep the dimenions here
	}

	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			EuclideanUtils.setNormToUse(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_TIME_DIFF)) {
			timeDiff = Integer.parseInt(propertyValue);
			if (timeDiff < 0) {
				throw new RuntimeException("Time difference must be >= 0. Flip data1 and data2 around if required.");
			}
		}
	}

	public void addObservations(double[][] source, double[][] destination) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void addObservations(double[][] source, double[][] destination, int startTime, int numTimeSteps) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination, boolean[] sourceValid, boolean[] destValid) throws Exception {
		// Will need to handle timeDiff here
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination, boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		// Will need to handle timeDiff here
		throw new RuntimeException("Not implemented yet");
	}

	public void startAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void finaliseAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] observations1,
			double[][] observations2) throws Exception {
		if (observations1.length != observations2.length) {
			throw new Exception("Time steps for observations2 " +
					observations2.length + " does not match the length " +
					"of observations1 " + observations1.length);
		}
		if ((observations1[0].length == 0) || (observations2[0].length == 0)) {
			throw new Exception("Computing MI with a null set of data");
		}
		// Grab the data
		if (timeDiff == 0) {
			data1 = observations1;
			data2 = observations2;
		} else {
			// MI with time difference: I(x_{1,n}; x_{2,n+timeDiff})
			if (observations1.length < timeDiff) {
				throw new Exception("Computing MI with too few observations " +
						observations1.length + " given time diff " + timeDiff);
			}
			data1 = MatrixUtils.selectRows(observations1, 0, observations1.length - timeDiff);
			data2 = MatrixUtils.selectRows(observations2, timeDiff, observations2.length - timeDiff);
		}
		// Normalise it if required
		if (normalise) {
			// Take a copy since we're going to normalise it
			data1 = MatrixUtils.normaliseIntoNewArray(data1);
			data2 = MatrixUtils.normaliseIntoNewArray(data2);
		}
	}

	/**
	 * Compute the norms for each marginal time series
	 *
	 */
	protected void computeNorms() {
		int N = data1.length; // number of observations
		
		xNorms = new double[N][N];
		yNorms = new double[N][N];
		for (int t = 0; t < N; t++) {
			// Compute the norms from t to all other time points
			double[][] xyNormsForT = EuclideanUtils.computeNorms(data1, data2, t);
			for (int t2 = 0; t2 < N; t2++) {
				xNorms[t][t2] = xyNormsForT[t2][0];
				yNorms[t][t2] = xyNormsForT[t2][1];
			}
		}
	}
	
	public abstract double computeAverageLocalOfObservations() throws Exception;

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
	public abstract double computeAverageLocalOfObservations(int[] reordering) throws Exception;

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
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(data1.length, numPermutationsToCheck);
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
		double actualMI = mi;
		
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
		mi = actualMI;

		// And return the significance
		measDistribution.pValue = (double) countWhereMiIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = mi;
		return measDistribution;
	}

	public abstract double[] computeLocalOfPreviousObservations() throws Exception;

	public double[] computeLocalUsingPreviousObservations(double[][] states1, double[][] states2) throws Exception {
		// If implemented, will need to incorporate any time difference here.
		throw new Exception("Local method not implemented yet");
	}
	
	public abstract String printConstants(int N) throws Exception ;

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return mi;
	}
	
	public int getNumObservations() {
		return data1.length;
	}
}
