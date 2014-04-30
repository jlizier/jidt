package infodynamics.measures.continuous.kraskov;

import java.util.Random;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.MutualInfoMultiVariateCommon;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential mutual information of two given multivariate sets of
 *  observations,
 *  using Kraskov-Grassberger estimation (see Kraskov et al., below).
 *  Two child classes {@link MutualInfoCalculatorMultiVariateKraskov1} and
 *  {@link MutualInfoCalculatorMultiVariateKraskov2}
 *  actually implement the two algorithms in the Kraskov et al. paper</p>
 *  
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct see {@link MutualInfoCalculatorMultiVariateKraskov1#MutualInfoCalculatorMultiVariateKraskov1} </li>
 * 		<li>Set properties using {@link #setProperty(String, String)}</li>
 *		<li>Initialise: {@link #initialise(int, int)}</li>
 * 		<li>Provide the observations to the calculator using:
 * 			{@link #setObservations(double[][], double[][])}
 * 			a sequence of:
 * 			{@link #startAddObservations()},
 *          multiple calls to either {@link #addObservations(double[][], double[][])}
 *          or {@link #addObservations(double[][], double[][], int, int)}, and then
 *          {@link #finaliseAddObservations()}.</li>
 * 		<li>Compute the required information-theoretic results, primarily:
 * 			{@link #computeAverageLocalOfObservations()} to return the average MI
 *          entropy based on the supplied observations; or other calls to compute
 *          local values or statistical significance.</li>
 * 	</ol>
 * </p>
 * 
 * <p>
 * TODO Add fast nearest neighbour searches to the child classes
 * </p>
 * 
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * @see {@link http://dx.doi.org/10.1103/PhysRevE.69.066138}
 * @author Joseph Lizier, <a href="mailto:joseph.lizier at gmail.com">joseph.lizier at gmail.com</a>
 */
public abstract class MutualInfoCalculatorMultiVariateKraskov
	extends MutualInfoMultiVariateCommon
	implements MutualInfoCalculatorMultiVariate {

	/**
	 * we compute distances to the kth neighbour
	 */
	protected int k;
	
	protected EuclideanUtils normCalculator;
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
	/**
	 * Property for an amount of random Gaussian noise to be
	 *  added to the data.
	 */
	public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";
	
	protected boolean normalise = true;
	protected boolean addNoise = false;
	protected double noiseLevel = 0.0;

	public MutualInfoCalculatorMultiVariateKraskov() {
		super();
		k = 1; // by default
		normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
	}

	public void initialise(int sourceDimensions, int destDimensions) {
		super.initialise(sourceDimensions, destDimensions);
		xNorms = null;
		yNorms = null;
	}

	/**
	 * Sets properties for the calculator.
	 * Valid properties include:
	 * <ul>
	 *  <li>{@link #PROP_K} - number of neighbouring points in joint kernel space</li>
	 * 	<li>{@link #PROP_NORM_TYPE}</li> - normalization type to apply to 
	 * 		working out the norms between the points in each marginal space.
	 * 		Options are defined by {@link EuclideanUtils#setNormToUse(String)} -
	 * 		default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 *  <li>{@link #PROP_NORMALISE} - whether to normalise the individual
	 *      variables (true by default)</li>
	 *  <li>{@link #PROP_ADD_NOISE} - an amount of random noise to add to
	 *       each variable, to avoid having neighbourhoods with artificially
	 *       large counts. The amount is added in before any normalisation. 
	 *       (Recommended by Kraskov. MILCA uses 1e-8; but adds in
	 *       a random amount of noise in [0,noiseLevel) ). Default 0.</li>
	 * </ul>
	 * or any properties set in {@link MutualInfoMultiVariateCommon#setProperty(String, String)}.
	 * 
	 * @param propertyName
	 * @param propertyValue
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normCalculator.setNormToUse(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
			addNoise = true;
			noiseLevel = Double.parseDouble(propertyValue);
		} else {
			// No property was set here
			propertySet = false;
			// try the superclass:
			super.setProperty(propertyName, propertyValue);
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.MutualInfoMultiVariateCommon#finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		// Allow the parent to generate the data for us first
		super.finaliseAddObservations();
		
		if (addNoise) {
			Random random = new Random();
			// Add Gaussian noise of std dev noiseLevel to the data
			for (int r = 0; r < sourceObservations.length; r++) {
				for (int c = 0; c < dimensionsSource; c++) {
					sourceObservations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
				for (int c = 0; c < dimensionsDest; c++) {
					destObservations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
			}
		}
		
		// Then normalise the data if required
		if (normalise) {
			// Take a copy since we're going to normalise it
			sourceObservations = MatrixUtils.normaliseIntoNewArray(sourceObservations);
			destObservations = MatrixUtils.normaliseIntoNewArray(destObservations);
		}
	}

	/**
	 * Compute the norms for each marginal time series
	 *
	 */
	protected void computeNorms() {
		int N = sourceObservations.length; // number of observations
		
		xNorms = new double[N][N];
		yNorms = new double[N][N];
		for (int t = 0; t < N; t++) {
			// Compute the norms from t to all other time points
			double[][] xyNormsForT = normCalculator.computeNorms(sourceObservations, destObservations, t);
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

	public abstract double[] computeLocalOfPreviousObservations() throws Exception;

	public double[] computeLocalUsingPreviousObservations(double[][] states1, double[][] states2) throws Exception {
		// TODO If implemented, will need to incorporate any time difference here.
		// Will also need to handle normalisation of the incoming data
		//  appropriately
		throw new Exception("Local method not implemented yet");
	}
	
	public abstract String printConstants(int N) throws Exception ;	
}
