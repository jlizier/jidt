package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.utils.EuclideanUtils;

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
 * 
 * @author Joseph Lizier
 */
public abstract class ConditionalMutualInfoCalculatorMultiVariateKraskov 
	extends ConditionalMutualInfoMultiVariateCommon
	implements Cloneable { // See comments on clonability below
	
	/**
	 * we compute distances to the kth neighbour in the joint space
	 */
	protected int k;
		
	protected EuclideanUtils normCalculator;
	// Storage for the norms from each observation to each other one
	protected double[][] xNorms;
	protected double[][] yNorms;
	protected double[][] zNorms;
	// Keep the norms each time (making reordering very quick)
	//  (Should only be set to false for testing)
	public static boolean tryKeepAllPairsNorms = true;
	public static int MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM = 2000;
	
	/**
	 * Property name for the number of nearest neighbours k to use in the Kraskov algorithm in
	 *  the full joint space.
	 */
	public final static String PROP_K = "k";
	/**
	 * Normalisation to apply to the marginal spaces.
	 */
	public final static String PROP_NORM_TYPE = "NORM_TYPE";

	public ConditionalMutualInfoCalculatorMultiVariateKraskov() {
		super();
		k = 1; // by default
		normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
	}

	public void initialise(int dimensions1, int dimensions2, int dimensionsCond) {
		super.initialise(dimensions1, dimensions2, dimensionsCond);
		
		xNorms = null;
		yNorms = null;
		zNorms = null;
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
	 * </ul>
	 * 
	 * @param propertyName name of the property to set
	 * @param propertyValue value to set on that property
	 */
	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normCalculator.setNormToUse(propertyValue);
		} else {
			// Assume this is a property for the common parent class
			super.setProperty(propertyName, propertyValue);
		}
	}

	/**
	 * Compute the norms for each time series
	 *
	 */
	protected void computeNorms() {
		int N = var1Observations.length; // number of observations
		
		xNorms = new double[N][N];
		yNorms = new double[N][N];
		zNorms = new double[N][N];
		for (int t = 0; t < N; t++) {
			// Compute the norms from t to all other time points
			double[][] xyzNormsForT = normCalculator.computeNorms(var1Observations,
					var2Observations, condObservations, t);
			for (int t2 = 0; t2 < N; t2++) {
				xNorms[t][t2] = xyzNormsForT[t2][0];
				yNorms[t][t2] = xyzNormsForT[t2][1];
				zNorms[t][t2] = xyzNormsForT[t2][2];
			}
		}
	}
	
	public double[] computeLocalUsingPreviousObservations(double[][] states1,
			double[][] states2, double[][] condStates) throws Exception {
		// If implemented, will need to incorporate any normalisation here
		//  (normalising the incoming data the same way the previously
		//   supplied observations were normalised).
		throw new Exception("Local method not implemented yet");
	}
	
	public abstract String printConstants(int N) throws Exception;

	// Note: no extra implementation of clone provided; we're simply
	//  allowing clone() to produce a shallow copy, which is find
	//  for the statistical significance calculation (none of the array
	//  data will be changed there.
	//
	// public ConditionalMutualInfoCalculatorMultiVariateKraskov clone() {
	//	return this;
	// }
}
