package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.MeasurementDistribution;
import infodynamics.utils.RandomGenerator;

import java.util.Hashtable;


/**
 * <p>Compute the Mutual Information between two vectors using the Kraskov estimation method.
 * Computes this using the multi-info (or integration) in the marginal spaces.
 * Two child classes actually implement the two algorithms in the Kraskov paper.</p>
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * @see http://dx.doi.org/10.1103/PhysRevE.69.066138
 * 
 * @author Joseph Lizier
 */
public abstract class MutualInfoCalculatorMultiVariateKraskovByMulti implements
		MutualInfoCalculatorMultiVariate {

	/**
	 * Storage for the properties ready to pass onto the underlying MI calculators
	 */
	private Hashtable<String,String> props;
	
	/**
	 * Properties for the underlying MultiInfoCalculatorKraskov.
	 * Added here so they can be accessed externally and the accessor doesn't need
	 *  to know that they're really part of the underlying multi-info calculators.
	 */
	public final static String PROP_K = MultiInfoCalculatorKraskov.PROP_K;
	public final static String PROP_NORM_TYPE = MultiInfoCalculatorKraskov.PROP_NORM_TYPE;
	public final static String PROP_TRY_TO_KEEP_ALL_PAIRS_NORM = MultiInfoCalculatorKraskov.PROP_TRY_TO_KEEP_ALL_PAIRS_NORM;

	/**
	 * MultiInfo calculator for the joint space
	 */
	protected MultiInfoCalculatorKraskov multiInfoJoint;
	/**
	 * MultiInfo calculator for marginal space 1
	 */
	protected MultiInfoCalculatorKraskov multiInfo1;
	/**
	 * MultiInfo calculator for marginal space 2
	 */
	protected MultiInfoCalculatorKraskov multiInfo2;

	private double[][] data1;
	private double[][] data2;
	private int dimensions1;
	private int dimensions2;
	private int numObservations;
	
	protected boolean debug;
	protected double mi;
	protected boolean miComputed;
	
	public MutualInfoCalculatorMultiVariateKraskovByMulti() {
		super();
		props = new Hashtable<String,String>();
		createMultiInfoCalculators();
	}

	/**
	 * Create the underlying Kraskov multi info calculators
	 *
	 */
	protected abstract void createMultiInfoCalculators();
	
	public void initialise(int dimensions1, int dimensions2) {
		mi = 0.0;
		miComputed = false;
		numObservations = 0;
		data1 = null;
		data2 = null;
		// Set the properties for the Kraskov multi info calculators
		for (String key : props.keySet()) {
			multiInfoJoint.setProperty(key, props.get(key));
			multiInfo1.setProperty(key, props.get(key));
			multiInfo2.setProperty(key, props.get(key));
		}
		// Initialise the Kraskov multi info calculators
		multiInfoJoint.initialise(dimensions1 + dimensions2);
		multiInfo1.initialise(dimensions1);
		multiInfo2.initialise(dimensions2);
		this.dimensions1 = dimensions1;
		this.dimensions2 = dimensions2;
	}

	/**
	 * Sets properties for the calculator.
	 * Valid properties include:
	 * <ul>
	 * 	<li>Any valid properties for MultiInfoCalculatorKraskov.setProperty</li>
	 * </ul>
	 * One should set MultiInfoCalculatorKraskov.PROP_K here, the number
	 *  of neighbouring points one should count up to in determining the joint kernel size. 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 */
	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_TIME_DIFF)) {
			int diff = Integer.parseInt(propertyValue);
			if (diff != 0) {
				throw new RuntimeException(PROP_TIME_DIFF + " property != 0 not implemented yet");
			}
		}
		// No other local properties here, so 
		// assume it was a property for the MI calculator
		props.put(propertyName, propertyValue);
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

	/**
	 * Set the observations from which to compute the mutual information
	 * 
	 * @param observations1
	 * @param observations2
	 */
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
		data1 = observations1;
		data2 = observations2;
		multiInfoJoint.setObservations(data1, data2);
		multiInfo1.setObservations(data1);
		multiInfo2.setObservations(data2);
		numObservations = data1.length;
	}

	/**
	 * 
	 * @return the average mutual information
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		double jointMultiInfo = multiInfoJoint.computeAverageLocalOfObservations();
		shareNormsIfPossible();
		// Now compute the marginal multi-infos
		double marginal1MultiInfo = multiInfo1.computeAverageLocalOfObservations();
		double marginal2MultiInfo = multiInfo2.computeAverageLocalOfObservations();
		// And return the mutual info
		mi = jointMultiInfo - marginal1MultiInfo - marginal2MultiInfo;
		if (debug) {
			System.out.println("jointMultiInfo=" + jointMultiInfo + " - marginal1MultiInfo=" +
					marginal1MultiInfo + " - marginal2MultiInfo=" + marginal2MultiInfo + 
					" = " + mi);
		}
		miComputed = true;
		return mi;
	}

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
		int[][] reorderingForJointSpace = null;
		int[][] reorderingFor2Space = null;
		if (reordering != null) {
			// We need to make the reordering for the second marginal data set apply 
			//  to each variable within that data set, and keep all variables in the first data
			//  set unchanged
			reorderingForJointSpace = new int[dimensions1 + dimensions2][reordering.length];
			reorderingFor2Space = new int[dimensions2][];
			for (int t = 0; t < numObservations; t++) {
				// Keep the first marginal space not reordered
				for (int v = 0; v < dimensions1; v++) {
					reorderingForJointSpace[v][t] = t;
				}
				// Reorder the second marginal space to match the requested reordering
				for (int v = 0; v < dimensions2; v++) {
					reorderingForJointSpace[v + dimensions1][t] = reordering[t];
				}
			}
			for (int v = 0; v < dimensions2; v++) {
				reorderingFor2Space[v] = reorderingForJointSpace[v + dimensions1];
			}
		}
		double jointMultiInfo = multiInfoJoint.computeAverageLocalOfObservations(reorderingForJointSpace);
		shareNormsIfPossible();
		// Now compute the marginal multi-info in space 1 without reordering
		double marginal1MultiInfo = multiInfo1.computeAverageLocalOfObservations();
		// Now compute the marginal multi-info in space 2 with the reordering applied
		double marginal2MultiInfo = multiInfo2.computeAverageLocalOfObservations(reorderingFor2Space);
		// And return the mutual info
		mi = jointMultiInfo - marginal1MultiInfo - marginal2MultiInfo;
		miComputed = true;
		return mi;
	}

	/**
	 * If the underlying joint space calculator has computed the norms, share them
	 *  with the marginal calculators 
	 *
	 */
	protected void shareNormsIfPossible() {
		if (multiInfoJoint.norms != null) {
			// Share the norms already computed for the joint space:
			if (multiInfo1.norms == null) {
				multiInfo1.norms = new double[dimensions1][][];
				for (int v = 0; v < dimensions1; v++) {
					multiInfo1.norms[v] = multiInfoJoint.norms[v];
				}
			}
			if (multiInfo2.norms == null) {
				multiInfo2.norms = new double[dimensions2][][];
				for (int v = 0; v < dimensions2; v++) {
					multiInfo2.norms[v] = multiInfoJoint.norms[dimensions1 + v];
				}
			}
		}
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
	public synchronized MeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
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
	public MeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception {
		int numPermutationsToCheck = newOrderings.length;
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		// Store the real observations and their MI:
		double actualMI = mi;
		
		MeasurementDistribution measDistribution = new MeasurementDistribution(numPermutationsToCheck);
		
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

	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] localJointMultiInfo = multiInfoJoint.computeLocalOfPreviousObservations();
		shareNormsIfPossible();
		// Now compute the marginal multi-infos
		double[] localMarginal1MultiInfo = multiInfo1.computeLocalOfPreviousObservations();
		MatrixUtils.subtractInPlace(localJointMultiInfo, localMarginal1MultiInfo);
		double[] localMarginal2MultiInfo = multiInfo2.computeLocalOfPreviousObservations();
		MatrixUtils.subtractInPlace(localJointMultiInfo, localMarginal2MultiInfo);
		// And return the mutual info
		mi = multiInfoJoint.getLastAverage() - multiInfo1.getLastAverage() - multiInfo2.getLastAverage();
		miComputed = true;
		return localJointMultiInfo;
	}

	public double[] computeLocalUsingPreviousObservations(double[][] states1, double[][] states2) throws Exception {
		throw new Exception("Local method not implemented yet");
	}
	
	public void setDebug(boolean debug) {
		this.debug = debug;
		multiInfoJoint.debug = debug;
		multiInfo1.debug = debug;
		multiInfo2.debug = debug;
	}

	public double getLastAverage() {
		return mi;
	}

	public String printConstants(int N) throws Exception {
		return multiInfoJoint.printConstants(N);
	}

	public int getNumObservations() {
		return numObservations;
	}
}
