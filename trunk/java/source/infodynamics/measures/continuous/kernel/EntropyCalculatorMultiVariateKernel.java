package infodynamics.measures.continuous.kernel;

import infodynamics.measures.continuous.EntropyCalculatorMultiVariate;

/**
 * Class to compute entropy for
 *  a multi-variate values, using kernel estimates.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public class EntropyCalculatorMultiVariateKernel implements EntropyCalculatorMultiVariate {

	private KernelEstimatorMultiVariate mvke = null;
	private int totalObservations = 0;
	// private int dimensions = 0;
	private boolean debug = false;
	private double[][] observations = null;
	private double lastEntropy;
	
	private boolean normalise = true;
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	/**
	 * Default value for epsilon
	 */
	private static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	private double epsilon = DEFAULT_EPSILON;
	public static final String EPSILON_PROP_NAME = "EPSILON";

	public EntropyCalculatorMultiVariateKernel() {
		mvke = new KernelEstimatorMultiVariate();
		mvke.setDebug(debug);
		mvke.setNormalise(normalise);
		lastEntropy = 0.0;
	}

	/**
	 * Initialises with the default value for epsilon
	 */
	public void initialise(int dimensions) {
		initialise(dimensions, epsilon);
	}

	public void initialise(int dimensions, double epsilon) {
		this.epsilon = epsilon;
		mvke.initialise(dimensions, epsilon);
		// this.dimensions = dimensions;
		lastEntropy = 0.0;
	}

	/**
	 * Set the observations for the PDFs.
	 * Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated). 
	 * 
	 * @param observations
	 */
	public void setObservations(double observations[][]) {
		mvke.setObservations(observations);
		totalObservations = observations.length;
		this.observations = observations;
	}

	public double computeAverageLocalOfObservations() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvke.getProbability(observations[b]);
			double cont = Math.log(prob);
			entropy -= cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + (-cont/Math.log(2.0)) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		lastEntropy = entropy / (double) totalObservations / Math.log(2.0);
		return lastEntropy;
	}

	
	public double[] computeLocalOfPreviousObservations() {
		return computeLocalUsingPreviousObservations(observations);
	}

	public double[] computeLocalUsingPreviousObservations(double states[][]) {
		double entropy = 0.0;
		double[] localEntropy = new double[states.length];
		for (int b = 0; b < states.length; b++) {
			double prob = mvke.getProbability(states[b]);
			double cont = -Math.log(prob);
			localEntropy[b] = cont;
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		entropy /= (double) totalObservations / Math.log(2.0);
		return localEntropy;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
		mvke.setDebug(debug);
	}

	public double getLastAverage() {
		return lastEntropy;
	}
	
	/**
	 * Allows the user to set properties for the underlying calculator implementation
	 * These can include:
	 * <ul>
	 * 		<li>{@link #EPSILON_PROP_NAME}</li>
	 * 		<li>{@link #NORMALISE_PROP_NAME}</li>
	 * </ul> 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;

		// TODO If we implement a dynamic correlation exclusion property,
		//  then we will need to call getProbability(double, int) instead of
		//  just getProbability(double) above.
		
		if (propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			epsilon = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
			mvke.setNormalise(normalise);
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	public int getNumObservations() throws Exception {
		return totalObservations;
	}

}
