package infodynamics.measures.continuous.kernel;

import infodynamics.measures.continuous.EntropyCalculator;

public class EntropyCalculatorKernel implements EntropyCalculator {

	protected KernelEstimatorSingleVariate svke = null;
	protected int totalObservations = 0;
	protected boolean debug = false;
	protected double[] observations;
	
	private boolean normalise = true;
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	/**
	 * Default value for epsilon
	 */
	public static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	private double epsilon = DEFAULT_EPSILON;
	public static final String EPSILON_PROP_NAME = "EPSILON";
	
	public EntropyCalculatorKernel() {
		svke = new KernelEstimatorSingleVariate();
		svke.setDebug(debug);
		svke.setNormalise(normalise);
	}

	public void initialise() {
		initialise(epsilon);
	}

	public void initialise(double epsilon) {
		this.epsilon = epsilon;
		svke.initialise(epsilon);
	}

	/**
	 * Set the observations for the PDFs.
	 * Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated). 
	 * 
	 * @param observations
	 */
	public void setObservations(double observations[]) {
		this.observations = observations;
		svke.setObservations(observations);
		totalObservations = observations.length;
	}
	
	public double computeAverageLocalOfObservations() {
		double entropy = 0.0;
		for (int t = 0; t < observations.length; t++) {
			double prob = svke.getProbability(observations[t]);
			double cont = Math.log(prob);
			entropy -= cont;
			if (debug) {
				System.out.println(t + ": p(" + observations[t] + ")= " +
						prob + " -> " + (cont/Math.log(2.0)) + " -> sum: " +
						(entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}
	
	public void setDebug(boolean debug) {
		this.debug = debug;
		if (svke != null) {
			svke.setDebug(debug);
		}
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
	 * @throws Exception
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
			svke.setNormalise(normalise);
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

}
