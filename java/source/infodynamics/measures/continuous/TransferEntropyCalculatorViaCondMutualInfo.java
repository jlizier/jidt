package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;

import java.util.Vector;


/**
 * 
 * <p>Transfer entropy calculator which is implemented using a 
 * Conditional Mutual Information calculator.
 * </p>
 * 
 * TODO Delete TransferEntropyCalculatorCommon once we've switched everything over to use this?
 * Might be useful to leave it after all, and move common functionality from here to there.
 * 
 * TODO Switch the properties like l etc into TransferEntropyCalculator
 * once all TE calculators follow this class.
 * 
 * @see "Schreiber, Physical Review Letters 85 (2) pp.461-464, 2000;
 *  <a href='http://dx.doi.org/10.1103/PhysRevLett.85.461'>download</a>
 *  (for definition of transfer entropy)"
 * @see "Lizier, Prokopenko and Zomaya, Physical Review E 77, 026110, 2008;
 * <a href='http://dx.doi.org/10.1103/PhysRevE.77.026110'>download</a>
 *  (for definition of <i>local</i> transfer entropy and qualification
 *  of naming it as <i>apparent</i> transfer entropy)"
 *  
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 */
public abstract class TransferEntropyCalculatorViaCondMutualInfo implements
		TransferEntropyCalculator {

	/**
	 * Underlying conditional mutual information calculator
	 */
	protected ConditionalMutualInfoCalculatorMultiVariate condMiCalc;
	/**
	 * Length of past destination history to consider (embedding length)
	 */
	protected int k = 1;
	/**
	 * Embedding delay to use between elements of the destination embeding vector.
	 * We're hard-coding a delay of 1 between the history vector and the next 
	 *  observation however.
	 */
	protected int k_tau = 1;
	/**
	 * Length of past source history to consider (embedding length)
	 */
	protected int l = 1;
	/**
	 * Embedding delay to use between elements of the source embeding vector.
	 */
	protected int l_tau = 1;
	/**
	 * Source-destination next observation delay
	 */
	protected int delay = 1;
	/**
	 * Time index of the last point in the destination embedding of the first
	 *  (destination past, source past, destination next) tuple than can be 
	 *  taken from any set of time-series observations. 
	 */
	protected int startTimeForFirstDestEmbedding;

	protected boolean debug = false;

	// TODO Move these properties up the TransferEntropyCalculator once all the calcs use them
	
	/**
	 * Embedding delay for the destination past history vector
	 */
	public static final String K_TAU_PROP_NAME = "k_TAU";
	/**
	 * Embedding delay for the destination past history vector
	 */
	public static final String L_PROP_NAME = "l_HISTORY";
	/**
	 * Embedding delay for the source past history vector
	 */
	public static final String L_TAU_PROP_NAME = "l_TAU";
	/**
	 * Source-destination delay
	 */
	public static final String DELAY_PROP_NAME = "DELAY";

	
	public TransferEntropyCalculatorViaCondMutualInfo(String condMiCalculatorClassName)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		@SuppressWarnings("unchecked")
		Class<ConditionalMutualInfoCalculatorMultiVariate> condMiClass = 
				(Class<ConditionalMutualInfoCalculatorMultiVariate>) Class.forName(condMiCalculatorClassName);
		ConditionalMutualInfoCalculatorMultiVariate condMiCalc = condMiClass.newInstance();
		construct(condMiCalc);
	}

	public TransferEntropyCalculatorViaCondMutualInfo(Class<ConditionalMutualInfoCalculatorMultiVariate> condMiCalcClass)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		ConditionalMutualInfoCalculatorMultiVariate condMiCalc = condMiCalcClass.newInstance();
		construct(condMiCalc);
	}

	/**
	 * Construct this calculator by passing in a constructed but not initialised
	 * underlying Conditional Mutual information calculator.
	 * 
	 * @param condMiCalc
	 */
	public TransferEntropyCalculatorViaCondMutualInfo(ConditionalMutualInfoCalculatorMultiVariate condMiCalc) {
		construct(condMiCalc);
	}
	
	protected void construct(ConditionalMutualInfoCalculatorMultiVariate condMiCalc) {
		this.condMiCalc = condMiCalc;
	}
	
	public void initialise() throws Exception {
		initialise(k, k_tau, l, l_tau, delay);
	}
	
	public void initialise(int k) throws Exception {
		initialise(k, k_tau, l, l_tau, delay);
	}
	
	public void initialise(int k, int l) throws Exception {
		initialise(k, k_tau, l, l_tau, delay);
	}

	public void initialise(int k, int l, int delay) throws Exception {
		initialise(k, k_tau, l, l_tau, delay);
	}

	/**
	 * Initialise the calculator
	 * 
	 * @param k Length of destination past history to consider
	 * @param k_tau embedding delay for the destination variable
	 * @param l length of source past history to consider
	 * @param l_tau embedding delay for the source variable
	 * @param delay time lag between last element of source and destination next value
	 */
	public void initialise(int k, int k_tau, int l, int l_tau, int delay) throws Exception {
		if (delay < 0) {
			throw new Exception("Cannot compute TE with source-destination delay < 0");
		}
		this.k = k;
		this.k_tau = k_tau;
		this.l = l;
		this.l_tau = l_tau;
		this.delay = delay;
		
		// Now check which point we can start taking observations from in any
		//  addObservations call. These two integers represent the last
		//  point of the destination embedding, in the cases where the destination
		//  embedding itself determines where we can start taking observations, or
		//  the case where the source embedding plus delay is longer and so determines
		//  where we can start taking observations.
		int startTimeBasedOnDestPast = (k-1)*k_tau;
		int startTimeBasedOnSourcePast = (l-1)*l_tau + delay - 1;
		startTimeForFirstDestEmbedding = Math.max(startTimeBasedOnDestPast, startTimeBasedOnSourcePast);

		condMiCalc.initialise(l, 1, k);
	}

	/**
	 * <p>Set the given property to the given value.
	 * These can include:
	 * <ul>
	 * 		<li>{@link #K_PROP_NAME}</li>
	 * 		<li>{@link #K_TAU_PROP_NAME}</li>
	 * 		<li>{@link #L_TAU_PROP_NAME}</li>
	 * 		<li>{@link #DELAY_PROP_NAME}</li>
	 * </ul>
	 * Or else it is assumed the property
	 *  is for the underlying {@link ConditionalMutualInfoCalculatorMultiVariate#setProperty(String, String)} implementation.
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(K_TAU_PROP_NAME)) {
			k_tau = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			l = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(L_TAU_PROP_NAME)) {
			l_tau = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(DELAY_PROP_NAME)) {
			delay = Integer.parseInt(propertyValue);
		} else {
			// No property was set on this class, assume it is for the underlying
			//  conditional MI calculator
			condMiCalc.setProperty(propertyName, propertyValue);
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @throws Exception
	 */
	public void setObservations(double[] source, double[] destination) throws Exception {
		
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to add here, the time series is too short
			throw new Exception("Not enough observations to set here given k, k_tau, l, l_tau and delay parameters");
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(destination, k, k_tau,
						startTimeForFirstDestEmbedding,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(destination, 1,
						startTimeForFirstDestEmbedding + 1,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(source, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						source.length - startTimeForFirstDestEmbedding - 1);
		condMiCalc.setObservations(currentSourcePastVectors, currentDestNextVectors, currentDestPastVectors);
	}

	public void startAddObservations() {
		condMiCalc.startAddObservations();
	}
	
	public void addObservations(double[] source, double[] destination) throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to add here, the time series is too short
			// Don't throw an exception, do nothing since more observations
			//  can be added later.
			return;
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(destination, k, k_tau,
						startTimeForFirstDestEmbedding,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(destination, 1,
						startTimeForFirstDestEmbedding + 1,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(source, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						source.length - startTimeForFirstDestEmbedding - 1);
		condMiCalc.addObservations(currentSourcePastVectors, currentDestNextVectors, currentDestPastVectors);
	}

	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length < startTime + numTimeSteps) {
			// There are not enough observations given the arguments here
			throw new Exception("Not enough observations to set here given startTime and numTimeSteps parameters");
		}
		addObservations(MatrixUtils.select(source, startTime, numTimeSteps),
					    MatrixUtils.select(destination, startTime, numTimeSteps));
	}

	public void finaliseAddObservations() throws Exception {
		condMiCalc.finaliseAddObservations();
	}

	public void setObservations(double[] source, double[] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		
		Vector<int[]> startAndEndTimePairs = computeStartAndEndTimePairs(sourceValid, destValid);
		
		// We've found the set of start and end times for this pair
		startAddObservations();
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservations(source, destination, startTime, endTime - startTime + 1);
		}
		finaliseAddObservations();
	}

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid series of both source and destinations.
	 * 
	 * Made public so it can be used if one wants to compute the number of
	 *  observations prior to setting the observations.
	 * 
	 * @param sourceValid
	 * @param destValid
	 * @return
	 * @throws Exception 
	 */
	public Vector<int[]> computeStartAndEndTimePairs(boolean[] sourceValid, boolean[] destValid) throws Exception {
		
		if (sourceValid.length != destValid.length) {
			throw new Exception("Validity arrays must be of same length");
		}
		
		int lengthOfDestPastRequired = (k-1)*k_tau + 1;
		int lengthOfSourcePastRequired = (l-1)*l_tau + 1;
		// int numSourcePointsBeforeDestStart = delay - 1 + lengthOfSourcePastRequired
		//									- lengthOfDestPastRequired;

		// Scan along the data avoiding invalid values
		int startTime = 0;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();

		// Simple solution -- this takes more complexity in time, but is 
		//  much faster to code:
		boolean previousWasOk = false;
		for (int t = startTimeForFirstDestEmbedding; t < destValid.length - 1; t++) {
			// Check the tuple with the history vector starting from
			//  t and running backwards
			if (previousWasOk) {
				// Just check the very next values of each:
				if (destValid[t + 1] && sourceValid[t + 1 - delay]) {
					// We can continue adding to this sequence
					continue;
				} else {
					// We need to shut down this sequence now
					previousWasOk = false;
					int[] timePair = new int[2];
					timePair[0] = startTime;
					timePair[1] = t; // Previous time step was last valid one
					startAndEndTimePairs.add(timePair);
					continue;
				}
			}
			// Otherwise we're trying to start a new sequence, so check all values
			if (!destValid[t + 1]) {
				continue;
			}
			boolean allOk = true;
			for (int tBack = 0; tBack < lengthOfDestPastRequired; tBack++) {
				if (!destValid[t - tBack]) {
					allOk = false;
					break;
				}
			}
			if (!allOk) {
				continue;
			}
			allOk = true;
			for (int tBack = delay - 1; tBack < delay - 1 + lengthOfSourcePastRequired; tBack++) {
				if (!sourceValid[t - tBack]) {
					allOk = false;
					break;
				}
			}
			if (!allOk) {
				continue;
			}
			// Postcondition: We've got a first valid tuple:
			startTime = t - startTimeForFirstDestEmbedding;
			previousWasOk = true;
		}
		// Now check if we were running a sequence and terminate it:
		if (previousWasOk) {
			// We need to shut down this sequence now
			previousWasOk = false;
			int[] timePair = new int[2];
			timePair[0] = startTime;
			timePair[1] = destValid.length - 1;
			startAndEndTimePairs.add(timePair);
		}

		return startAndEndTimePairs;
	}
	
	public double computeAverageLocalOfObservations() throws Exception {
		return condMiCalc.computeAverageLocalOfObservations();
	}

	/**
	 * Returns a time series of local TE values.
	 * Pads the first (k-1)*tau + 1 elements with zeros (since AIS is undefined here)
	 *  if only one time series of observations was used.
	 * Otherwise, local values for all separate series are concatenated, and without
	 *  padding of zeros at the start.
	 *  
	 * @return an array of local TE values of the previously submitted observations.
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] local = condMiCalc.computeLocalOfPreviousObservations();
		if (!condMiCalc.getAddedMoreThanOneObservationSet()) {
			double[] localsToReturn = new double[local.length + startTimeForFirstDestEmbedding + 1];
			System.arraycopy(local, 0, localsToReturn, startTimeForFirstDestEmbedding + 1, local.length);
			return localsToReturn;
		} else {
			return local;
		}
	}

	public double[] computeLocalUsingPreviousObservations(double[] newSourceObservations, double[] newDestObservations) throws Exception {
		if (newSourceObservations.length != newDestObservations.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					newSourceObservations.length, newDestObservations.length));
		}
		if (newDestObservations.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to compute for here
			return new double[newDestObservations.length];
		}
		double[][] newDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newDestObservations, k, k_tau,
						startTimeForFirstDestEmbedding,
						newDestObservations.length - startTimeForFirstDestEmbedding - 1);
		double[][] newDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(newDestObservations, 1,
						startTimeForFirstDestEmbedding + 1,
						newDestObservations.length - startTimeForFirstDestEmbedding - 1);
		double[][] newSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newSourceObservations, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						newSourceObservations.length - startTimeForFirstDestEmbedding - 1);
		double[] local = condMiCalc.computeLocalUsingPreviousObservations(
						newSourcePastVectors, newDestNextVectors, newDestPastVectors);
		// Pad the front of the array with zeros where local TE isn't defined:
		double[] localsToReturn = new double[local.length + startTimeForFirstDestEmbedding + 1];
		System.arraycopy(local, 0, localsToReturn, startTimeForFirstDestEmbedding + 1, local.length);
		return localsToReturn;

	}
	
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		return condMiCalc.computeSignificance(1, numPermutationsToCheck); // Reorder the source
	}

	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		return condMiCalc.computeSignificance(1, newOrderings); // Reorder the source
	}

	public double getLastAverage() {
		return condMiCalc.getLastAverage();
	}

	public int getNumObservations() throws Exception {
		return condMiCalc.getNumObservations();
	}
	
	public boolean getAddedMoreThanOneObservationSet() {
		return condMiCalc.getAddedMoreThanOneObservationSet();
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
		condMiCalc.setDebug(debug);
	}
}
