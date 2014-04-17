package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.ParsedProperties;

import java.util.Vector;

/**
 * <p>An abstract Conditional Transfer entropy calculator which is implemented using a 
 * Conditional Mutual Information calculator.
 * The Conditional Mutual Information calculator must be supplied at construction time.
 * </p>
 * 
 * <p>There are no abstract methods of this class, and conceivably it could be constructed
 * with a {@link ConditionalMutualInfoCalculatorMultiVariate} class supplied, 
 * however there are typically extra considerations for each estimator type.
 * As such, the children of this abstract class provide concrete implementations using various
 * estimator types; see e.g. {@link infodynamics.continuous.gaussian.ConditionalTransferEntropyCalculatorGaussian}.
 * </p>
 * 
 * @see "Schreiber, Physical Review Letters 85 (2) pp.461-464 (2000);
 *  <a href='http://dx.doi.org/10.1103/PhysRevLett.85.461'>download</a>
 *  (for definition of transfer entropy)"
 * @see "Lizier, Prokopenko and Zomaya, Physical Review E 77, 026110 (2008);
 * <a href='http://dx.doi.org/10.1103/PhysRevE.77.026110'>download</a>
 *  (for the extension to <i>conditional</i> transfer entropy 
 *  or <i>complete</i> where all other causal sources are conditioned on,
 *  and <i>local</i> transfer entropy)"
 * @see "Lizier, Prokopenko and Zomaya, Chaos 20, 3, 037109 (2010);
 * <a href='http://dx.doi.org/10.1063/1.3486801'>download</a>
 *  (for further clarification on <i>conditional</i> transfer entropy 
 *  or <i>complete</i> where all other causal sources are conditioned on)"
 *  
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 */
public abstract class ConditionalTransferEntropyCalculatorViaCondMutualInfo implements
		ConditionalTransferEntropyCalculator {

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
	 * Array of embedding lengths for each conditional variable.
	 *  Can be an empty array or null if there are no conditional variables.
	 */
	protected int[] condEmbedDims = null;
	/**
	 * Array of embedding delays for the conditional variables.
	 *  Must be same length as condEmbedDims array.
	 */
	protected int[] cond_taus = null;
	/**
	 * Array of time lags between last element of each conditional variable
	 *  and destination next value.
	 */
	protected int[] condDelays = null;
	
	/**
	 * Time index of the last point in the destination embedding of the first
	 *  (destination past, source past, destination next) tuple than can be 
	 *  taken from any set of time-series observations. 
	 */
	protected int startTimeForFirstDestEmbedding;

	/**
	 * The total dimensionality of our embedded conditional values
	 *  (sum of condEmbedDims)
	 */
	protected int dimOfConditionals = 0;
	
	protected boolean debug = false;
	
	/**
	 * Construct a conditional transfer entropy calculator using an instance of
	 * condMiCalculatorClassName as the underlying conditional mutual information calculator.
	 * 
	 * @param condMiCalculatorClassName name of the class which must implement
	 * 	{@link ConditionalMutualInfoCalculatorMultiVariate}
	 * @throws InstantiationException if the given class cannot be instantiated
	 * @throws IllegalAccessException if illegal access occurs while trying to create an instance
	 *   of the class
	 * @throws ClassNotFoundException if the given class is not found
	 */
	public ConditionalTransferEntropyCalculatorViaCondMutualInfo(String condMiCalculatorClassName)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		@SuppressWarnings("unchecked")
		Class<ConditionalMutualInfoCalculatorMultiVariate> condMiClass = 
				(Class<ConditionalMutualInfoCalculatorMultiVariate>) Class.forName(condMiCalculatorClassName);
		ConditionalMutualInfoCalculatorMultiVariate condMiCalc = condMiClass.newInstance();
		construct(condMiCalc);
	}

	/**
	 * Construct a conditional transfer entropy calculator using an instance of
	 * condMiCalcClass as the underlying conditional mutual information calculator.
	 * 
	 * @param condMiCalcClass the class which must implement
	 * 	{@link ConditionalMutualInfoCalculatorMultiVariate}
	 * @throws InstantiationException if the given class cannot be instantiated
	 * @throws IllegalAccessException if illegal access occurs while trying to create an instance
	 *   of the class
	 * @throws ClassNotFoundException if the given class is not found
	 */
	public ConditionalTransferEntropyCalculatorViaCondMutualInfo(Class<ConditionalMutualInfoCalculatorMultiVariate> condMiCalcClass)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		ConditionalMutualInfoCalculatorMultiVariate condMiCalc = condMiCalcClass.newInstance();
		construct(condMiCalc);
	}

	/**
	 * Construct this calculator by passing in a constructed but not initialised
	 * underlying Conditional Mutual information calculator.
	 * 
	 * @param condMiCalc An instantiated conditional mutual information calculator.
	 * @throws Exception if the supplied calculator has not yet been instantiated.
	 */
	public ConditionalTransferEntropyCalculatorViaCondMutualInfo(ConditionalMutualInfoCalculatorMultiVariate condMiCalc) throws Exception {
		if (condMiCalc == null) {
			throw new Exception("Conditional MI calculator used to construct ConditionalTransferEntropyCalculatorViaCondMutualInfo " +
					" must have already been instantiated.");
		}
		construct(condMiCalc);
	}
	
	/**
	 * Internal method to set the conditional mutual information calculator.
	 * Can be overridden if anything else needs to be done with it by the child classes.
	 * 
	 * @param condMiCalc
	 */
	protected void construct(ConditionalMutualInfoCalculatorMultiVariate condMiCalc) {
		this.condMiCalc = condMiCalc;
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorCommon#initialise()
	 */
	public void initialise() throws Exception {
		initialise(k, k_tau, l, l_tau, delay, condEmbedDims, cond_taus, condDelays);
	}
	
	public void initialise(int k) throws Exception {
		initialise(k, k_tau, l, l_tau, delay, condEmbedDims, cond_taus, condDelays);
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#initialise(int, int, int)
	 */
	public void initialise(int k, int l, int condEmbedDim) throws Exception {
		if (condEmbedDim == 0) {
			// No conditional variables:
			initialise(k, 1, l, 1, 1, null, null, null);
		} else {
			// We have a conditional variable:
			int[] condEmbedDimsArray = new int[1];
			condEmbedDimsArray[0] = condEmbedDim;
			int[] cond_taus = new int[1];
			cond_taus[0] = 1;
			int[] cond_delays = new int[1];
			cond_delays[0] = 1;
			initialise(k, 1, l, 1, 1, condEmbedDimsArray, cond_taus, cond_delays);
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#initialise(int, int, int, int, int, int, int, int)
	 */
	public void initialise(int k, int k_tau, int l, int l_tau, int delay,
			int condEmbedDim, int cond_tau, int condDelay) throws Exception {
		if (condEmbedDim == 0) {
			// No conditional variables:
			initialise(k, k_tau, l, l_tau, delay, null, null, null);
		} else {
			// We have a conditional variable:
			int[] condEmbedDimsArray = new int[1];
			condEmbedDimsArray[0] = condEmbedDim;
			int[] cond_taus = new int[1];
			cond_taus[0] = cond_tau;
			int[] cond_delays = new int[1];
			cond_delays[0] = condDelay;
			initialise(k, k_tau, l, l_tau, delay, condEmbedDimsArray, cond_taus, cond_delays);
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#initialise(int, int, int, int, int, int[], int[], int[])
	 */
	public void initialise(int k, int k_tau, int l, int l_tau, int delay,
			int[] condEmbedDims, int[] cond_taus, int[] condDelays)
			throws Exception {
		
		// First, check consistency:
		if (delay < 0) {
			throw new Exception("Cannot compute TE with source-destination delay < 0");
		}
		if (condEmbedDims == null) {
			// Allow this if all conditional parameter arrays null or 0 length
			condEmbedDims = new int[0];
		}
		if (cond_taus == null) {
			// Allow this if all conditional parameter arrays null or 0 length
			cond_taus = new int[0];
		}
		if (condDelays == null) {
			// Allow this if all conditional parameter arrays null or 0 length
			condDelays = new int[0];
		}
		if ((condEmbedDims.length != cond_taus.length) ||
			(condEmbedDims.length != condDelays.length)) {
			throw new Exception("condEmbedDims, cond_taus and condDelays must have" +
					" same length in argument to ConditionalTransferEntropyCalculatorViaCondMutualInfo.initialise()");
		}
		for (int i = 0; i < condDelays.length; i++) {
			if (condDelays[i] < 0) {
				throw new Exception("Cannot compute TE with conditional-destination delay < 0");
			}
		}
		
		// Next, store the parameters.
		this.k = k;
		this.k_tau = k_tau;
		this.l = l;
		this.l_tau = l_tau;
		this.delay = delay;
		this.condEmbedDims = condEmbedDims;
		this.cond_taus = cond_taus;
		this.condDelays = condDelays;
		
		// Now check which point we can start taking observations from in any
		//  addObservations call. These two integers represent the last
		//  point of the destination embedding, in the cases where the destination
		//  embedding itself determines where we can start taking observations, or
		//  the case where the source embedding plus delay is longer and so determines
		//  where we can start taking observations, or the case where
		//  the conditional embeding plus delay is longer and so determines
		//  where we can start taking observations
		int startTimeBasedOnDestPast = (k-1)*k_tau;
		int startTimeBasedOnSourcePast = (l-1)*l_tau + delay - 1;
		int startTimeBasedOnCondPast = 0;
		dimOfConditionals = 0;
		for (int i = 0; i < condDelays.length; i++) {
			// Check what the start time would be based on this conditional variable
			int startTimeBasedOnThisConditional =
					(condEmbedDims[i]-1)*cond_taus[i] + condDelays[i] - 1;
			if (startTimeBasedOnThisConditional > startTimeBasedOnCondPast) {
				startTimeBasedOnCondPast = startTimeBasedOnThisConditional;
			}
			// And while we're looping compute the total dimension of conditionals
			dimOfConditionals += condEmbedDims[i];
		}
		startTimeForFirstDestEmbedding = Math.max(startTimeBasedOnDestPast,
					Math.max(startTimeBasedOnSourcePast, startTimeBasedOnCondPast));

		condMiCalc.initialise(l, 1, k + dimOfConditionals);
	}

	/**
	 * <p>Set the given property to the given value.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * These can include the following properties from {@link ConditionalTransferEntropyCalculator}:
	 * <ul>
	 * 		<li>{@link #K_PROP_NAME}</li>
	 * 		<li>{@link #K_TAU_PROP_NAME}</li>
	 * 		<li>{@link #L_PROP_NAME}</li>
	 * 		<li>{@link #L_TAU_PROP_NAME}</li>
	 * 		<li>{@link #DELAY_PROP_NAME}</li>
	 * 		<li>{@link #COND_EMBED_LENGTHS_PROP_NAME} - supplied as a comma separated integer list</li>
	 * 		<li>{@link #COND_EMBED_DELAYS_PROP_NAME} - supplied as a comma separated integer list</li>
	 * 		<li>{@link #COND_DELAYS_PROP_NAME} - supplied as a comma separated integer list</li>
	 * </ul>
	 * Note that you can pass in any of the last three properties as different length
	 *  arrays to the current values of the others; this will not be checked here (to allow
	 *  them to be changed here sequentially), but will be checked at the next initialisation
	 *  of the object.
	 * </p>
	 * 
	 * <p>Otherwise, it is assumed the property
	 *  is for the underlying {@link ConditionalMutualInfoCalculatorMultiVariate#setProperty(String, String)} implementation.
	 * </p>
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
		} else if (propertyName.equalsIgnoreCase(COND_EMBED_LENGTHS_PROP_NAME)) {
			condEmbedDims = ParsedProperties.parseStringArrayOfInts(propertyValue);
		} else if (propertyName.equalsIgnoreCase(COND_EMBED_DELAYS_PROP_NAME)) {
			cond_taus = ParsedProperties.parseStringArrayOfInts(propertyValue);
		} else if (propertyName.equalsIgnoreCase(COND_DELAYS_PROP_NAME)) {
			condDelays = ParsedProperties.parseStringArrayOfInts(propertyValue);
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
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#setObservations(double[], double[], double[][])
	 */
	public void setObservations(double[] source, double[] destination,
			double[][] conditionals) throws Exception {
		startAddObservations();
		addObservations(source, destination, conditionals);
		finaliseAddObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#setObservations(double[], double[], double[])
	 */
	public void setObservations(double[] source, double[] destination,
			double[] conditionals) throws Exception {
		if (condEmbedDims.length != 1) {
			throw new Exception("Cannot call setObservations(double[], double[], double[]) when the " +
					"conditional TE calculator was not initialised for one conditional variable");
		}
		startAddObservations();
		addObservations(source, destination, conditionals);
		finaliseAddObservations();
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#startAddObservations()
	 */
	public void startAddObservations() {
		condMiCalc.startAddObservations();
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#addObservations(double[], double[], double[][])
	 */
	public void addObservations(double[] source, double[] destination,
			double[][] conditionals) throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (conditionals == null) {
			if (condEmbedDims.length > 0) {
				throw new Exception(String.format("No conditionals supplied (expected %d-dimensional conditionals", condEmbedDims.length));
			} else {
				// This is allowed; make a dummy set of conditionals
				conditionals = new double[destination.length][0];
			}
		}
		if (conditionals.length != destination.length) {
			throw new Exception(String.format("Conditionals and destination lengths (%d and %d) must match!",
					conditionals.length, destination.length));
		}
		// Postcondition -- all time series have same length
		if (source.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to add here, the time series is too short
			// Don't throw an exception, do nothing since more observations
			//  can be added later.
			return;
		}
		if (conditionals[0].length != condEmbedDims.length) {
			throw new Exception(String.format("Number of conditional variables %d does not " +
					"match the initialised number %d", conditionals[0].length, condEmbedDims.length));
		}
		// All parameters are as expected
		double[][][] embeddedVectorsForCondMI =
				embedSourceDestAndConditionalsForCondMI(source, destination, conditionals);
		
		condMiCalc.addObservations(embeddedVectorsForCondMI[0],
				embeddedVectorsForCondMI[1], embeddedVectorsForCondMI[2]);
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#addObservations(double[], double[], double[])
	 */
	public void addObservations(double[] source, double[] destination,
			double[] conditionals) throws Exception {
		if (condEmbedDims.length != 1) {
			throw new Exception("Cannot call addObservations(double[], double[], double[]) when the " +
					"conditional TE calculator was not initialised for one conditional variable");
		}
		double[][] conditionalsIn2D = null;
		if (conditionals != null) {
			// This isn't incredibly efficient, but is easy to code and doesn't cost more
			//  than an increase in the linear time multiplier.
			conditionalsIn2D = new double[conditionals.length][1];
			MatrixUtils.copyIntoColumn(conditionalsIn2D, 0, conditionals);
		}
		addObservations(source, destination, conditionalsIn2D);
	}

	/**
	 * Internal method to take (pre-screened) vectors for a source, destination
	 *  and conditional variables, and embed them using the given embedding 
	 *  parameters, as well as combining the destination past and conditionals,
	 *  making all ready for a conditional MI calculation.
	 * 
	 * @param source source time series observations
	 * @param destination destination time series observations
	 * @param conditionals 2D array of conditional time series observations.
	 * @return double[][][] returnValue; where returnValue[0] is the embedded
	 *  source vectors, returnValue[1] is the destination next values,
	 *  and returnValue[2] is the joined embedded destination past
	 *  and conditionals.
	 * @throws Exception
	 */
	protected double[][][] embedSourceDestAndConditionalsForCondMI(double[] source, double[] destination,
			double[][] conditionals) throws Exception {
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
		// Now combine the destination past vectors with the conditionals:
		double[][] currentCombinedConditionalVectors =
				new double[currentSourcePastVectors.length][k + dimOfConditionals];
		MatrixUtils.arrayCopy(currentDestPastVectors, 0, 0,
				currentCombinedConditionalVectors, 0, 0,
				currentDestPastVectors.length, k);
		int nextColumnToCopyInto = k;
		for (int i = 0; i < condEmbedDims.length; i++) {
			// Extract the embedding for conditional variable i
			double[][] currentThisConditonalVectors = 
					MatrixUtils.makeDelayEmbeddingVector(conditionals, i,
							condEmbedDims[i], this.cond_taus[i],
							startTimeForFirstDestEmbedding + 1 - condDelays[i],
							conditionals.length - startTimeForFirstDestEmbedding - 1);
			// And add this embedding to our set of conditional variables
			MatrixUtils.arrayCopy(currentThisConditonalVectors, 0, 0,
					currentCombinedConditionalVectors, 0, nextColumnToCopyInto,
					currentThisConditonalVectors.length, condEmbedDims[i]);
			nextColumnToCopyInto += condEmbedDims[i];
		}
		
		double[][][] returnSet = new double[3][][];
		returnSet[0] = currentSourcePastVectors;
		returnSet[1] = currentDestNextVectors;
		returnSet[2] = currentCombinedConditionalVectors;
		return returnSet;
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#addObservations(double[], double[], double[][], int, int)
	 */
	public void addObservations(double[] source, double[] destination,
			double[][] conditionals, int startTime, int numTimeSteps)
			throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (conditionals == null) {
			if (condEmbedDims.length > 0) {
				throw new Exception(String.format("No conditionals supplied (expected %d-dimensional conditionals", condEmbedDims.length));
			} else {
				// This is allowed; make a dummy set of conditionals
				conditionals = new double[destination.length][0];
			}
		}
		if (conditionals.length != destination.length) {
			throw new Exception(String.format("Conditionals and destination lengths (%d and %d) must match!",
					conditionals.length, destination.length));
		}
		// Postcondition -- all time series have same length
		if (source.length < startTime + numTimeSteps) {
			// There are not enough observations given the arguments here
			throw new Exception("Not enough observations to set here given startTime and numTimeSteps parameters");
		}
		addObservations(MatrixUtils.select(source, startTime, numTimeSteps),
					    MatrixUtils.select(destination, startTime, numTimeSteps),
					    MatrixUtils.selectRows(conditionals, startTime, numTimeSteps));
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#addObservations(double[], double[], double[], int, int)
	 */
	public void addObservations(double[] source, double[] destination,
			double[] conditionals, int startTime, int numTimeSteps) throws Exception {
		if (condEmbedDims.length != 1) {
			throw new Exception("Cannot call addObservations(double[], double[], double[], int, int) when the " +
					"conditional TE calculator was not initialised for one conditional variable");
		}
		double[][] conditionalsIn2D = null;
		if (conditionals != null) {
			// This isn't incredibly efficient, but is easy to code and doesn't cost more
			//  than an increase in the linear time multiplier.
			conditionalsIn2D = new double[conditionals.length][1];
			MatrixUtils.copyIntoColumn(conditionalsIn2D, 0, conditionals);
		}
		addObservations(source, destination, conditionalsIn2D, startTime, numTimeSteps);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#finaliseAddObservations()
	 */
	public void finaliseAddObservations() throws Exception {
		condMiCalc.finaliseAddObservations();
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#setObservations(double[], double[], double[][], boolean[], boolean[], boolean[][])
	 */
	public void setObservations(double[] source, double[] destination,
			double[][] conditionals, boolean[] sourceValid,
			boolean[] destValid, boolean[][] conditionalsValid)
			throws Exception {
		
		Vector<int[]> startAndEndTimePairs = 
				computeStartAndEndTimePairs(sourceValid, destValid, conditionalsValid);
		
		// We've found the set of start and end times for this pair
		startAddObservations();
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservations(source, destination, conditionals, startTime, endTime - startTime + 1);
		}
		finaliseAddObservations();
	}

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid tuples of source, destinations and conditionals.
	 *   (i.e. all points within the
	 *  embedding vectors must be valid, even if the invalid points won't be included
	 *  in any tuples)
	 * 
	 * Made public so it can be used if one wants to compute the number of
	 *  observations prior to setting the observations.
	 * 
	 * @param sourceValid
	 * @param destValid
	 * @return
	 * @throws Exception 
	 */
	public Vector<int[]> computeStartAndEndTimePairs(boolean[] sourceValid,
			boolean[] destValid, boolean[][] condValid) throws Exception {
		
		if (sourceValid.length != destValid.length) {
			throw new Exception("Validity arrays must be of same length");
		}
		if (condValid.length != destValid.length) {
			throw new Exception("Validity arrays must be of same length");
		}
		
		int lengthOfDestPastRequired = (k-1)*k_tau + 1;
		int lengthOfSourcePastRequired = (l-1)*l_tau + 1;
		int[] lengthOfConditionalsPastsRequired = new int[condEmbedDims.length];
		for (int i = 0; i < condEmbedDims.length; i++) {
			lengthOfConditionalsPastsRequired[i] = (condEmbedDims[i]-1)*cond_taus[i] + 1;
		}

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
				boolean nextCondsValid = true;
				for (int i = 0; i < condEmbedDims.length; i++) {
					nextCondsValid &= condValid[t + 1 - condDelays[i]][i];
				}
				if (nextCondsValid && destValid[t + 1] && sourceValid[t + 1 - delay]) {
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
			// allOk == true at this point
			for (int tBack = delay - 1; tBack < delay - 1 + lengthOfSourcePastRequired; tBack++) {
				if (!sourceValid[t - tBack]) {
					allOk = false;
					break;
				}
			}
			if (!allOk) {
				continue;
			}
			// allOk == true at this point
			for (int i = 0; i < condEmbedDims.length; i++) {
				for (int tBack = condDelays[i] - 1; tBack < condDelays[i] - 1 + lengthOfConditionalsPastsRequired[i]; tBack++) {
					if (!condValid[t - tBack][i]) {
						allOk = false;
						break;
					}
				}
				if (!allOk) {
					continue;
				}
			}
			// allOk == true at this point
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
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorCommon#computeAverageLocalOfObservations()
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		return condMiCalc.computeAverageLocalOfObservations();
	}

	/**
	 * Returns a time series of local conditional TE values.
	 * Pads the first {@link #startTimeForFirstDestEmbedding} elements with zeros (since local TE is undefined here)
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
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#computeLocalUsingPreviousObservations(double[], double[], double[][])
	 */
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations,
			double[][] newCondObservations) throws Exception {
		
		if (newSourceObservations.length != newDestObservations.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					newSourceObservations.length, newDestObservations.length));
		}
		if (newCondObservations == null) {
			if (condEmbedDims.length > 0) {
				throw new Exception(String.format("No conditionals supplied (expected %d-dimensional conditionals)", condEmbedDims.length));
			} else {
				// This is allowed; make a dummy set of conditionals
				newCondObservations = new double[newDestObservations.length][0];
			}
		}
		if (newCondObservations.length != newDestObservations.length) {
			throw new Exception(String.format("Conditionals and destination lengths (%d and %d) must match!",
					newCondObservations.length, newDestObservations.length));
		}
		// Postcondition -- all time series have same length
		if (newCondObservations[0].length != condEmbedDims.length) {
			throw new Exception(String.format("Number of conditional variables %d does not " +
					"match the initialised number %d", newCondObservations[0].length, condEmbedDims.length));
		}
		if (newDestObservations.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to compute for here
			return new double[newDestObservations.length];
		}
		// All parameters are as expected
		double[][][] embeddedVectorsForCondMI =
				embedSourceDestAndConditionalsForCondMI(newSourceObservations,
						newDestObservations, newCondObservations);

		double[] local = condMiCalc.computeLocalUsingPreviousObservations(
				embeddedVectorsForCondMI[0], embeddedVectorsForCondMI[1],
				embeddedVectorsForCondMI[2]);
		// Pad the front of the array with zeros where local TE isn't defined:
		double[] localsToReturn = new double[local.length + startTimeForFirstDestEmbedding + 1];
		System.arraycopy(local, 0, localsToReturn, startTimeForFirstDestEmbedding + 1, local.length);
		return localsToReturn;
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalTransferEntropyCalculator#computeLocalUsingPreviousObservations(double[], double[], double[])
	 */
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations,
			double[] newCondObservations) throws Exception {
		
		if (condEmbedDims.length != 1) {
			throw new Exception("Cannot call computeLocalUsingPreviousObservations(double[], double[], double[]) when the " +
					"conditional TE calculator was not initialised for one conditional variable");
		}
		double[][] conditionalsIn2D = null;
		if (newCondObservations != null) {
			// This isn't incredibly efficient, but is easy to code and doesn't cost more
			//  than an increase in the linear time multiplier.
			conditionalsIn2D = new double[newCondObservations.length][1];
			MatrixUtils.copyIntoColumn(conditionalsIn2D, 0, newCondObservations);
		}
		return computeLocalUsingPreviousObservations(newSourceObservations,
				newDestObservations, conditionalsIn2D);
	}
	
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		return condMiCalc.computeSignificance(1, numPermutationsToCheck); // Reorder the source vectors
	}

	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		return condMiCalc.computeSignificance(1, newOrderings); // Reorder the source vectors
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
