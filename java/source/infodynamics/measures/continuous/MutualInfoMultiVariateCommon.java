package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

import java.util.Iterator;
import java.util.Vector;

/**
 * <p>Base class for implementations of {@link MutualInfoCalculatorMultiVariate},
 * e.g. kernel estimation, Kraskov style extensions.
 * It implements some common code to be used across mutual information calculators
 * </p>
 * 
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com
 *
 */
public abstract class MutualInfoMultiVariateCommon implements
		MutualInfoCalculatorMultiVariate {

	/**
	 * Number of dimenions for each of our multivariate data sets
	 */
	protected int dimensionsDest = 1;
	protected int dimensionsSource = 1;
	
	/**
	 * The set of source observations, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link addObservations(double[][], double[][])} functions.
	 */
	protected double[][] sourceObservations;

	/**
	 * The set of destination observations, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link addObservations(double[][], double[][])} functions.
	 */
	protected double[][] destObservations;

	/**
	 * Total number of observations supplied.
	 * Only valid after {@link #finaliseAddObservations()} is called.
	 */
	protected int totalObservations = 0;

	/**
	 * Store the last computed average
	 */
	protected double lastAverage;

	/**
	 * Track whether we've computed the average for the supplied
	 *  observations yet
	 */
	protected boolean miComputed;

	/**
	 * Whether to report debug messages or not
	 */
	protected boolean debug;

	/**
	 * Time difference from the source to the destination observations
	 *  (i.e. destination lags the source by this time:
	 *  	we compute I(x_{1,n}; x_{2,n+timeDiff})
	 * (Note that our internal sourceObservations and destObservations
	 *  are adjusted so that there is no timeDiff between them).
	 */
	protected int timeDiff;

	/**
	 * Storage for source observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfSourceObservations;
	/**
	 * Storage for destination observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfDestinationObservations;
	
	protected boolean addedMoreThanOneObservationSet;

	/**
	 * Clear any previously supplied probability distributions and prepare
	 * the calculator to be used again.
	 * 
	 * @param sourceDimensions number of joint variables in the source
	 * @param destDimensions number of joint variables in the destination
	 */
	public void initialise(int sourceDimensions, int destDimensions) {
		dimensionsSource = sourceDimensions;
		dimensionsDest = destDimensions;
		lastAverage = 0.0;
		totalObservations = 0;
		miComputed = false;
		sourceObservations = null;
		destObservations = null;
	}

	/**
	 * <p>Set the given property to the given value.
	 * These can include:
	 * <ul>
	 * 		<li>{@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}</li>
	 * </ul>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_TIME_DIFF)) {
			timeDiff = Integer.parseInt(propertyValue);
			if (timeDiff < 0) {
				throw new Exception("Time difference must be >= 0. Flip data1 and data2 around if required.");
			}
		} else {
			// No property was set here
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/**
	 * Provide the complete set of observations to use to compute the
	 *  mutual information.
	 * One cannot use the {@link #addObservations(double[][], double[][])}
	 *  style methods after this without calling 
	 *  {@link #initialise(int, int)} again first.
	 * 
	 * @param source time series of multivariate source observations (first
	 *  index is time, second is variable number)
	 * @param destination time series of multivariate destination observations (first
	 *  index is time, second is variable number)
	 */
	public void setObservations(double[][] source, double[][] destination) throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	/**
	 * Elect to add in the observations from several disjoint time series.
	 *
	 */
	public void startAddObservations() {
		vectorOfSourceObservations = new Vector<double[][]>();
		vectorOfDestinationObservations = new Vector<double[][]>();
	}
	
	/**
	 * Add some more observations.
	 * Note that the arrays source and destination must not be over-written by the user
	 *  until after finaliseAddObservations() has been called.
	 * 
	 * @param source
	 * @param destination
	 */
	public void addObservations(double[][] source, double[][] destination) throws Exception {
		if (vectorOfSourceObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length <= timeDiff) {
			// we won't be taking any observations here
			return;
		}
		if (source[0].length != dimensionsSource) {
			throw new Exception("Number of joint variables in source data " +
					"does not match the initialised value");
		}
		if (destination[0].length != dimensionsDest) {
			throw new Exception("Number of joint variables in destination data " +
					"does not match the initialised value");
		}
		vectorOfSourceObservations.add(source);
		vectorOfDestinationObservations.add(destination);
		if (vectorOfSourceObservations.size() > 1) {
			addedMoreThanOneObservationSet = true;
		}
	}

	/**
	 * Add some more observations but only for the selected time steps.
	 * 
	 * @param source
	 * @param destination
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 */
	public void addObservations(double[][] source, double[][] destination,
			int startTime, int numTimeSteps) throws Exception {
		if (vectorOfSourceObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (numTimeSteps <= timeDiff) {
			// We won't be taking any observations here
			return;
		}
		double[][] sourceToAdd = new double[numTimeSteps][];
		System.arraycopy(source, startTime, sourceToAdd, 0, numTimeSteps);
		vectorOfSourceObservations.add(sourceToAdd);
		double[][] destToAdd = new double[numTimeSteps][];
		System.arraycopy(destination, startTime, destToAdd, 0, numTimeSteps);
		vectorOfDestinationObservations.add(destToAdd);
		if (vectorOfSourceObservations.size() > 1) {
			addedMoreThanOneObservationSet = true;
		}
	}

	/**
	 * Set the observations, but only where both the source
	 *  and destination time steps are valid.
	 * 
	 * 
	 * @param source
	 * @param destination
	 * @param sourceValid whether the values for the source are valid at
	 *  each time step
	 * @param destValid whether the values for the destination are valid at
	 *  each time step
	 */
	public void setObservations(double[][] source, double[][] destination,
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
	 * Set the observations, but only where all variables for both the source
	 *  and destination time steps are valid.
	 * 
	 * 
	 * @param source
	 * @param destination
	 * @param sourceValid whether each variable for the source are valid at
	 *  each time step
	 * @param destValid whether each variable for the destination are valid at
	 *  each time step
	 */
	public void setObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception {

		boolean[] allSourceValid = MatrixUtils.andRows(sourceValid);
		boolean[] allDestValid = MatrixUtils.andRows(destValid);
		setObservations(source, destination, allSourceValid, allDestValid);
	}

	/**
	 * Finalise the addition of multiple observation sets.
	 * 
	 * This default implementation simply puts all of the observations into
	 *  the {@link #sourceObservations} and {@link #destObservations} arrays.
	 * Usually child implementations will override this, call this implementation
	 *  to perform the common processing, then perform their own processing.
	 * @throws Exception Allow child classes to throw an exception if there
	 *  is an issue detected specific to that calculator.
	 * 
	 */
	public void finaliseAddObservations() throws Exception {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[][] destination : vectorOfDestinationObservations) {
			totalObservations += destination.length - timeDiff;
		}
		destObservations = new double[totalObservations][dimensionsDest];
		sourceObservations = new double[totalObservations][dimensionsSource];
		
		// Construct the joint vectors from the given observations
		//  (removing redundant data which is outside any timeDiff)
		int startObservation = 0;
		Iterator<double[][]> iterator = vectorOfDestinationObservations.iterator();
		for (double[][] source : vectorOfSourceObservations) {
			double[][] destination = iterator.next();
			// Copy the data from these given observations into our master 
			//  array, aligning them incorporating the timeDiff:
			MatrixUtils.arrayCopy(source, 0, 0,
					sourceObservations, startObservation, 0,
					source.length - timeDiff, dimensionsSource);
			MatrixUtils.arrayCopy(destination, timeDiff, 0,
					destObservations, startObservation, 0,
					destination.length - timeDiff, dimensionsDest);
			startObservation += destination.length - timeDiff;
		}
		// We don't need to keep the vectors of observation sets anymore:
		vectorOfSourceObservations = null;
		vectorOfDestinationObservations = null;
	}
	
	/**
	 * <p>Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 * This is in the spirit of Chavez et al. (see below), 
	 *  which was performed for Transfer entropy.</p>
	 * 
	 * <p>We permute the source variable against the destination
	 *  (though in theory this doesn't matter for this function call).
	 * </p>
	 * 
	 * @param numPermutationsToCheck
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 * @see "Chavez et. al., 'Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals', Journal of Neuroscience Methods 124 (2003) 113-128"
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(
				sourceObservations.length, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}

	/**
	 * <p>Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 * This is in the spirit of Chavez et al. (see below), 
	 *  which was performed for Transfer entropy.</p>
	 * 
	 * <p>We provide a simple implementation which would be suitable for
	 *  any child class, though the child class may prefer to make its
	 *  own implementation to make class-specific optimisations.
	 *  Child classes must implement {@link java.lang.Cloneable}
	 *  for this method to be callable for them, and indeed implement
	 *  the clone() method in a way that protects their structure
	 *  from alteration by surrogate data being supplied to it.</p>
	 *  
	 * <p>We permute the source variable against the destination
	 *  to be consistent with the description in {@link ChannelCalculator#computeSignificance(int[][])}
	 *  (though in theory this doesn't matter for this function call).
	 * </p>
	 * 
	 * @param newOrderings the specific new orderings to use
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 * @see "Chavez et. al., 'Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals', Journal of Neuroscience Methods 124 (2003) 113-128"
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays)
		MutualInfoMultiVariateCommon miSurrogateCalculator =
				(MutualInfoMultiVariateCommon) this.clone();
		
		double[] surrogateMeasurements = new double[numPermutationsToCheck];
		
		// Now compute the MI for each set of shuffled data:
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Generate a new re-ordered source data
			double[][] shuffledSourceData = 
					MatrixUtils.extractSelectedTimePointsReusingArrays(
							sourceObservations, newOrderings[i]);
			// Perform new initialisations
			miSurrogateCalculator.initialise(dimensionsSource, dimensionsDest);
			// Set new observations
			miSurrogateCalculator.setObservations(shuffledSourceData, destObservations);
			// Compute the MI
			surrogateMeasurements[i] = miSurrogateCalculator.computeAverageLocalOfObservations();
			if (debug){
				System.out.println("New MI was " + surrogateMeasurements[i]);
			}
		}
		
		return new EmpiricalMeasurementDistribution(surrogateMeasurements, lastAverage);
	}

	/**
	 * <p>Compute the mutual information if the first (source) variable were
	 *  ordered as per the ordering specified in newOrdering.</p>
	 * 
	 * <p>We provide a simple implementation which would be suitable for
	 *  any child class, though the child class may prefer to make its
	 *  own implementation to make class-specific optimisations.
	 *  Child classes must implement {@link java.lang.Cloneable}
	 *  for this method to be callable for them, and indeed implement
	 *  the clone() method in a way that protects their structure
	 *  from alteration by surrogate data being supplied to it.</p>
	 *  
	 * @param newOrdering array of time indices with which to reorder the data
	 * @return a surrogate MI evaluated for the given ordering of the source variable
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[] newOrdering)
			throws Exception {
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays)
		MutualInfoMultiVariateCommon miSurrogateCalculator =
				(MutualInfoMultiVariateCommon) this.clone();
		
		// Generate a new re-ordered source data
		double[][] shuffledSourceData =
				MatrixUtils.extractSelectedTimePointsReusingArrays(
					sourceObservations, newOrdering);
		// Perform new initialisations
		miSurrogateCalculator.initialise(dimensionsSource, dimensionsDest);
		// Set new observations
		miSurrogateCalculator.setObservations(shuffledSourceData, destObservations);
		// Compute the MI
		return miSurrogateCalculator.computeAverageLocalOfObservations();
	}

	/**
	 * Set whether debug messages will be displayed
	 * 
	 * @param debug debug setting
	 */
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	/**
	 * @return the previously computed average mutual information
	 */
	public double getLastAverage() {
		return lastAverage;
	}

	/**
	 * @return the number of supplied observations
	 * @throws Exception if child class computes MI without explicit observations
	 */
	public int getNumObservations() throws Exception {
		return totalObservations;
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
	 */
	public Vector<int[]> computeStartAndEndTimePairs(boolean[] sourceValid, boolean[] destValid) {
		// Scan along the data avoiding invalid values
		int startTime = 0;
		int endTime = 0;
		boolean lookingForStart = true;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();
		for (int t = 0; t < destValid.length; t++) {
			if (lookingForStart) {
				// Precondition: startTime holds a candidate start time
				//  (source value is at startTime == t - timeDiff)
				if (destValid[t] && sourceValid[t - timeDiff]) {
					// This point is OK at the source and destination
					// Set a candidate endTime
					endTime = t;
					lookingForStart = false;
					if (t == destValid.length - 1) {
						// we need to terminate now
						int[] timePair = new int[2];
						timePair[0] = startTime;
						timePair[1] = endTime;
						startAndEndTimePairs.add(timePair);
						// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
					}
				} else {
					// We need to keep looking.
					// Move the potential start time to the next point
					startTime++;
				}
			} else {
				// Precondition: startTime holds the start time for this set, 
				//  endTime holds a candidate end time
				// Check if we can include the current time step
				boolean terminateSequence = false;
				if (destValid[t] && sourceValid[t - timeDiff]) {
					// We can extend
					endTime = t;
				} else {
					terminateSequence = true;
				}
				if (t == destValid.length - 1) {
					// we need to terminate the sequence anyway
					terminateSequence = true;
				}
				if (terminateSequence) {
					// This section is done
					int[] timePair = new int[2];
					timePair[0] = startTime;
					timePair[1] = endTime;
					startAndEndTimePairs.add(timePair);
					// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
					lookingForStart = true;
					startTime = t + 1;
				}
			}
		}
		return startAndEndTimePairs;
	}	
}
