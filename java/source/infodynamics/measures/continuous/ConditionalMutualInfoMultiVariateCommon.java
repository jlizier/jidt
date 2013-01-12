package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

import java.util.Iterator;
import java.util.Vector;

/**
 * <p>Base class for implementations of {@link ConditionalMutualInfoCalculatorMultiVariate},
 * e.g. kernel estimation, Kraskov style extensions.
 * It implements some common code to be used across conditional 
 * mutual information calculators
 * </p>
 * 
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com
 *
 */
public abstract class ConditionalMutualInfoMultiVariateCommon implements
		ConditionalMutualInfoCalculatorMultiVariate {

	/**
	 * Number of dimenions for each of our multivariate data sets
	 */
	protected int dimensionsVar1 = 1;
	protected int dimensionsVar2 = 1;
	protected int dimensionsCond = 1;
	
	/**
	 * The set of observations for var1, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link #addObservations(double[][], double[][], double[][])} functions.
	 */
	protected double[][] var1Observations;

	/**
	 * The set of observations for var2, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link #addObservations(double[][], double[][], double[][])} functions.
	 */
	protected double[][] var2Observations;

	/**
	 * The set of observations for the conditional, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link #addObservations(double[][], double[][], double[][])} functions.
	 */
	protected double[][] condObservations;

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
	protected boolean condMiComputed;

	/**
	 * Whether to report debug messages or not
	 */
	protected boolean debug;

	/**
	 * Storage for var1 observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfVar1Observations;
	/**
	 * Storage for var2 observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfVar2Observations;
	/**
	 * Storage for conditional variable observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfCondObservations;
	
	protected boolean addedMoreThanOneObservationSet;

	/**
	 * Clear any previously supplied probability distributions and prepare
	 * the calculator to be used again.
	 * 
	 * @param var1Dimensions number of joint variables in variable 1
	 * @param var2Dimensions number of joint variables in variable 2
	 * @param condDimensions number of joint variables in the conditional
	 */
	public void initialise(int var1Dimensions, int var2Dimensions, int condDimensions) {
		dimensionsVar1 = var1Dimensions;
		dimensionsVar2 = var2Dimensions;
		dimensionsCond = condDimensions;
		lastAverage = 0.0;
		totalObservations = 0;
		condMiComputed = false;
		var1Observations = null;
		var2Observations = null;
		condObservations = null;
	}

	// No properties to set on this abstract calculator

	/**
	 * Provide the complete set of observations to use to compute the
	 *  mutual information.
	 * One cannot use the
	 *  {@link #addObservations(double[][], double[][], double[][])}
	 *  style methods after this without calling 
	 *  {@link #initialise(int, int, int)} again first.
	 * 
	 * @param var1 time series of multivariate variable 1 observations (first
	 *  index is time, second is variable number)
	 * @param var2 time series of multivariate variable 2 observations (first
	 *  index is time, second is variable number)
	 * @param cond time series of multivariate conditional variable observations (first
	 *  index is time, second is variable number)
	 */
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond) throws Exception {
		startAddObservations();
		addObservations(var1, var2, cond);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	/**
	 * Elect to add in the observations from several disjoint time series.
	 *
	 */
	public void startAddObservations() {
		vectorOfVar1Observations = new Vector<double[][]>();
		vectorOfVar2Observations = new Vector<double[][]>();
		vectorOfCondObservations = new Vector<double[][]>();
	}
	
	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond) throws Exception {
		if (vectorOfVar1Observations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if ((var1.length != var2.length) || (var1.length != cond.length)) {
			throw new Exception(String.format("Observation vector lengths (%d, %d and %d) must match!",
					var1.length, var2.length, cond.length));
		}
		if (var1[0].length != dimensionsVar1) {
			throw new Exception("Number of joint variables in var1 data " +
					"does not match the initialised value");
		}
		if (var2[0].length != dimensionsVar2) {
			throw new Exception("Number of joint variables in var2 data " +
					"does not match the initialised value");
		}
		if (cond[0].length != dimensionsCond) {
			throw new Exception("Number of joint variables in cond data " +
					"does not match the initialised value");
		}
		vectorOfVar1Observations.add(var1);
		vectorOfVar2Observations.add(var2);
		vectorOfCondObservations.add(cond);
		if (vectorOfVar1Observations.size() > 1) {
			addedMoreThanOneObservationSet = true;
		}
	}

	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond,
			int startTime, int numTimeSteps) throws Exception {
		if (vectorOfVar1Observations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		double[][] var1ToAdd = new double[numTimeSteps][];
		System.arraycopy(var1, startTime, var1ToAdd, 0, numTimeSteps);
		double[][] var2ToAdd = new double[numTimeSteps][];
		System.arraycopy(var2, startTime, var2ToAdd, 0, numTimeSteps);
		double[][] condToAdd = new double[numTimeSteps][];
		System.arraycopy(cond, startTime, condToAdd, 0, numTimeSteps);
		addObservations(var1ToAdd, var2ToAdd, condToAdd);
	}

	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[] var1Valid, boolean[] var2Valid,
			boolean[] condValid) throws Exception {
		
		Vector<int[]> startAndEndTimePairs =
				computeStartAndEndTimePairs(var1Valid, var2Valid, condValid);
		
		// We've found the set of start and end times for this pair
		startAddObservations();
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservations(var1, var2, cond, startTime, endTime - startTime + 1);
		}
		finaliseAddObservations();
	}

	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[][] var1Valid, boolean[][] var2Valid,
			boolean[][] condValid) throws Exception {

		boolean[] allVar1Valid = MatrixUtils.andRows(var1Valid);
		boolean[] allVar2Valid = MatrixUtils.andRows(var2Valid);
		boolean[] allCondValid = MatrixUtils.andRows(condValid);
		setObservations(var1, var2, cond, allVar1Valid, allVar2Valid, allCondValid);
	}

	/**
	 * Finalise the addition of multiple observation sets.
	 * 
	 * This default implementation simply puts all of the observations into
	 *  the {@link #var1Observations}, {@link #var2Observations} 
	 *  and {@link #condObservations} arrays.
	 * Usually child implementations will override this, call this implementation
	 *  to perform the common processing, then perform their own processing.
	 *  
	 * @throws Exception Allow child classes to throw an exception if there
	 *  is an issue detected specific to that calculator.
	 * 
	 */
	public void finaliseAddObservations() throws Exception {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[][] var2 : vectorOfVar2Observations) {
			totalObservations += var2.length;
		}
		var1Observations = new double[totalObservations][dimensionsVar1];
		var2Observations = new double[totalObservations][dimensionsVar2];
		condObservations = new double[totalObservations][dimensionsCond];
		
		int startObservation = 0;
		Iterator<double[][]> iteratorVar2 = vectorOfVar2Observations.iterator();
		Iterator<double[][]> iteratorCond = vectorOfCondObservations.iterator();
		for (double[][] var1 : vectorOfVar1Observations) {
			double[][] var2 = iteratorVar2.next();
			double[][] cond = iteratorCond.next();
			// Copy the data from these given observations into our master 
			//  array, aligning them incorporating the timeDiff:
			MatrixUtils.arrayCopy(var1, 0, 0,
					var1Observations, startObservation, 0,
					var1.length, dimensionsVar1);
			MatrixUtils.arrayCopy(var2, 0, 0,
					var2Observations, startObservation, 0,
					var2.length, dimensionsVar2);
			MatrixUtils.arrayCopy(cond, 0, 0,
					condObservations, startObservation, 0,
					cond.length, dimensionsCond);
			startObservation += var2.length;
		}
		// We don't need to keep the vectors of observation sets anymore:
		vectorOfVar1Observations = null;
		vectorOfVar2Observations = null;
		vectorOfCondObservations = null;
	}
	
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// Use var1 length (all variables have same length) even though
		//  we may be randomising the other variable:
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(
				var1Observations.length, numPermutationsToCheck);
		return computeSignificance(variableToReorder, newOrderings);
	}

	/**
	 * <p>As per {@link #computeSignificance(int, int)} but supplies
	 *  the re-orderings of the observations of the named variable.</p>
	 * 
	 * <p>We provide a simple implementation which would be suitable for
	 *  any child class, though the child class may prefer to make its
	 *  own implementation to make class-specific optimisations.
	 *  Child classes must implement {@link java.lang.Cloneable}
	 *  for this method to be callable for them, and indeed implement
	 *  the clone() method in a way that protects their structure
	 *  from alteration by surrogate data being supplied to it.</p>
	 * 
	 * @param variableToReorder 1 for variable 1, 2 for variable 2
	 * @param newOrderings first index is permutation number, i.e. newOrderings[i]
	 * 		is an array of 1 permutation of 0..n-1, where there were n observations.
	 * 		If the length of each permutation in newOrderings
	 * 		is not equal to numObservations, an Exception is thrown. 
	 * @param newOrderings the specific new orderings to use
	 * @return
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		if (!condMiComputed) {
			computeAverageLocalOfObservations();
		}
		
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays - child classes should override this)
		ConditionalMutualInfoMultiVariateCommon miSurrogateCalculator =
				(ConditionalMutualInfoMultiVariateCommon) this.clone();
		
		double[] surrogateMeasurements = new double[numPermutationsToCheck];
		
		// Now compute the MI for each set of shuffled data:
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Generate a new re-ordered source data
			double[][] shuffledData = 
					MatrixUtils.extractSelectedTimePointsReusingArrays(
							(variableToReorder == 1) ? var1Observations : var2Observations,
							newOrderings[i]);
			// Perform new initialisations
			miSurrogateCalculator.initialise(
					dimensionsVar1, dimensionsVar2, dimensionsCond);
			// Set new observations
			if (variableToReorder == 1) {
				miSurrogateCalculator.setObservations(shuffledData,
						var2Observations, condObservations);
			} else {
				miSurrogateCalculator.setObservations(var1Observations,
						shuffledData, condObservations);
			}
			// Compute the MI
			surrogateMeasurements[i] = miSurrogateCalculator.computeAverageLocalOfObservations();
			if (debug){
				System.out.println("New MI was " + surrogateMeasurements[i]);
			}
		}
		
		return new EmpiricalMeasurementDistribution(surrogateMeasurements, lastAverage);
	}

	/**
	 * <p>Compute the mutual information if the given variable were
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
	 * @param variableToReorder 1 for variable 1, 2 for variable 2
	 * @param newOrdering array of time indices with which to reorder the data
	 * @return a surrogate MI evaluated for the given ordering of the source variable
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int variableToReorder, int[] newOrdering)
			throws Exception {
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays - child class should override this)
		ConditionalMutualInfoMultiVariateCommon miSurrogateCalculator =
				(ConditionalMutualInfoMultiVariateCommon) this.clone();
		
		// Generate a new re-ordered source data
		double[][] shuffledData =
				MatrixUtils.extractSelectedTimePointsReusingArrays(
					(variableToReorder == 1) ? var1Observations : var2Observations,
					newOrdering);
		// Perform new initialisations
		miSurrogateCalculator.initialise(
				dimensionsVar1, dimensionsVar2, dimensionsCond);
		// Set new observations
		if (variableToReorder == 1) {
			miSurrogateCalculator.setObservations(shuffledData,
					var2Observations, condObservations);
		} else {
			miSurrogateCalculator.setObservations(var1Observations,
					shuffledData, condObservations);
		}
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
	 *  valid series of all variables.
	 * 
	 * Made public so it can be used if one wants to compute the number of
	 *  observations prior to setting the observations.
	 * 
	 * @param var1Valid
	 * @param var2Valid
	 * @return
	 */
	public Vector<int[]> computeStartAndEndTimePairs(
			boolean[] var1Valid, boolean[] var2Valid, boolean[] var3Valid) {
		// Scan along the data avoiding invalid values
		int startTime = 0;
		int endTime = 0;
		boolean lookingForStart = true;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();
		for (int t = 0; t < var2Valid.length; t++) {
			if (lookingForStart) {
				// Precondition: startTime holds a candidate start time
				//  (var1 value is at startTime == t)
				if (var1Valid[t] && var2Valid[t] && var3Valid[t]) {
					// This point is OK at the variables
					// Set a candidate endTime
					endTime = t;
					lookingForStart = false;
					if (t == var1Valid.length - 1) {
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
				if (var1Valid[t] && var2Valid[t] && var3Valid[t]) {
					// We can extend
					endTime = t;
				} else {
					terminateSequence = true;
				}
				if (t == var2Valid.length - 1) {
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
