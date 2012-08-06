package infodynamics.measures.continuous;

import infodynamics.utils.MatrixUtils;

import java.util.Iterator;
import java.util.Vector;

/**
 * <p>Base class for implementations of the mutual information,
 * e.g. kernel estimation, Kraskov style extensions.
 * It implements some common code to be used across mutual information calculators
 * </p>
 * 
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com
 *
 * @see T. M. Cover and J. A. Thomas, "Elements of Information
Theory" (John Wiley & Sons, New York, 1991).
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

	protected double lastAverage;

	protected boolean debug;

	/**
	 * Time difference from the source to the destination
	 *  (i.e. destination lags the source by this time:
	 *  	we compute I(x_{1,n}; x_{2,n+timeDiff})
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
		lastAverage = 0;
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
	 * @throws Exception
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		if (propertyName.equalsIgnoreCase(PROP_TIME_DIFF)) {
			timeDiff = Integer.parseInt(propertyValue);
			if (timeDiff < 0) {
				throw new RuntimeException("Time difference must be >= 0. Flip data1 and data2 around if required.");
			}
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
	 * Add some more observations.
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
	 * Finalise the addition of multiple observation sets.
	 * 
	 * This default implementation simply puts all of the observations into
	 *  the {@link #sourceObservations} and {@link #destObservations} arrays.
	 * Usually child implementations will override this, call this implementation
	 *  to perform the common processing, then perform their own processing.
	 * 
	 */
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		int totalObservations = 0;
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
			//  array
			MatrixUtils.arrayCopy(source, 0, 0,
					sourceObservations, startObservation, 0,
					source.length - timeDiff, dimensionsSource);
			MatrixUtils.arrayCopy(destination, timeDiff, 0,
					destObservations, startObservation, 0,
					destination.length - timeDiff, dimensionsDest);
			startObservation += destination.length - timeDiff;
		}
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

}
