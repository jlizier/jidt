package infodynamics.measures.continuous;

import infodynamics.utils.MatrixUtils;

import java.util.Vector;


/**
 * <p>Base class for implementations of the transfer entropy (see Schreiber, PRL, 2000)
 *  and local transfer entropy (see Lizier et al, PRE, 2008).
 *  We use the term <i>apparent</i> transfer entropy to mean that
 *  we compute the transfer that appears to come from a single
 *  source variable, without examining any other potential sources
 *  (see Lizier et al, PRE, 2008).</p>
 * 
 * <p>Specifically, this provides base implementations of the transfer entropy for 
 * <i>continuous</i>-valued variables (univariate).
 * It provides common code used in multiple child implementations.</p>
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
public abstract class TransferEntropyCommon implements
		TransferEntropyCalculator {

	protected final static double log2 = Math.log(2.0);
	
	/**
	 * Length of past history to consider
	 */
	protected int k;
	protected int totalObservations = 0;
	protected boolean debug = false;
	protected double lastAverage;

	/**
	 * Storage for source observations for addObservsations
	 */
	protected Vector<double[]> vectorOfSourceObservations;
	/**
	 * Storage for destination observations for addObservsations
	 */
	protected Vector<double[]> vectorOfDestinationObservations;
	
	protected boolean addedMoreThanOneObservationSet;

	/**
	 * Initialise the calculator, and call initialise() to complete
	 * 
	 * @param k Length of past history to consider
	 */
	public void initialise(int k) throws Exception {
		this.k = k;
		addedMoreThanOneObservationSet = false;
		initialise();
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
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}

	/**
	 * Elect to add in the observations from several disjoint time series.
	 *
	 */
	public void startAddObservations() {
		vectorOfSourceObservations = new Vector<double[]>();
		vectorOfDestinationObservations = new Vector<double[]>();
	}
	
	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. this
	 *  will not join them up to examine k
	 *  consecutive values in time (which would be an incorrect transfer entropy
	 *  calculation, since the end of one observation set should not 
	 *  necessarily be followed by the start of another).</p>
	 *  
	 * <p>Note that the arrays source and destination must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @throws Exception
	 */
	public void addObservations(double[] source, double[] destination) throws Exception {
		if (vectorOfSourceObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length <= k) {
			// we won't be taking any observations here
			return;
		}
		vectorOfSourceObservations.add(source);
		vectorOfDestinationObservations.add(destination);
	}

	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. this
	 *  will not join them up to examine k
	 *  consecutive values in time (which would be an incorrect transfer entropy
	 *  calculation, since the end of one observation set should not 
	 *  necessarily be followed by the start of another).</p>
	 *  
	 * <p>Note that the arrays source and destination must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 * @throws Exception
	 */
	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception {
		if (vectorOfSourceObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (numTimeSteps <= k) {
			// We won't be taking any observations here
			return;
		}
		double[] sourceToAdd = new double[numTimeSteps];
		System.arraycopy(source, startTime, sourceToAdd, 0, numTimeSteps);
		vectorOfSourceObservations.add(sourceToAdd);
		double[] destToAdd = new double[numTimeSteps];
		System.arraycopy(destination, startTime, destToAdd, 0, numTimeSteps);
		vectorOfDestinationObservations.add(destToAdd);
	}

	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param sourceValid time series (with time indices the same as source)
	 *  indicating whether the source at that point is valid.
	 * @param destValid time series (with time indices the same as destination)
	 *  indicating whether the destination at that point is valid.
	 * @throws Exception
	 */
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
				if (destValid[t]) {
					// This point is OK at the destination
					if (t - startTime < k) {
						// We're still checking the past history only, so
						continue;
					} else {
						// We've got the full past history ok, so check the source also
						if (sourceValid[t - 1]) {
							// source is good to go also
							// set a candidate endTime
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
							// source was not good to go, so try moving along one time point
							startTime++;
						}
					}
				} else {
					// We need to keep looking.
					// Move the potential start time to the next point
					startTime = t + 1;
				}
			} else {
				// Precondition: startTime holds the start time for this set, 
				//  endTime holds a candidate end time
				// Check if we can include the current time step
				boolean terminateSequence = false;
				if (destValid[t] && sourceValid[t - 1]) {
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
	
	/**
	 * Generate a vector for each time step, containing the past k states of the destination.
	 * Does not include a vector for the first k time steps.
	 * 
	 * @param destination
	 * @return array of vectors for each time step
	 */
	protected double[][] makeJointVectorForPast(double[] destination) {
		try {
			// We want one less delay vector here - we don't need the last k point,
			//  because there is no next state for these.
			return MatrixUtils.makeDelayEmbeddingVector(destination, k, k-1, destination.length - k);
		} catch (Exception e) {
			// The parameters for the above call should be fine, so we don't expect to
			//  throw an Exception here - embed in a RuntimeException if it occurs 
			throw new RuntimeException(e);
		}
	}
	
	/**
	 * Generate a vector for each time step, containing the past k states of
	 *  the destination, and the current state.
	 * Does not include a vector for the first k time steps.
	 * 
	 * @param destination
	 * @return
	 */
	protected double[][] makeJointVectorForNextPast(double[] destination) {
		// We want all delay vectors here
		return MatrixUtils.makeDelayEmbeddingVector(destination, k+1);
	}
	
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	public double getLastAverage() {
		return lastAverage;
	}

	public int getNumObservations() {
		return totalObservations;
	}
	
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public boolean getAddedMoreThanOneObservationSet() {
		return addedMoreThanOneObservationSet;
	}
}
