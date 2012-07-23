package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariate;
import infodynamics.utils.MatrixUtils;

import java.util.Iterator;
import java.util.Vector;


/**
 * <p>Compute the Transfer Entropy using the Kraskov estimation method.<p>
 * <p>This calculator extends the {@link TransferEntropyCalculatorKraskovByMulti TransferEntropyCalculatorKraskov}
 *  by allowing multivariate sources and destinations. See 
 *  {@link TransferEntropyCalculatorKraskovByMulti TransferEntropyCalculatorKraskov}
 *  for further comments on the implementation of Transfer entropy via the Kraskov 
 *  Mutual Information estimation.</p>
 * </p>
 * 
 * 
 * @author Joseph Lizier; joseph.lizier at gmail.com
 *
 */
public class TransferEntropyCalculatorMultiVariateKraskovByMulti
	extends TransferEntropyCalculatorKraskovByMulti implements TransferEntropyCalculatorMultiVariate {

	private int destDimensions;
	private int sourceDimensions;

	/**
	 * Storage for source observations for addObservsations
	 */
	private Vector<double[][]> vectorOfJointSourceObservations;
	/**
	 * Storage for destination observations for addObservsations
	 */
	private Vector<double[][]> vectorOfJointDestinationObservations;

	public void initialise(int k) throws Exception {
		// Assume the user only wants 1 source and 1 destination dimension
		initialise(k, 1, 1);
	}

	public void initialise(int k, int sourceDimensions, int destDimensions) throws Exception {
		super.initialise(k);
		initialise(sourceDimensions, destDimensions);
	}
	
	/**
	 * Initialise using existing or default value of k
	 */
	public void initialise(int sourceDimensions, int destDimensions) throws Exception {
		this.destDimensions = destDimensions;
		this.sourceDimensions = sourceDimensions;
	}

	@Override
	protected void initialiseKraskovCalculators() throws Exception {
		mickPastToSource.initialise(k * destDimensions, sourceDimensions);
		mickNextPastToSource.initialise((k+1) * destDimensions, sourceDimensions);
	}

	/**
	 * Set the observations to compute the probabilities from 
	 * 
	 * @param source
	 * @param destination
	 */
	public void setObservations(double[][] source, double[][] destination) throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}
	
	@Override
	public void startAddObservations() {
		vectorOfJointSourceObservations = new Vector<double[][]>();
		vectorOfJointDestinationObservations = new Vector<double[][]>();
	}

	/**
	 * Add observations of a single-dimensional source and destination pair.
	 * 
	 * Only allow this call if source and destination dimenions were 1.
	 * 
	 * @param source
	 * @param destination
	 */
	@Override
	public void addObservations(double[] source, double[] destination) throws Exception {
		double[][] sourceMatrix = new double[source.length][1];
		MatrixUtils.copyIntoColumn(sourceMatrix, 0, source);
		double[][] destMatrix = new double[destination.length][1];
		MatrixUtils.copyIntoColumn(destMatrix, 0, destination);
		addObservations(sourceMatrix, destMatrix);
	}

	/**
	 * Add observations of a single-dimensional source and destination pair.
	 * 
	 * @param source
	 * @param destination
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 */
	@Override
	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception {
		double[][] sourceMatrix = new double[numTimeSteps][1];
		MatrixUtils.copyIntoColumn(sourceMatrix, 0, 0, source, startTime, numTimeSteps);
		double[][] destMatrix = new double[destination.length][1];
		MatrixUtils.copyIntoColumn(destMatrix, 0, 0, destination, startTime, numTimeSteps);
		addObservations(sourceMatrix, destMatrix);
	}

	/**
	 * Add observations of the joint source and destinations
	 * 
	 * @param source
	 * @param destination
	 * @throws Exception
	 */
	public void addObservations(double[][] source, double[][] destination) throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		int thisSourceDimensions = source[0].length;
		int thisDestDimensions = destination[0].length;
		if ((thisDestDimensions != destDimensions) || (thisSourceDimensions != sourceDimensions)) {
			throw new Exception("Cannot add observsations for source and destination variables " +
					" of " + thisSourceDimensions + " and " + thisDestDimensions +
					" dimensions respectively for TE calculator set up for " + sourceDimensions + " " +
					destDimensions + " source and destination dimensions respectively");
		}
		if (vectorOfJointSourceObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		vectorOfJointSourceObservations.add(source);
		vectorOfJointDestinationObservations.add(destination);
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
		double[][] sourceToAdd = new double[numTimeSteps][source[0].length];
		System.arraycopy(source, startTime, sourceToAdd, 0, numTimeSteps);
		double[][] destToAdd = new double[numTimeSteps][destination[0].length];
		System.arraycopy(destination, startTime, destToAdd, 0, numTimeSteps);
		addObservations(sourceToAdd, destToAdd);
	}

	/**
	 * Flag that the observations are complete, probability distribution functions can now be built.
	 *
	 */
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[][] destination : vectorOfJointDestinationObservations) {
			totalObservations += destination.length - k;
		}
		jointPastVectors = new double[totalObservations][k * destDimensions];
		jointNextAndPastVectors = new double[totalObservations][(k+1) * destDimensions];
		sourceVectors = new double[totalObservations][sourceDimensions];

		// Construct the joint vectors from the given observations
		int startObservation = 0;
		Iterator<double[][]> iterator = vectorOfJointDestinationObservations.iterator();
		for (double[][] source : vectorOfJointSourceObservations) {
			double[][] destination = iterator.next();
			double[][] currentDestPastVectors = makeJointVectorForPast(destination);
			MatrixUtils.arrayCopy(currentDestPastVectors, 0, 0,
					jointPastVectors, startObservation, 0, currentDestPastVectors.length,
					k * destDimensions);
			double[][] currentDestNextPastVectors = makeJointVectorForNextPast(destination);
			MatrixUtils.arrayCopy(currentDestNextPastVectors, 0, 0,
					jointNextAndPastVectors, startObservation, 0,
					currentDestNextPastVectors.length, (k + 1) * destDimensions);
			MatrixUtils.arrayCopy(source, k-1, 0, sourceVectors, startObservation, 0,
					source.length - k, sourceDimensions);
			startObservation += destination.length - k;
		}
		
		// Now set the joint vectors in the kernel estimators
		try {
			mickPastToSource.setObservations(jointPastVectors, sourceVectors);
			mickNextPastToSource.setObservations(jointNextAndPastVectors, sourceVectors);
		} catch (Exception e) {
			// The above should not throw an exception since they were constructed here
			//  of the same time length, so wrap in a runtime exception
			throw new RuntimeException(e);
		}

		// Store whether there was more than one observation set:
		addedMoreThanOneObservationSet = vectorOfJointDestinationObservations.size() > 1;
		
		// And clear the vector of observations
		vectorOfJointSourceObservations = null;
		vectorOfJointDestinationObservations = null;
	}

	/**
	 * If mickPastToSource can utilise anything from mickPastNextToSource after the latter
	 * has run computeAverageLocalOfObservations, arrange that here 
	 *
	 */
	protected void shareDataBetweenUnderlyingCalculators() {
		if (! MutualInfoCalculatorMultiVariateKraskovByMulti.class.isInstance(mickNextPastToSource)) {
			// We don't know of what to share for other calculator types.
			// Subclasses may know and can over-ride this method.
			return;
		}
		MutualInfoCalculatorMultiVariateKraskovByMulti micmvkNextPastToSource =
			(MutualInfoCalculatorMultiVariateKraskovByMulti) mickNextPastToSource;
		MutualInfoCalculatorMultiVariateKraskovByMulti micmvkPastToSource =
			(MutualInfoCalculatorMultiVariateKraskovByMulti) mickPastToSource;
		if (micmvkNextPastToSource.multiInfoJoint.norms != null) {
			// Share the norms already computed with the other mutual info calculator.
			// Just assign the joint norms for it, it will filter them through to the
			//  marginals itself.
			if (micmvkPastToSource.multiInfoJoint.norms == null) {
				micmvkPastToSource.multiInfoJoint.norms =
					new double[destDimensions * k + sourceDimensions][][];
				// Point to the norms for the past variables
				for (int t = 0; t < k; t++) {
					for (int d = 0; d < destDimensions; d++) {
						micmvkPastToSource.multiInfoJoint.norms[t*destDimensions + d] =
							micmvkNextPastToSource.multiInfoJoint.norms[(1+t)*destDimensions + d];
					}
				}
				// Point to the norms for the source variable
				for (int s = 0; s < sourceDimensions; s++) {
					micmvkPastToSource.multiInfoJoint.norms[k*destDimensions + s] =
						micmvkNextPastToSource.multiInfoJoint.norms[(k + 1)*destDimensions + s];
				}
			}
		}
	}

	/**
	 * Generate a vector for each time step, containing the past k states of the destination.
	 * Note that each state of the destination is a joint vector of destDimensions variables.
	 * Does not include a vector for the first k time steps.
	 * 
	 * @param destination
	 * @return array of vectors for each time step
	 */
	private double[][] makeJointVectorForPast(double[][] destination) {
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
	 * Sets the observations to compute the PDFs from.
	 * Cannot be called in conjunction with start/add/finaliseAddObservations.
	 * destValid is a time series (with time indices the same as destination)
	 *  indicating whether the destination at that point is valid.
	 * sourceValid is the same for the source
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param sourceValid
	 * @param destValid
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

	public void setObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		
		boolean[] jointSourceValid = MatrixUtils.andRows(sourceValid);
		boolean[] jointDestValid = MatrixUtils.andRows(destValid);
		setObservations(source, destination, jointSourceValid, jointDestValid);
	}

	/**
	 * Generate a vector for each time step, containing the past k states of
	 *  the destination, and the current state.
	 * Does not include a vector for the first k time steps.
	 * 
	 * @param destination
	 * @return
	 */
	protected double[][] makeJointVectorForNextPast(double[][] destination) {
		// We want all delay vectors here
		return MatrixUtils.makeDelayEmbeddingVector(destination, k+1);
	}

	public void setDebug(boolean debug) {
		super.setDebug(debug);
		mickNextPastToSource.setDebug(debug);
		mickPastToSource.setDebug(debug);
	}
}
