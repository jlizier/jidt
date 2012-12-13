package infodynamics.measures.continuous.kernel;

import infodynamics.utils.MatrixUtils;

import java.util.Iterator;
import java.util.Vector;

/**
 * 
 * <p>
 * Implements a transfer entropy calculator using kernel estimation.
 * (see Schreiber, PRL 85 (2) pp.461-464, 2000)</p> 
 * 
 * <p>
 * This calculator handles multi-variate source and destination variables,
 *  and should only be used to add observation tuples, i.e.
 *  (source, destination next state, destination past)
 *  one at a time. This allows the user to specify the variable that
 *  should be used as the destination past, for advanced applications.
 * </p>
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct</li>
 * 		<li>SetProperty() for each property</li>
 *		<li>initialise()</li>
 * 		<li>setObservations(), or [startAddObservations(),
 * 			(addObservations|addSingleObservation)*, finaliseAddObservations()]
 *   Note: If not using setObservations(), the results from computeLocal or getSignificance
 *    are not guaranteed to be particularly sensible.</li> 
 * 		<li>computeAverageLocalOfObservations() or ComputeLocalOfPreviousObservations()</li>
 * 	</ol>
 * </p>
 * 
 * <p>
 * TODO Implement dynamic correlation exclusion with multiple observation sets. (see the
 *  way this is done in Plain calculator).
 * TODO Think about added error-trapping code to make sure the user only makes one type of addObservations call.
 * </p>
 * 
 * @author Joseph Lizier
 * @see For transfer entropy: Schreiber, PRL 85 (2) pp.461-464, 2000; http://dx.doi.org/10.1103/PhysRevLett.85.461
 * @see For local transfer entropy: Lizier et al, PRE 77, 026110, 2008; http://dx.doi.org/10.1103/PhysRevE.77.026110
 *
 */
public class TransferEntropyCalculatorMultiVariateSingleObservationsKernel
		extends TransferEntropyCalculatorMultiVariateKernel {

	/**
	 * Storage for destination history observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfJointDestinationPastObservations;

	protected int destPastDimensions = 1;

	public TransferEntropyCalculatorMultiVariateSingleObservationsKernel() {
		super();
	}

	/**
	 * Initialises the calculator
	 * 
	 * @param epsilon kernel width
	 */
	public void initialise(double epsilon) throws Exception {
		this.epsilon = epsilon;
		initialise(1, 1); // assume 1 dimension in source and dest
	}

	public void initialise(int sourceDimensions, int destDimensions) throws Exception {
		this.destDimensions = destDimensions;
		this.sourceDimensions = sourceDimensions;
		this.destPastDimensions = destDimensions; // assume same
		super.initialise(1); // Feeds k=1 to super and calls initialise();
	}

	public void initialiseAllDimensions(int sourceDimensions,
			int destDimensions, int destPastDimensions) throws Exception {
		this.destDimensions = destDimensions;
		this.sourceDimensions = sourceDimensions;
		this.destPastDimensions = destPastDimensions;
		
		// Mimic super.initialise(1) (it would replace dest and source
		//  dimensions if we're not careful)
		addedMoreThanOneObservationSet = false;
		k = 1;
		// Mimic super.initialise() (it would use k * destDimenions in the kernel estimator
		//  for destPast instead of destPastDimensions if we're not careful)
		teKernelEstimator.initialise(destPastDimensions,
				sourceDimensions, epsilon, epsilon);
		nextStateKernelEstimator.initialise(destDimensions, epsilon);
		destPastVectors = null;
		destNextVectors = null;
		sourceVectors = null;
		localProbNextCondPast = null;
	}

	/**
	 * Set the observations to compute the probabilities from 
	 * 
	 * @param source
	 * @param destination
	 */
	public void setObservations(double[][] source, double[][] destination,
			double[][] destinationPast) throws Exception {
		startAddObservations();
		addObservations(source, destination, destinationPast);
		finaliseAddObservations();
	}
	
	@Override
	public void startAddObservations() {
		vectorOfJointDestinationPastObservations = new Vector<double[][]>();
		super.startAddObservations();
	}

	/**
	 * Add observations of the joint source and destinations
	 * 
	 * @param source
	 * @param destination
	 * @param destinationPast
	 * @throws Exception
	 */
	public void addObservations(double[][] source, double[][] destination,
			double[][] destinationPast) throws Exception {
		if (destinationPast.length != destination.length) {
			throw new Exception(String.format("Destination past and destination lengths (%d and %d) must match!",
					destinationPast.length, destination.length));
		}
		int thisDestPastDimensions = destinationPast[0].length;
		if ((thisDestPastDimensions != destPastDimensions)) {
			throw new Exception("Cannot add observsations for destination past variables " +
					" of " + thisDestPastDimensions +
					" dimensions for TE calculator set up for " + destPastDimensions +
					" destination past dimensions");
		}
		if (vectorOfJointDestinationPastObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		vectorOfJointDestinationPastObservations.add(destinationPast);
		super.addObservations(source, destination);
	}

	/**
	 * Add a single observation of the joint source, destinations and
	 * destination past
	 * 
	 * @param source
	 * @param destination
	 * @param destinationPast
	 * @throws Exception
	 */
	public void addSingleObservation(double[] source, double[] destination,
			double[] destinationPast) throws Exception {
		int thisSourceDimensions = source.length;
		int thisDestDimensions = destination.length;
		int thisDestPastDimensions = destinationPast.length;
		if ((thisDestDimensions != destDimensions) ||
				(thisSourceDimensions != sourceDimensions) ||
				(thisDestPastDimensions != destPastDimensions)) {
			throw new Exception("Cannot add observsations for source, destination and destPast variables " +
					" of " + thisSourceDimensions + " and " + thisDestDimensions + " and " +
					thisDestPastDimensions +
					" dimensions respectively for TE calculator set up for " + sourceDimensions + ", " +
					destDimensions + " and " + destPastDimensions +
					" source, destination and destPast dimensions respectively");
		}
		if (vectorOfJointDestinationPastObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		// Now make the multidimensional arrays to add in
		double[][] sourceContainer = new double[1][];
		double[][] destContainer = new double[1][];
		double[][] destPastContainer = new double[1][];
		sourceContainer[0] = source;
		destContainer[0] = destination;
		destPastContainer[0] = destinationPast;
		vectorOfJointDestinationPastObservations.add(destPastContainer);
		super.addObservations(sourceContainer, destContainer);
	}

	/**
	 * Flag that the observations are complete, probability distribution functions can now be built.
	 *
	 */
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[][] destination : vectorOfJointDestinationObservations) {
			// No need t jump k values in, since we've got the destination
			//  past values held separately
			totalObservations += destination.length;
		}
		destPastVectors = new double[totalObservations][destPastDimensions];
		destNextVectors = new double[totalObservations][destDimensions];
		sourceVectors = new double[totalObservations][sourceDimensions];
		
		// Construct the joint vectors from the given observations
		int startObservation = 0;
		Iterator<double[][]> iterator = vectorOfJointDestinationObservations.iterator();
		Iterator<double[][]> iteratorDestPast = vectorOfJointDestinationPastObservations.iterator();
		for (double[][] source : vectorOfJointSourceObservations) {
			double[][] destination = iterator.next();
			double[][] destinationPast = iteratorDestPast.next();
			// Add in all observations - no need to offset by k since
			//  we've got the destination past held separately.
			MatrixUtils.arrayCopy(destinationPast, 0, 0,
					destPastVectors, startObservation, 0, destinationPast.length,
					destPastDimensions);
			MatrixUtils.arrayCopy(destination, 0, 0,
					destNextVectors, startObservation, 0,
					destination.length, destDimensions);
			MatrixUtils.arrayCopy(source, 0, 0,
					sourceVectors, startObservation, 0,
					source.length, sourceDimensions);
			startObservation += destination.length;
		}
		
		// Now set the joint vectors in the kernel estimators
		teKernelEstimator.setObservations(destPastVectors, destNextVectors, sourceVectors);

		// Store whether there was more than one observation set:
		addedMoreThanOneObservationSet = vectorOfJointDestinationObservations.size() > 1;
		
		// And clear the vector of observations
		vectorOfJointSourceObservations = null;
		vectorOfJointDestinationObservations = null;
		vectorOfJointDestinationPastObservations = null;
	}

}
