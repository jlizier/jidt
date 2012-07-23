package infodynamics.measures.continuous.kernel;

import infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariate;
import infodynamics.measures.continuous.TransferEntropyCommon;
import infodynamics.measures.continuous.kernel.TransferEntropyKernelCounts;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.MeasurementDistribution;
import infodynamics.utils.RandomGenerator;

import java.util.Iterator;
import java.util.Vector;


/**
 * 
 * <p>
 * Implements a transfer entropy calculator using kernel estimation.
 * (see Schreiber, PRL 85 (2) pp.461-464, 2000)</p> 
 * 
 * <p>
 * This calculator handles multi-variate source and destination variables.
 * </p>
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct</li>
 * 		<li>SetProperty() for each property</li>
 *		<li>initialise()</li>
 * 		<li>setObservations(), or [startAddObservations(), addObservations()*, finaliseAddObservations()]
 *   Note: If not using setObservations(), the results from computeLocal or getSignificance
 *    are not likely to be particularly sensible.</li> 
 * 		<li>computeAverageLocalOfObservations() or ComputeLocalOfPreviousObservations()</li>
 * 	</ol>
 * </p>
 * 
 * <p>
 * TODO Implement dynamic correlation exclusion with multiple observation sets. (see the
 *  way this is done in Plain calculator).
 * </p>
 * 
 * @author Joseph Lizier
 * @see For transfer entropy: Schreiber, PRL 85 (2) pp.461-464, 2000; http://dx.doi.org/10.1103/PhysRevLett.85.461
 * @see For local transfer entropy: Lizier et al, PRE 77, 026110, 2008; http://dx.doi.org/10.1103/PhysRevE.77.026110
 *
 */
public class TransferEntropyCalculatorMultiVariateKernel
	extends TransferEntropyCommon implements TransferEntropyCalculatorMultiVariate {

	protected KernelEstimatorTransferEntropyMultiVariate teKernelEstimator = null;
	// Keep a kernel estimator for the next state, in case we wish to compute
	//  Active info storage also:
	protected KernelEstimatorMultiVariate nextStateKernelEstimator = null;

	/**
	 * Storage for source observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfJointSourceObservations;
	/**
	 * Storage for destination observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfJointDestinationObservations;

	// Keep joint vectors so we don't need to regenerate them
	protected double[][] destPastVectors;
	protected double[][] destNextVectors;
	protected double[][] sourceVectors;
	
	protected int destDimensions = 1;
	protected int sourceDimensions = 1;
	
	// Store the local conditional probability of next on past state as 
	//  computed during a local TE computation, in case the caller wants to
	//  compute the local active info storage next.
	protected double[] localProbNextCondPast;
	
	protected boolean normalise = true;
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	protected boolean dynCorrExcl = false;
	protected int dynCorrExclTime = 100;
	public static final String DYN_CORR_EXCL_TIME_NAME = "DYN_CORR_EXCL";
	
	protected boolean forceCompareToAll = false;
	public static final String FORCE_KERNEL_COMPARE_TO_ALL = "FORCE_KERNEL_COMPARE_TO_ALL";
	
	/**
	 * Default value for epsilon
	 */
	public static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	protected double epsilon = DEFAULT_EPSILON;
	public static final String EPSILON_PROP_NAME = "EPSILON";
	
	/**
	 * Creates a new instance of the kernel-estimate style transfer entropy calculator
	 *
	 */
	public TransferEntropyCalculatorMultiVariateKernel() {
		super();
		teKernelEstimator = new KernelEstimatorTransferEntropyMultiVariate();
		teKernelEstimator.setNormalise(normalise);
		nextStateKernelEstimator = new KernelEstimatorMultiVariate();
		nextStateKernelEstimator.setNormalise(normalise);
	}

	/**
	 * Initialises the calculator with the existing value for epsilon
	 * 
	 * @param k history length
	 */
	public void initialise(int k) throws Exception {
		initialise(k, epsilon);
	}
	
	/**
	 * Initialises the calculator
	 * 
	 * @param k history length
	 * @param epsilon kernel width
	 */
	public void initialise(int k, double epsilon) throws Exception {
		this.epsilon = epsilon;
		initialise(k, 1, 1); // assume 1 dimension in source and dest
	}

	public void initialise(int k, int sourceDimensions, int destDimensions) throws Exception {
		this.destDimensions = destDimensions;
		this.sourceDimensions = sourceDimensions;
		super.initialise(k); // calls initialise();
	}
	
	public void initialise(int sourceDimensions, int destDimensions) throws Exception {
		this.destDimensions = destDimensions;
		this.sourceDimensions = sourceDimensions;
		super.initialise(k); // calls initialise();
	}

	/**
	 * Initialise using default or existing values for k and epsilon
	 */
	public void initialise() {
		teKernelEstimator.initialise(k * destDimensions,
				sourceDimensions, epsilon, epsilon);
		nextStateKernelEstimator.initialise(destDimensions, epsilon);
		destPastVectors = null;
		destNextVectors = null;
		sourceVectors = null;
		localProbNextCondPast = null;
	}
	
	/**
	 * Set properties for the transfer entropy calculator.
	 * These can include:
	 * <ul>
	 * 		<li>K_PROP_NAME</li>
	 * 		<li>EPSILON_PROP_NAME</li>
	 * 		<li>NORMALISE_PROP_NAME</li>
	 * 		<li>DYN_CORR_EXCL_TIME_NAME</li>
	 * 		<li>FORCE_KERNEL_COMPARE_TO_ALL</li>
	 * </ul> 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception 
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		super.setProperty(propertyName, propertyValue);
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			epsilon = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
			teKernelEstimator.setNormalise(normalise);
			nextStateKernelEstimator.setNormalise(normalise);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
			if (dynCorrExcl) {
				teKernelEstimator.setDynamicCorrelationExclusion(dynCorrExclTime);
				nextStateKernelEstimator.setDynamicCorrelationExclusion(dynCorrExclTime);
			} else {
				teKernelEstimator.clearDynamicCorrelationExclusion();
				nextStateKernelEstimator.clearDynamicCorrelationExclusion();
			}
		} else if (propertyName.equalsIgnoreCase(FORCE_KERNEL_COMPARE_TO_ALL)) {
			forceCompareToAll = Boolean.parseBoolean(propertyValue);
			teKernelEstimator.setForceCompareToAll(forceCompareToAll);
			nextStateKernelEstimator.setForceCompareToAll(forceCompareToAll);
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
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
	public void setObservations(double[][] source, double[][] destination, boolean[] sourceValid, boolean[] destValid) throws Exception {
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
		destPastVectors = new double[totalObservations][k * destDimensions];
		destNextVectors = new double[totalObservations][destDimensions];
		sourceVectors = new double[totalObservations][sourceDimensions];
		
		// Construct the joint vectors from the given observations
		int startObservation = 0;
		Iterator<double[][]> iterator = vectorOfJointDestinationObservations.iterator();
		for (double[][] source : vectorOfJointSourceObservations) {
			double[][] destination = iterator.next();
			double[][] currentDestPastVectors = makeJointVectorForPast(destination);
			MatrixUtils.arrayCopy(currentDestPastVectors, 0, 0,
					destPastVectors, startObservation, 0, currentDestPastVectors.length,
					k * destDimensions);
			MatrixUtils.arrayCopy(destination, k, 0,
					destNextVectors, startObservation, 0,
					destination.length - k, destDimensions);
			MatrixUtils.arrayCopy(source, k - 1, 0,
					sourceVectors, startObservation, 0,
					source.length - k, sourceDimensions);
			startObservation += destination.length - k;
		}
		
		// Now set the joint vectors in the kernel estimators
		teKernelEstimator.setObservations(destPastVectors, destNextVectors, sourceVectors);

		// Store whether there was more than one observation set:
		addedMoreThanOneObservationSet = vectorOfJointDestinationObservations.size() > 1;
		
		// And clear the vector of observations
		vectorOfJointSourceObservations = null;
		vectorOfJointDestinationObservations = null;
	}
	
	/**
	 * <p>Computes the average Transfer Entropy for the previously supplied observations</p> 
	 * 
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		
		double te = 0.0;
		if (debug) {
			MatrixUtils.printMatrix(System.out, destPastVectors);
		}
		for (int b = 0; b < totalObservations; b++) {
			TransferEntropyKernelCounts kernelCounts = teKernelEstimator.getCount(destPastVectors[b],
					destNextVectors[b], sourceVectors[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (kernelCounts.countNextPastSource > 0) {
				logTerm = ((double) kernelCounts.countNextPastSource / (double) kernelCounts.countPastSource) /
							((double) kernelCounts.countNextPast / (double) kernelCounts.countPast);
				cont = Math.log(logTerm);
			}
			te += cont;
			if (debug) {
				System.out.println(b + ": " + destPastVectors[b][0] + " (" +
						kernelCounts.countNextPastSource + " / " + kernelCounts.countPastSource + ") / (" +
						kernelCounts.countNextPast + " / " + kernelCounts.countPast + ") = " + 
						logTerm + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
		}
		lastAverage = te / (double) totalObservations / Math.log(2.0);
		return lastAverage;
	}

	/**
	 * <p>Computes the average Transfer Entropy for the previously supplied observations,
	 *  using the Grassberger correction for the point count k: log_e(k) ~= digamma(k).</p>
	 * <p>Kaiser and Schreiber, Physica D 166 (2002) pp. 43-62 suggest (on p. 57) that for the TE
	 *   though the adverse correction of the bias correction is worse than the correction
	 *   itself (because the probabilities being multiplied/divided are not independent), 
	 *   so recommend not to use this method.
	 * </p>
	 * <p>It is implemented here for testing purposes only.</p>
	 * 
	 */
	public double computeAverageLocalOfObservationsWithCorrection() throws Exception {
		
		double te = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			TransferEntropyKernelCounts kernelCounts = teKernelEstimator.getCount(destPastVectors[b],
					destNextVectors[b], sourceVectors[b], b);
			double cont = 0.0;
			if (kernelCounts.countNextPastSource > 0) {
				cont = MathsUtils.digamma(kernelCounts.countNextPastSource) - 
						MathsUtils.digamma(kernelCounts.countPastSource) -
						MathsUtils.digamma(kernelCounts.countNextPast) +
						MathsUtils.digamma(kernelCounts.countPast);
			}
			te += cont;
			/*
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
			*/
		}
		// Average it, and convert results to bytes
		lastAverage = te / (double) totalObservations / Math.log(2.0);
		return lastAverage;
	}

	/**
	 * Computes the local transfer entropies for the previous supplied observations.
	 * 
	 * Where more than one time series has been added, the array
	 *  contains the local values for each tuple in the order in
	 *  which they were added.
	 * 
	 * If there was only a single time series added, the array
	 *  contains k zero values before the local values.
	 *  (This means the length of the return array is the same
	 *  as the length of the input time series).
	 * 
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		return computeLocalUsingPreviousObservations(null, null, true);
	}

	/**
	 * Comptues local transfer entropies for the given observations, using the previously supplied
	 * observations to compute the PDFs.
	 * I don't think it's such a good idea to do this for continuous variables (e.g. where
	 * one can get kernel estimates for probabilities of zero now) but I've implemented
	 * it anyway. I guess getting kernel estimates of zero here is no different than what
	 * can occur with dynamic correlation exclusion.
	 * 
	 * @param source
	 * @param destination
	 * @return
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] source, double[][] destination) throws Exception {
		return computeLocalUsingPreviousObservations(source, destination, false);
	}
	
	/**
	 * Returns the local TE at every time point.
	 * 
	 * @param source
	 * @param destination
	 * @param isPreviousObservations
	 * @return
	 * @throws Exception
	 */
	private double[] computeLocalUsingPreviousObservations(double[][] source,
			double[][] destination, boolean isPreviousObservations) throws Exception {
		
		double[][] newDestPastVectors;
		double[][] newDestNextValues;
		double[][] newSourceValues;

		if (isPreviousObservations) {
			// We've already computed the joint vectors for these observations
			newDestPastVectors = destPastVectors;
			newDestNextValues = destNextVectors;
			newSourceValues = sourceVectors;
		} else {
			// We need to compute a new set of joint vectors
			newDestPastVectors = makeJointVectorForPast(destination);
			newDestNextValues = new double[destination.length - k][destDimensions];
			MatrixUtils.arrayCopy(destination, k, 0,
					newDestNextValues, 0, 0,
					destination.length - k, destDimensions);
			newSourceValues = new double[source.length - k][sourceDimensions];
			MatrixUtils.arrayCopy(source, k - 1, 0,
					newSourceValues, 0, 0,
					source.length - k, sourceDimensions);
		}

		double te = 0.0;
		int numLocalObservations = newDestPastVectors.length;
		double[] localTE;
		int offset = 0;
		if (isPreviousObservations && addedMoreThanOneObservationSet) {
			// We're returning the local values for a set of disjoint
			//  observations. So we don't add k zeros to the start
			localTE = new double[numLocalObservations];
			offset = 0;
		} else {
			localTE = new double[numLocalObservations + k];
			offset = k;
		}
		localProbNextCondPast = new double[numLocalObservations];
		double avKernelCount = 0;
		TransferEntropyKernelCounts kernelCounts;
		for (int b = 0; b < numLocalObservations; b++) {
			// System.out.print("Observation number " + String.valueOf(b) + "\n");
			if (isPreviousObservations) {
				kernelCounts = teKernelEstimator.getCount(
						newDestPastVectors[b],
						newDestNextValues[b], newSourceValues[b], b);
			} else {
				kernelCounts = teKernelEstimator.getCount(
						newDestPastVectors[b],
						newDestNextValues[b], newSourceValues[b], -1);
			}
			avKernelCount += kernelCounts.countNextPastSource;
			double logTerm = 0.0;
			double local = 0.0;
			if (kernelCounts.countPast > 0) {
				// Store this ratio for a potential active info calculation later
				localProbNextCondPast[b] = (double) kernelCounts.countNextPast / (double) kernelCounts.countPast;
			}
			if (kernelCounts.countNextPastSource > 0) {
				logTerm = ((double) kernelCounts.countNextPastSource / (double) kernelCounts.countPastSource) /
							localProbNextCondPast[b];
				local = Math.log(logTerm);
			}
			localTE[offset + b] = local;
			te += local;
			/*
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (local/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
			*/
		}
		avKernelCount = avKernelCount / (double) numLocalObservations;
		if (debug) {
			System.out.printf("Average kernel count was %.3f\n", avKernelCount);
		}
		lastAverage = te / (double) numLocalObservations / Math.log(2.0);
		return localTE;
	}

	/**
	 * Compute the significance of obtaining the given average TE from the given observations
	 * 
	 * 	This is as per Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128.
	 *
	 * Basically, we shuffle the source observations against the destination tuples.
	 * This keeps the marginal PDFs the same (including the entropy rate of the destination)
	 *  but destroys any correlation between the source and state change of the destination.
	 * 
	 * @param numPermutationsToCheck number of new orderings of the source values to compare against
	 * @return
	 */
	public MeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(totalObservations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * As per {@link computeSignificance(int) computeSignificance()} but supplies
	 *  the re-orderings of the observations of the source variables.
	 *  
	 * 
	 * @param newOrderings first index is permutation number, i.e. newOrderings[i]
	 * 		is an array of 1 permutation of 0..n-1, where there were n observations.
	 * @return
	 * @throws Exception
	 */
	public MeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		
		double actualTE = computeAverageLocalOfObservations();
		
		// Space for the source observations:
		double[][] oldSourceValues = sourceVectors;
		
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		MeasurementDistribution measDistribution = new MeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the source in the destPastSourceVectors 
			//  and destNextPastSourceVectors vectors
			sourceVectors = MatrixUtils.extractSelectedTimePoints(oldSourceValues, newOrderings[p]);
			
			// Make the equivalent operations of intialise
			teKernelEstimator.initialise(k * destDimensions,
					sourceDimensions, epsilon, epsilon);
			// Make the equivalent operations of setObservations:
			teKernelEstimator.setObservations(destPastVectors, destNextVectors, sourceVectors);
			// And get a TE value for this realisation:
			double newTe = computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newTe;
			if (newTe >= actualTE) {
				countWhereTeIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the local variables:
		lastAverage = actualTE;
		sourceVectors = oldSourceValues;
		// And set the kernel estimator back to their previous state
		teKernelEstimator.initialise(k * destDimensions,
				sourceDimensions, epsilon, epsilon);
		teKernelEstimator.setObservations(destPastVectors, destNextVectors, sourceVectors);
		
		// And return the significance
		measDistribution.pValue = (double) countWhereTeIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualTE;
		return measDistribution;
	}

	/**
	 * Computes the local active info storage for the previous supplied observations.
	 * 
	 * Where more than one time series has been added, the array
	 *  contains the local values for each tuple in the order in
	 *  which they were added.
	 * 
	 * If there was only a single time series added, the array
	 *  contains k zero values before the local values.
	 *  (This means the length of the return array is the same
	 *  as the length of the input time series).
	 *  
	 * Precondition: The user must have computed the local TEs first for these
	 *  vectors.
	 * 
	 */
	public double[] computeLocalActiveOfPreviousObservations() throws Exception {
		return computeLocalActiveUsingPreviousObservations(null, true);
	}

	/**
	 * Comptues local active info storage for the given observations, using the previously supplied
	 * observations to compute the PDFs.
	 * I don't think it's such a good idea to do this for continuous variables (e.g. where
	 * one can get kernel estimates for probabilities of zero now) but I've implemented
	 * it anyway. I guess getting kernel estimates of zero here is no different than what
	 * can occur with dynamic correlation exclusion.
	 * 
	 * Precondition: The user must have computed the local TEs first for these
	 *  vectors
	 * 
	 * @param source
	 * @param destination
	 * @return
	 * @throws Exception
	 */
	public double[] computeLocalActiveUsingPreviousObservations(double[][] destination) throws Exception {
		return computeLocalActiveUsingPreviousObservations(destination, false);
	}
	
	/**
	 * Returns the local active info at every time point.
	 * 
	 * Precondition: The user must have computed the local TEs first for these
	 *   vectors.
	 * 
	 * @param source
	 * @param destination
	 * @param isPreviousObservations
	 * @return
	 * @throws Exception
	 */
	private double[] computeLocalActiveUsingPreviousObservations(
			double[][] destination, boolean isPreviousObservations) throws Exception {

		// Precondition: the local TE must have already been computed
		if (localProbNextCondPast == null) {
			throw new RuntimeException("A local TE must have been computed before " +
					"the local active info storage can be computed by TransferEntropyCalculatorMultiVariateKernel");
		}
		
		double[][] newDestNextValues;

		// Set the observations on the kernel estimator:
		nextStateKernelEstimator.setObservations(destNextVectors);
		
		// Now set which observations we're going to compute the local 
		//  active info of:
		if (isPreviousObservations) {
			// We've already computed the joint vectors for these observations
			newDestNextValues = destNextVectors;
		} else {
			// We need to compute a new set of joint vectors
			newDestNextValues = new double[destination.length - k][destDimensions];
			MatrixUtils.arrayCopy(destination, k, 0,
					newDestNextValues, 0, 0,
					destination.length - k, destDimensions);
		}

		int numLocalObservations = newDestNextValues.length;
		double[] localActive;
		int offset = 0;
		if (isPreviousObservations && addedMoreThanOneObservationSet) {
			// We're returning the local values for a set of disjoint
			//  observations. So we don't add k zeros to the start
			localActive = new double[numLocalObservations];
			offset = 0;
		} else {
			localActive = new double[numLocalObservations + k];
			offset = k;
		}
		double nextStateProb;
		for (int b = 0; b < numLocalObservations; b++) {
			if (isPreviousObservations) {
				nextStateProb = nextStateKernelEstimator.getProbability(
						newDestNextValues[b], b);
			} else {
				nextStateProb = nextStateKernelEstimator.getProbability(
						newDestNextValues[b], -1);
			}
			double logTerm = 0.0;
			double local = 0.0;
			if (localProbNextCondPast[b] > 0) {
				logTerm = localProbNextCondPast[b] / nextStateProb;
				local = Math.log(logTerm);
			}
			localActive[offset + b] = local;
			/*
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (local/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
			*/
		}
		return localActive;
	}
	
	public void setDebug(boolean debug) {
		super.setDebug(debug);
		teKernelEstimator.setDebug(debug);
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
}
