package infodynamics.measures.continuous.kernel;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

public class MutualInfoCalculatorMultiVariateKernel implements
	MutualInfoCalculatorMultiVariate {

	KernelEstimatorMultiVariate mvke1 = null;
	KernelEstimatorMultiVariate mvke2 = null;
	KernelEstimatorMultiVariate mvkeJoint = null;

	private int totalObservations = 0;
	// private int dimensions1 = 0;
	// private int dimensions2 = 0;
	private boolean debug = false;
	private double[][] observations1; 
	private double[][] observations2; 
	private double lastAverage;
	private boolean miComputed;
	
	private boolean normalise = true;
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	private boolean dynCorrExcl = false;
	private int dynCorrExclTime = 100;
	public static final String DYN_CORR_EXCL_TIME_NAME = "DYN_CORR_EXCL";
	
	private boolean forceCompareToAll = false;
	public static final String FORCE_KERNEL_COMPARE_TO_ALL = "FORCE_KERNEL_COMPARE_TO_ALL";
	
	/**
	 * Default value for epsilon
	 */
	public static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	private double epsilon = DEFAULT_EPSILON;
	public static final String EPSILON_PROP_NAME = "EPSILON";
	
	public MutualInfoCalculatorMultiVariateKernel() {
		mvke1 = new KernelEstimatorMultiVariate();
		mvke2 = new KernelEstimatorMultiVariate();
		mvkeJoint = new KernelEstimatorMultiVariate();
		mvke1.setNormalise(normalise);
		mvke2.setNormalise(normalise);
		mvkeJoint.setNormalise(normalise);
	}

	/**
	 * Initialise using a default epsilon
	 * 
	 * @param dimensions1
	 * @param dimensions2
	 */
	public void initialise(int dimensions1, int dimensions2) {
		initialise(dimensions1, dimensions2, epsilon);
	}

	public void initialise(int dimensions1, int dimensions2, double epsilon) {
		this.epsilon = epsilon;
		mvke1.initialise(dimensions1, epsilon);
		mvke2.initialise(dimensions2, epsilon);
		mvkeJoint.initialise(dimensions1 + dimensions2, epsilon);
		// this.dimensions1 = dimensions1;
		// this.dimensions2 = dimensions2;
		lastAverage = 0.0;
		miComputed = false;
	}

	public void addObservations(double[][] source, double[][] destination) throws Exception {
		// TODO If we ever implement these (which will require changing the kernel
		//  estimators) we will need to throw an exception if dynamic correlation
		//  exclusion was set.
		throw new RuntimeException("Not implemented yet");
	}

	public void addObservations(double[][] source, double[][] destination, int startTime, int numTimeSteps) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination, boolean[] sourceValid, boolean[] destValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination, boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void startAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void finaliseAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	/**
	 * Set the observations for the PDFs.
	 * Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated). 
	 * 
	 * @param observations
	 */
	public void setObservations(double observations1[][], double observations2[][]) throws Exception {
		mvke1.setObservations(observations1);
		mvke2.setObservations(observations2);
		// This call will throw an exception for us if the length of observations1 and 2 
		//  are not the same
		mvkeJoint.setObservations(observations1, observations2);
		totalObservations = observations1.length;
		this.observations1 = observations1;
		this.observations2 = observations2;
	}
	
	/**
	 * Compute the MI from the observations we were given
	 * 
	 * @return
	 */
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob1 = mvke1.getProbability(observations1[b], b);
			double prob2 = mvke2.getProbability(observations2[b], b);
			double probJoint = mvkeJoint.getProbability(observations1[b], observations2[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (probJoint > 0.0) {
				// If we have counted joint correlations, we must have marginals for each
				logTerm = probJoint / (prob1 * prob2);
				cont = Math.log(logTerm);
			}
			mi += cont;
			if (debug) {
				System.out.printf("%d: (%.5f, %.5f, %.5f) %.5f -> %.5f -> %.5f\n",
						b, prob1, prob2, probJoint, logTerm, cont, mi);
			}
		}
		lastAverage = mi / (double) totalObservations / Math.log(2.0);
		miComputed = true;
		return lastAverage;
	}

	/**
	 * Compute the MI if data were reordered.
	 * 
	 * @param newOrdering
	 * @return MI under the reordering scheme
	 */
	public double computeAverageLocalOfObservations(int[] newOrdering) throws Exception {
		// Store the real observations and their MI:
		double actualMI = lastAverage;
		double[][] originalData2 = observations2;
		double[][] data2;
		
		// Generate a new re-ordered data2
		data2 = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData2, newOrdering);
		observations2 = data2;
		// Perform new initialisations
		mvkeJoint.initialise(observations1[0].length + originalData2[0].length, epsilon);
		// Set new observations
		mvkeJoint.setObservations(observations1, data2);
		// Compute the MI
		double newMI = computeAverageLocalOfObservations();
		
		// Restore the actual MI and the observations
		lastAverage = actualMI;
		observations2 = originalData2;
		mvkeJoint.initialise(observations1[0].length + originalData2[0].length, epsilon);
		mvkeJoint.setObservations(observations1, originalData2);
		return newMI;
	}

	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param numPermutationsToCheck
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public synchronized EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(observations1.length, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}

	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param newOrderings the specific new orderings to use
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception {
		int numPermutationsToCheck = newOrderings.length;
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		// Store the real observations and their MI:
		double actualMI = lastAverage;
		double[][] originalData1 = observations1;
		double[][] originalData2 = observations2;
		double[][] data2;
		
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		
		int countWhereMiIsMoreSignificantThanOriginal = 0;
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Generate a new re-ordered data2
			data2 = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData2, newOrderings[i]);
			observations2 = data2;
			// Perform new initialisations
			mvkeJoint.initialise(originalData1[0].length + originalData2[0].length, epsilon);
			// Set new observations
			mvkeJoint.setObservations(originalData1, data2);
			// Compute the MI
			double newMI = computeAverageLocalOfObservations();
			measDistribution.distribution[i] = newMI;
			if (debug){
				System.out.println("New MI was " + newMI);
			}
			if (newMI >= actualMI) {
				countWhereMiIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the actual MI and the observations
		lastAverage = actualMI;
		observations2 = originalData2;
		mvkeJoint.initialise(originalData1[0].length + originalData2[0].length, epsilon);
		mvkeJoint.setObservations(originalData1, originalData2);
		// And return the significance
		measDistribution.pValue = (double) countWhereMiIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualMI;
		
		return measDistribution;
	}

	/**
	 * Extra utility method to return the joint entropy
	 * 
	 * @return
	 */
	public double computeAverageJointEntropy() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvkeJoint.getProbability(observations1[b], observations2[b], b);
			double cont = 0.0;
			if (prob > 0.0) {
				cont = - Math.log(prob);
			}
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + cont/Math.log(2.0) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * Extra utility method to return the entropy of the first set of joint variables
	 * 
	 * @return
	 */
	public double computeAverageEntropyOfObservation1() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvke1.getProbability(observations1[b], b);
			double cont = 0.0;
			// Comparing the prob to 0.0 should be fine - it would have to be
			//  an impossible number of samples for us to hit machine resolution here.
			if (prob > 0.0) {
				cont = -Math.log(prob);
			}
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + cont/Math.log(2.0) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * Extra utility method to return the entropy of the second set of joint variables
	 * 
	 * @return
	 */
	public double computeAverageEntropyOfObservation2() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvke2.getProbability(observations2[b], b);
			double cont = 0.0;
			if (prob > 0.0) {
				cont = -Math.log(prob);
			}
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + cont/Math.log(2.0) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * Extra utility method to return the information distance
	 * 
	 * @return
	 */
	public double computeAverageInfoDistanceOfObservations() {
		double infoDistance = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob1 = mvke1.getProbability(observations1[b], b);
			double prob2 = mvke2.getProbability(observations2[b], b);
			double probJoint = mvkeJoint.getProbability(observations1[b], observations2[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (probJoint > 0.0) {
				logTerm = (prob1 * prob2) / (probJoint * probJoint);
				cont = Math.log(logTerm);
			}
			infoDistance += cont;
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (infoDistance/Math.log(2.0)));
			}
		}
		return infoDistance / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * Compute the local MI values for the previous observations.
	 * 
	 * @return
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		return computeLocalUsingPreviousObservations(observations1, observations2, true);
	}


	/**
	 * Compute the local MI values for these given values, using the previously provided
	 *  observations to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 * 
	 * @param states1
	 * @param states2
	 * @return
	 */
	public double[] computeLocalUsingPreviousObservations(double states1[][], double states2[][]) {
		return computeLocalUsingPreviousObservations(states1, states2, false);
	}
	
	/**
	 * Internal method implementing local computation
	 * 
	 * @param states1
	 * @param states2
	 * @param isOurPreviousObservations
	 * @return
	 */
	protected double[] computeLocalUsingPreviousObservations(double states1[][],
			double states2[][], boolean isOurPreviousObservations) {
		
		double mi = 0.0;
		int timeSteps = states1.length;
		double[] localMi = new double[timeSteps];
		double prob1, prob2, probJoint;
		for (int b = 0; b < timeSteps; b++) {
			if (isOurPreviousObservations) {
				// We've been called with our previous observations, so we
				//  can pass the time step through for dynamic correlation exclusion
				prob1 = mvke1.getProbability(states1[b], b);
				prob2 = mvke2.getProbability(states2[b], b);
				probJoint = mvkeJoint.getProbability(states1[b], states2[b], b);
			} else {
				// We don't know whether these were our previous observation or not
				//  so we don't do dynamic correlation exclusion
				prob1 = mvke1.getProbability(states1[b]);
				prob2 = mvke2.getProbability(states2[b]);
				probJoint = mvkeJoint.getProbability(states1[b], states2[b]);
			}
			double logTerm = 0.0;
			localMi[b] = 0.0;
			if (probJoint > 0.0) {
				// By necessity prob1 and prob2 will be > 0.0
				logTerm = probJoint / (prob1 * prob2);
				localMi[b] = Math.log(logTerm) / Math.log(2.0);
			}
			mi += localMi[b];
			if (debug) {
				System.out.printf("%d: (%.5f, %.5f, %.5f) %.5f -> %.5f -> %.5f\n",
						b, prob1, prob2, probJoint, logTerm, localMi[b], mi);
			}
		}
		lastAverage = mi / (double) totalObservations;
		miComputed = true;
		return localMi;
	}

	
	/**
	 * Compute the local joint entropy values of the previously provided
	 *  observations.
	 *  
	 * @param states1
	 * @param states2
	 * @return
	 */
	public double[] computeLocalJointEntropyOfPreviousObservations() throws Exception {
		return computeLocalJointEntropyUsingPreviousObservations(observations1,
				observations2, true);
	}

	/**
	 * Compute the local joint entropy values for these given values, using the previously provided
	 *  observations to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 *  
	 * @param states1
	 * @param states2
	 * @return
	 */
	public double[] computeLocalJointEntropyUsingPreviousObservations(double states1[][], double states2[][]) {
		return computeLocalJointEntropyUsingPreviousObservations(states1,
				states2, false);
	}
	
	/**
	 * Internal implementation
	 * 
	 * @param states1
	 * @param states2
	 * @param isOurPreviousObservations
	 * @return
	 */
	private double[] computeLocalJointEntropyUsingPreviousObservations(
			double states1[][], double states2[][], boolean isOurPreviousObservations) {
		int timeSteps = states1.length;
		double[] localJoint = new double[timeSteps];
		double prob;
		for (int b = 0; b < totalObservations; b++) {
			if (isOurPreviousObservations) {
				prob = mvkeJoint.getProbability(observations1[b], observations2[b], b);
			} else {
				prob = mvkeJoint.getProbability(observations1[b], observations2[b]);
			}
			localJoint[b] = 0.0;
			if (prob > 0.0) {
				localJoint[b] = - Math.log(prob) / Math.log(2.0);
			}
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + localJoint[b]);
			}
		}
		return localJoint;
	}

	/**
	 * Compute the local entropy values for the previously provided
	 *  observations for VARIABLE 1 to compute the probabilities.
	 * 
	 * @param states1
	 *  
	 * @return
	 */
	public double[] computeLocalEntropy1OfPreviousObservations() {
		return computeLocalEntropyFromPreviousObservations(observations1, 1, true);
	}

	/**
	 * Compute the local entropy values for these given values, using the previously provided
	 *  observations for VARIABLE 1 to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 * 
	 * @param states1
	 *  
	 * @return
	 */
	public double[] computeLocalEntropy1UsingPreviousObservations(double[][] states) {
		return computeLocalEntropyFromPreviousObservations(states, 1, false);
	}

	/**
	 * Compute the local entropy values for the previously provided
	 *  observations for VARIABLE 1 to compute the probabilities.
	 * 
	 * @param states2
	 *  
	 * @return
	 */
	public double[] computeLocalEntropy2OfPreviousObservations() {
		return computeLocalEntropyFromPreviousObservations(observations2, 2, true);
	}

	/**
	 * Compute the local entropy values for these given values, using the previously provided
	 *  observations for VARIABLE 2 to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 *  
	 * @param states1
	 * @param states2
	 * @return
	 */
	public double[] computeLocalEntropy2UsingPreviousObservations(double states[][]) {
		return computeLocalEntropyFromPreviousObservations(states, 2, false);
	}

	/**
	 * Utility function to implement computeLocalEntropy1FromPreviousObservations
	 *  and computeLocalEntropy2FromPreviousObservations
	 *  
	 * @param states
	 * @param useProbsForWhichVar use 1 for variable 1, 2 for variable 2
	 * @param isOurPreviousObservations
	 * @return
	 */
	private double[] computeLocalEntropyFromPreviousObservations(
			double states[][], int useProbsForWhichVar, boolean isOurPreviousObservations) {
		int timeSteps = states.length;
		double[] localEntropy = new double[timeSteps];
		double prob;
		for (int b = 0; b < totalObservations; b++) {
			if (useProbsForWhichVar == 1) {
				if (isOurPreviousObservations) {
					prob = mvke1.getProbability(states[b], b);
				} else {
					prob = mvke1.getProbability(states[b]);
				}
			} else {
				if (isOurPreviousObservations) {
					prob = mvke2.getProbability(states[b], b);
				} else {
					prob = mvke2.getProbability(states[b]);
				}
			}
			localEntropy[b] = 0.0;
			if (prob > 0.0) {
				localEntropy[b] = - Math.log(prob) / Math.log(2.0);
			}
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + localEntropy[b]);
			}
		}
		return localEntropy;
	}

	/**
	 * Compute the local Info distance values for the previously provided
	 *  observations to compute the probabilities.
	 * 
	 * @return
	 */
	public double[] computeLocalInfoDistanceOfPreviousObservations() {
		return computeLocalInfoDistanceUsingPreviousObservations(observations1,
				observations2, true);
	}

	/**
	 * Compute the local Info distance values for these given values, using the previously provided
	 *  observations to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 * 
	 * @return
	 */
	public double[] computeLocalInfoDistanceUsingPreviousObservations(double[][] states1, double[][] states2) {
		return computeLocalInfoDistanceUsingPreviousObservations(states1,
				states2, false);
	}

	
	protected double[] computeLocalInfoDistanceUsingPreviousObservations(
			double[][] states1, double[][] states2, boolean isOurPreviousObservations) {
		int timeSteps = states1.length;
		double[] localInfoDistance = new double[timeSteps];
		double prob1, prob2, probJoint;
		for (int b = 0; b < timeSteps; b++) {
			if (isOurPreviousObservations) {
				prob1 = mvke1.getProbability(states1[b], b);
				prob2 = mvke2.getProbability(states2[b], b);
				probJoint = mvkeJoint.getProbability(states1[b], states2[b], b);
			} else {
				prob1 = mvke1.getProbability(states1[b]);
				prob2 = mvke2.getProbability(states2[b]);
				probJoint = mvkeJoint.getProbability(states1[b], states2[b]);
			}
			double logTerm = 0.0;
			localInfoDistance[b] = 0.0;
			if (probJoint > 0.0) {
				logTerm = (prob1 * prob2) / (probJoint * probJoint);
				localInfoDistance[b] = Math.log(logTerm) / Math.log(2.0);
			}
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + localInfoDistance[b]);
			}
		}
		return localInfoDistance;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return lastAverage;
	}

	/**
	 * Set properties for the mutual information calculator.
	 * These can include:
	 * <ul>
	 * 		<li>{@link #EPSILON_PROP_NAME}</li>
	 * 		<li>{@link #NORMALISE_PROP_NAME}</li>
	 * 		<li>{@link #DYN_CORR_EXCL_TIME_NAME}</li>
	 * 		<li>{@link #FORCE_KERNEL_COMPARE_TO_ALL}</li>
	 * </ul>
	 * 
	 * Note that dynamic correlation exclusion may have unexpected results if multiple
	 *  observation sets have been added. This is because multiple observation sets
	 *  are treated as though they are from a single time series, so observations from
	 *  near the end of observation set i will be excluded from comparison to 
	 *  observations near the beginning of observation set (i+1). 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 */
	public void setProperty(String propertyName, String propertyValue) {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			epsilon = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
			mvke1.setNormalise(normalise);
			mvke2.setNormalise(normalise);
			mvkeJoint.setNormalise(normalise);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
			if (dynCorrExcl) {
				mvke1.setDynamicCorrelationExclusion(dynCorrExclTime);
				mvke2.setDynamicCorrelationExclusion(dynCorrExclTime);
				mvkeJoint.setDynamicCorrelationExclusion(dynCorrExclTime);
			} else {
				mvke1.clearDynamicCorrelationExclusion();
				mvke2.clearDynamicCorrelationExclusion();
				mvkeJoint.clearDynamicCorrelationExclusion();
			}
		} else if (propertyName.equalsIgnoreCase(FORCE_KERNEL_COMPARE_TO_ALL)) {
			forceCompareToAll = Boolean.parseBoolean(propertyValue);
			mvke1.setForceCompareToAll(forceCompareToAll);
			mvke2.setForceCompareToAll(forceCompareToAll);
			mvkeJoint.setForceCompareToAll(forceCompareToAll);
		} else if (propertyName.equalsIgnoreCase(PROP_TIME_DIFF)) {
			int diff = Integer.parseInt(propertyValue);
			if (diff != 0) {
				throw new RuntimeException(PROP_TIME_DIFF + " property != 0 not implemented yet");
			}
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	public int getNumObservations() {
		return totalObservations;
	}
}
