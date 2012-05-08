package infodynamics.measures.continuous.kernel;

import infodynamics.measures.continuous.MultiInfoCalculator;
import infodynamics.utils.MatrixUtils;

import java.util.Vector;
import java.util.Random;


public class MultiInfoCalculatorKernel implements
	MultiInfoCalculator {

	KernelEstimatorSingleVariate[] svkeMarginals = null;
	KernelEstimatorMultiVariate mvkeJoint = null;

	private int dimensions = 0;
	private int totalObservations = 0;
	private boolean debug = false;
	private double[][] observations;
	private Vector<double[]> individualObservations;
	private double lastAverage;
	
	private boolean normalise = true;
	public static final String NORMALISE_PROP_NAME = "NORMALISE";

	private boolean dynCorrExcl = false;
	private int dynCorrExclTime = 100;
	public static final String DYN_CORR_EXCL_TIME_NAME = "DYN_CORR_EXCL";

	private boolean underSample = false;
	private double samplingFactor = 0.1;
	private Random rand;
	public static final String SAMPLING_FACTOR_PROP_NAME = "SAMPLING_FACTOR";

	/**
	 * Default value for epsilon
	 */
	public static final double DEFAULT_EPSILON = 0.25;
	private double epsilon = DEFAULT_EPSILON;
	public static final String EPSILON_PROP_NAME = "EPSILON";

	public MultiInfoCalculatorKernel() {
		mvkeJoint = new KernelEstimatorMultiVariate();
	}

	/**
	 * Initialise using a default epsilon
	 * 
	 * @param dimensions1
	 */
	public void initialise(int dimensions) {
		initialise(dimensions, epsilon);
	}

	public void initialise(int dimensions, double epsilon) {
		this.epsilon = epsilon;
		observations = null;
		if (this.dimensions != dimensions) {
			// Need to create a new array of marginal kernel estimators 
			this.dimensions = dimensions;
			svkeMarginals = new KernelEstimatorSingleVariate[dimensions];
			for (int i = 0; i < dimensions; i++) {
				svkeMarginals[i] = new KernelEstimatorSingleVariate();
				svkeMarginals[i].setNormalise(normalise);
				if (dynCorrExcl) {
					svkeMarginals[i].setDynamicCorrelationExclusion(dynCorrExclTime);
				} else {
					svkeMarginals[i].clearDynamicCorrelationExclusion();
				}
			}
		}
		// Initialise the marginal kernel estimators
		for (int i = 0; i < dimensions; i++) {
			svkeMarginals[i].initialise(epsilon);
		}
		// Initialise the joint kernel estimator
		mvkeJoint.initialise(dimensions, epsilon);
		lastAverage = 0.0;
	}

	/**
	 * Set the observations for the PDFs.
	 * Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated). 
	 * 
	 * @param observations
	 */
	public void setObservations(double observations[][]) throws Exception {
		if (observations[0].length != dimensions) {
			throw new Exception("Incorrect number of dimensions " + observations[0].length +
					" in supplied observations (expected " + dimensions + ")");
		}
		for (int d = 0; d < dimensions; d++) {
			svkeMarginals[d].setObservations(MatrixUtils.selectColumn(observations, d));
		}
		mvkeJoint.setObservations(observations);
		totalObservations = observations.length;
		this.observations = observations;
	}
	
	/**
	 * User elects to set observations one by one rather than in one go.
	 * Will need to call endIndividualObservations before calling any of the
	 *  compute functions, otherwise the previous observations will be used.
	 */
	public void startIndividualObservations() {
		if (dynCorrExcl) {
			// We have not properly implemented dynamic correlation exclusion for
			//  multiple observation sets, so throw an error
			throw new RuntimeException("Addition of multiple observation sets is not currently " +
					"supported with property DYN_CORR_EXCL set");
		}
		individualObservations = new Vector<double[]>();
	}
	
	public void addObservation(double observation[]) {
		if (underSample && (rand.nextDouble() >= samplingFactor)) {
			// Don't take this sample
			return;
		}
		individualObservations.add(observation);
	}
	
	public void endIndividualObservations() throws Exception {
		double[][] data = new double[individualObservations.size()][];
		for (int t = 0; t < data.length; t++) {
			data[t] = individualObservations.elementAt(t);
		}
		// Allow vector to be reclaimed
		individualObservations = null;
		setObservations(data);
	}

	/**
	 * Compute the MI from the observations we were given
	 * 
	 * @return
	 */
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			if (debug) {
				System.out.print(b + ": ");
			}
			double marginalProbProducts = 1.0;
			for (int d = 0; d < dimensions; d++) {
				double marginalProb = svkeMarginals[d].getProbability(observations[b][d], b);
				if (debug) {
					System.out.print(observations[b][d] + " p=" + marginalProb + ", ");
				}
				marginalProbProducts *= marginalProb;
			}
			double probJoint = mvkeJoint.getCount(observations[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (probJoint > 0.0) {
				// TODO Should probably check that marginalProbProducts has not 
				//  gone to zero (this is possible with several multiplications
				//  of 1/N, though is unlikely). I'm not sure what we would do if
				//  it had gone to zero though ... ignore this value?
				logTerm = probJoint / marginalProbProducts;
				cont = Math.log(logTerm);
			}
			mi += cont;
			if (debug) {
				System.out.println(", p(joint) = " + probJoint
						+ " -> " + logTerm + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (mi/Math.log(2.0)));
			}
		}
		lastAverage = mi / (double) totalObservations / Math.log(2.0);
		return lastAverage;
	}

	/**
	 * Extra utility method to return the joint entropy
	 * 
	 * @return
	 */
	public double computeAverageJointEntropy() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvkeJoint.getCount(observations[b], b);
			double cont = 0.0;
			if (prob > 0.0) {
				cont = - Math.log(prob);
			}
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob
						+ " -> " + cont/Math.log(2.0) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * Extra utility method to return the entropy of the first set of joint variables
	 * 
	 * @return
	 */
	public double computeAverageMarginalEntropy(int variableIndex) {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = svkeMarginals[variableIndex].getProbability(observations[b][variableIndex], b);
			double cont = 0.0;
			if (prob > 0.0) {
				cont = -Math.log(prob);
			}
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob
						+ " -> " + cont/Math.log(2.0) + " -> sum: " + (entropy/Math.log(2.0)));
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
			double marginalProbProducts = 1.0;
			for (int d = 0; d < dimensions; d++) {
				marginalProbProducts *= svkeMarginals[d].getProbability(observations[b][d], b);
			}
			double probJoint = mvkeJoint.getProbability(observations[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (probJoint > 0.0) {
				// TODO Should probably check that marginalProbProducts has not 
				//  gone to zero (this is possible with several multiplications
				//  of 1/N, though is unlikely). I'm not sure what we would do if
				//  it had gone to zero though ... ignore this value?
				//  It's not easy to ignore since we can't say p log p -> 0 here
				//  because marginalProbProducts is not multiplying out the front
				logTerm = marginalProbProducts / (probJoint * probJoint);
				cont = Math.log(logTerm);
			}
			infoDistance += cont;
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (cont/Math.log(2.0)) +
						" -> sum: " + (infoDistance/Math.log(2.0)));
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
		return computeLocalUsingPreviousObservations(observations, true);
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
	public double[] computeLocalUsingPreviousObservations(double states[][]) {
		return computeLocalUsingPreviousObservations(states, false);
	}
	
	/**
	 * Internal implementation
	 */
	protected double[] computeLocalUsingPreviousObservations(double states[][],
			boolean isOurPreviousObservations) {
		double mi = 0.0;
		int timeSteps = states.length;
		double[] localMi = new double[timeSteps];
		double probJoint;
		for (int b = 0; b < timeSteps; b++) {
			double marginalProbProducts = 1.0;
			for (int d = 0; d < dimensions; d++) {
				if (isOurPreviousObservations) {
					marginalProbProducts *= svkeMarginals[d].getProbability(states[b][d], b);
				} else {
					marginalProbProducts *= svkeMarginals[d].getProbability(states[b][d]);
				}
			}
			if (isOurPreviousObservations) {
				probJoint = mvkeJoint.getProbability(states[b], b);
			} else {
				probJoint = mvkeJoint.getProbability(states[b]);
			}
			double logTerm = 0.0;
			localMi[b] = 0.0;
			if (probJoint > 0.0) {
				// TODO Should probably check that marginalProbProducts has not 
				//  gone to zero (this is possible with several multiplications
				//  of 1/N, though is unlikely). I'm not sure what we would do if
				//  it had gone to zero though ... ignore this value?
				logTerm = probJoint / marginalProbProducts;
				localMi[b] = Math.log(logTerm) / Math.log(2.0);
			}
			mi += localMi[b];
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + localMi[b] + " -> sum: " + mi);
			}
		}
		lastAverage = mi / (double) totalObservations;
		return localMi;
	}

	
	public double[] computeLocalJointEntropyOfPreviousObservations() throws Exception {
		return computeLocalJointEntropyUsingPreviousObservations(observations, true);
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
	public double[] computeLocalJointEntropyUsingPreviousObservations(double states[][]) {
		return computeLocalJointEntropyUsingPreviousObservations(states, false);
	}
	
	/**
	 * Internal implementation
	 * 
	 * @param states
	 * @param isOurPreviousObservations
	 * @return
	 */
	protected double[] computeLocalJointEntropyUsingPreviousObservations(
			double states[][], boolean isOurPreviousObservations) {

		int timeSteps = states.length;
		double[] localJoint = new double[timeSteps];
		double probJoint;
		for (int b = 0; b < totalObservations; b++) {
			if (isOurPreviousObservations) {
				probJoint = mvkeJoint.getProbability(states[b], b);
			} else {
				probJoint = mvkeJoint.getProbability(states[b]);
			}
			localJoint[b] = 0.0;
			if (probJoint > 0.0) {
				localJoint[b] = - Math.log(probJoint) / Math.log(2.0);
			}
			if (debug) {
				System.out.println(b + ": " + probJoint + " -> " + localJoint[b]);
			}
		}
		return localJoint;
	}

	
	public double[] computeLocalMarginalEntropyOfPreviousObservations(int variableIndex) {
		return computeLocalMarginalEntropyUsingPreviousObservations(observations, variableIndex, true);
	}

	/**
	 * 
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 *  
	 * @param states
	 * @param useProbsForWhichVar use 1 for variable 1, 2 for variable 2
	 * @return
	 */
	public double[] computeLocalMarginalEntropyUsingPreviousObservations(double states[][], int variableIndex) {
		return computeLocalMarginalEntropyUsingPreviousObservations(states, variableIndex, false);
	}
	
	/**
	 * Internal implementation
	 * 
	 * @param states
	 * @param variableIndex
	 * @param isOurPreviousObservations
	 * @return
	 */
	protected double[] computeLocalMarginalEntropyUsingPreviousObservations(
			double states[][], int variableIndex, boolean isOurPreviousObservations) {
		int timeSteps = states.length;
		double[] localEntropy = new double[timeSteps];
		double prob;
		for (int b = 0; b < totalObservations; b++) {
			if (isOurPreviousObservations) {
				prob = svkeMarginals[variableIndex].getProbability(states[b][variableIndex], b);
			} else {
				prob = svkeMarginals[variableIndex].getProbability(states[b][variableIndex]);
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
		return computeLocalInfoDistanceUsingPreviousObservations(observations, true);
	}

	/**
	 * Compute the local Info distance values for these given values, using the previously provided
	 *  observations to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 * 
	 * @return
	 */
	public double[] computeLocalInfoDistanceUsingPreviousObservations(double[][] states) {
		return computeLocalInfoDistanceUsingPreviousObservations(states, false);
	}

	/**
	 * Internal implementation
	 * 
	 * @param states
	 * @param isOurPreviousObservations
	 * @return
	 */
	protected double[] computeLocalInfoDistanceUsingPreviousObservations(
			double[][] states, boolean isOurPreviousObservations) {
		
		int timeSteps = states.length;
		double[] localInfoDistance = new double[timeSteps];
		double probJoint;
		for (int b = 0; b < timeSteps; b++) {
			double marginalProbProducts = 1.0;
			for (int d = 0; d < dimensions; d++) {
				if (isOurPreviousObservations) {
					marginalProbProducts *= svkeMarginals[d].getProbability(states[b][d], b);
				} else {
					marginalProbProducts *= svkeMarginals[d].getProbability(states[b][d]);
				}
			}
			if (isOurPreviousObservations) {
				probJoint = mvkeJoint.getProbability(states[b], b);
			} else {
				probJoint = mvkeJoint.getProbability(states[b]);
			}
			double logTerm = 0.0;
			localInfoDistance[b] = 0.0;
			if (probJoint > 0.0) {
				// TODO Should probably check that marginalProbProducts has not 
				//  gone to zero (this is possible with several multiplications
				//  of 1/N, though is unlikely). I'm not sure what we would do if
				//  it had gone to zero though ... ignore this value?
				//  It's not easy to ignore since we can't say p log p -> 0 here
				//  because marginalProbProducts is not multiplying out the front
				logTerm = marginalProbProducts / (probJoint * probJoint);
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

	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			epsilon = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
			for (int d = 0; d < dimensions; d++) {
				svkeMarginals[d].setNormalise(normalise);
			}
			mvkeJoint.setNormalise(normalise);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
			if (dynCorrExcl) {
				for (int d = 0; d < dimensions; d++) {
					svkeMarginals[d].setDynamicCorrelationExclusion(dynCorrExclTime);
				}
				mvkeJoint.setDynamicCorrelationExclusion(dynCorrExclTime);
			} else {
				for (int d = 0; d < dimensions; d++) {
					svkeMarginals[d].clearDynamicCorrelationExclusion();
				}
				mvkeJoint.clearDynamicCorrelationExclusion();
			}
		} else if (propertyName.equalsIgnoreCase(SAMPLING_FACTOR_PROP_NAME)) {
			// Use less than 100 % of samples in the addObservation method
			samplingFactor = Double.parseDouble(propertyValue);
			underSample = (samplingFactor < 1.0);
			rand = new Random();
		} else {
			// Property not recognised
			return;
		}
		if (debug) {
			System.out.println("Set property " + propertyName +
					" to " + propertyValue);
		}
	}
}
