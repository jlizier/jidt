package infodynamics.measures.spiking.integration;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Vector;
import java.util.ArrayList;

//import infodynamics.measures.continuous.kraskov.EuclideanUtils;
import infodynamics.measures.spiking.TransferEntropyCalculatorSpiking;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.KdTree;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.NeighbourNodeData;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.UnivariateNearestNeighbourSearcher;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.ParsedProperties;

/**
 * Computes the transfer entropy between a pair of spike trains, using an
 * integration-based measure in order to match the theoretical form of TE
 * between such spike trains.
 * 
 * <p>
 * Usage paradigm is as per the interface
 * {@link TransferEntropyCalculatorSpiking}
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 *         <a href="http://lizier.me/joseph/">www</a>)
 */

/*
 * TODO
 * This implementation of the estimator does not implement dynamic exclusion windows. Such windows make sure
 * that history embeddings that overlap are not considered in nearest-neighbour searches (as this breaks the
 * independece assumption). Getting dynamic exclusion windows working will probably require modifications to the
 * KdTree class.
 */
public class TransferEntropyCalculatorSpikingIntegration implements TransferEntropyCalculatorSpiking {

	/**
	 * Number of past destination interspike intervals to consider (akin to embedding length)
	 */
	protected int k = 1;
	/**
	 * Number of past source interspike intervals to consider (akin to embedding length)
	 */
	protected int l = 1;

	/**
	 * Property name for number of interspike intervals for the conditional variables
	 */
	public static final String COND_EMBED_LENGTHS_PROP_NAME = "COND_EMBED_LENGTHS";
	/**
	 * Array of history interspike interval embedding lengths for the conditional variables.
	 *  Can be an empty array or null if there are no conditional variables.
	 */
	protected int[] condEmbedDims = new int[] {};

	/**
	 * Number of nearest neighbours to search for in the full joint space
	 */
	protected int Knns = 4;

	/**
	 * Storage for source observations supplied via
	 * {@link #addObservations(double[], double[])} etc.
	 */
	protected Vector<double[]> vectorOfSourceSpikeTimes = null;

	/**
	 * Storage for destination observations supplied via
	 * {@link #addObservations(double[], double[])} etc.
	 */
	protected Vector<double[]> vectorOfDestinationSpikeTimes = null;

	/**
	 * Storage for conditional observations supplied via
	 * {@link #addObservations(double[], double[])} etc.
	 */
	protected Vector<double[][]> vectorOfConditionalSpikeTimes = null;

	Vector<double[]> conditioningEmbeddingsFromSpikes = null;
	Vector<double[]> jointEmbeddingsFromSpikes = null;
	Vector<double[]> conditioningEmbeddingsFromSamples = null;
	Vector<double[]> jointEmbeddingsFromSamples = null;
	Vector<Double> processTimeLengths = null;

	protected KdTree kdTreeJointAtSpikes = null;
	protected KdTree kdTreeJointAtSamples = null;
	protected KdTree kdTreeConditioningAtSpikes = null;
	protected KdTree kdTreeConditioningAtSamples = null;

	public static final String KNNS_PROP_NAME = "Knns";

	/**
	 * Property name for an amount of random Gaussian noise to be added to the data
	 * (default is 1e-8, matching the MILCA toolkit).
	 */
	public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";	
	
	/**
	 * Whether to add an amount of random noise to the incoming data
	 */
	protected boolean addNoise = true;
	/**
	 * Amount of random Gaussian noise to add to the incoming data
	 */
	protected double noiseLevel = (double) 1e-8;

	/**
	 * Stores whether we are in debug mode
	 */
	protected boolean debug = false;

	/**
	 * Property name for the number of random sample points to use as a multiple
	 * of the number of target spikes.
	 */
	public static final String PROP_SAMPLE_MULTIPLIER = "NUM_SAMPLES_MULTIPLIER";
	protected double numSamplesMultiplier = 1.0;
	/**
	 * Property name for the number of random sample points to use in the construction of the surrogates as a multiple
	 * of the number of target spikes.
	 */
	public static final String PROP_SURROGATE_SAMPLE_MULTIPLIER = "SURROGATE_NUM_SAMPLES_MULTIPLIER";
	protected double surrogateNumSamplesMultiplier = 1.0;

	/**
	 * Property for the number of nearest neighbours to use in the construction of the surrogates
	 */
	public static final String PROP_K_PERM = "K_PERM";
	protected int kPerm = 10;

	
	/**
	 * Property name for what type of norm to use between data points
	 *  for each marginal variable -- Options are defined by 
	 *  {@link KdTree#setNormType(String)} and the
	 *  default is {@link EuclideanUtils#NORM_EUCLIDEAN}.
	 */
	public final static String PROP_NORM_TYPE = "NORM_TYPE";
	protected int normType = EuclideanUtils.NORM_EUCLIDEAN;
	
	public TransferEntropyCalculatorSpikingIntegration() {
		super();
	}
	

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#initialise(
	 * int)
	 */
	@Override
	public void initialise() throws Exception {
		initialise(k, l);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#initialise(
	 * int)
	 */
	@Override
	public void initialise(int k) throws Exception {
		initialise(k, this.l);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#initialise(
	 * int, int)
	 */
	@Override
	public void initialise(int k, int l) throws Exception {
		if ((k < 1) || (l < 1)) {
			throw new Exception("Zero history length not supported");
		}
		this.k = k;
		this.l = l;
		vectorOfSourceSpikeTimes = null;
		vectorOfDestinationSpikeTimes = null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#setProperty(
	 * java.lang.String, java.lang.String)
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			int kTemp = Integer.parseInt(propertyValue);
			if (kTemp < 1) {
				throw new Exception ("Invalid k value less than 1.");
			} else {
				k = kTemp;
			}
		} else if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			int lTemp = Integer.parseInt(propertyValue);
			if (lTemp < 1) {
				throw new Exception ("Invalid l value less than 1.");
			} else {
				l = lTemp;
			}
		} else if (propertyName.equalsIgnoreCase(COND_EMBED_LENGTHS_PROP_NAME)) {
			int[] condEmbedDimsTemp = ParsedProperties.parseStringArrayOfInts(propertyValue);
			for (int dim : condEmbedDimsTemp) {
				if (dim < 1) {
					throw new Exception ("Invalid conditional embedding value less than 1.");
				}
			}
			condEmbedDims = condEmbedDimsTemp;
		} else if (propertyName.equalsIgnoreCase(KNNS_PROP_NAME)) {
			Knns = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_K_PERM)) {
			kPerm = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
			if (propertyValue.equals("0") || propertyValue.equalsIgnoreCase("false")) {
				addNoise = false;
				noiseLevel = 0;
			} else {
				addNoise = true;
				noiseLevel = Double.parseDouble(propertyValue);
			}
			
		} else if (propertyName.equalsIgnoreCase(PROP_SAMPLE_MULTIPLIER)) {
			double tempNumSamplesMultiplier = Double.parseDouble(propertyValue);
			if (tempNumSamplesMultiplier <= 0) {
				throw new Exception ("Num samples multiplier must be greater than 0.");
			} else {
				numSamplesMultiplier = tempNumSamplesMultiplier;
			}
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normType = KdTree.validateNormType(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_SURROGATE_SAMPLE_MULTIPLIER)) {
			double tempSurrogateNumSamplesMultiplier = Double.parseDouble(propertyValue);
			if (tempSurrogateNumSamplesMultiplier <= 0) {
				throw new Exception ("Surrogate Num samples multiplier must be greater than 0.");
			} else {
				surrogateNumSamplesMultiplier = tempSurrogateNumSamplesMultiplier;
			}
		} else {
			// No property was set on this class
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(
					this.getClass().getSimpleName() + ": Set property " + propertyName + " to " + propertyValue);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#getProperty(
	 * java.lang.String)
	 */
	@Override
	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			return Integer.toString(k);
		} else if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			return Integer.toString(l);
		} else if (propertyName.equalsIgnoreCase(KNNS_PROP_NAME)) {
			return Integer.toString(Knns);
		} else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
			return Double.toString(noiseLevel);
		} else if (propertyName.equalsIgnoreCase(PROP_SAMPLE_MULTIPLIER)) {
			return Double.toString(numSamplesMultiplier);
		} else {
			// No property matches for this class
			return null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#
	 * setObservations(double[], double[])
	 */
	@Override
	public void setObservations(double[] source, double[] destination) throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}

	public void setObservations(double[] source, double[] destination, double[][] conditionals) throws Exception {
		startAddObservations();
		addObservations(source, destination, conditionals);
		finaliseAddObservations();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#
	 * startAddObservations()
	 */
	@Override
	public void startAddObservations() {
		vectorOfSourceSpikeTimes = new Vector<double[]>();
		vectorOfDestinationSpikeTimes = new Vector<double[]>();
		vectorOfConditionalSpikeTimes = new Vector<double[][]>();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#
	 * addObservations(double[], double[])
	 */
	@Override
	public void addObservations(double[] source, double[] destination) throws Exception {
		// Store these observations in our vector for now
		vectorOfSourceSpikeTimes.add(source);
		vectorOfDestinationSpikeTimes.add(destination);
	}

	public void addObservations(double[] source, double[] destination, double[][] conditionals) throws Exception {
		// Store these observations in our vector for now
		vectorOfSourceSpikeTimes.add(source);
		vectorOfDestinationSpikeTimes.add(destination);
		vectorOfConditionalSpikeTimes.add(conditionals);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#
	 * finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {

		conditioningEmbeddingsFromSpikes = new Vector<double[]>();
		jointEmbeddingsFromSpikes = new Vector<double[]>();
		conditioningEmbeddingsFromSamples = new Vector<double[]>();
		jointEmbeddingsFromSamples = new Vector<double[]>();
		processTimeLengths = new Vector<Double>();

		// Send all of the observations through:
		Iterator<double[]> sourceIterator = vectorOfSourceSpikeTimes.iterator();
		Iterator<double[][]> conditionalIterator = null;
		if (vectorOfConditionalSpikeTimes.size() > 0) {
			conditionalIterator = vectorOfConditionalSpikeTimes.iterator();
		}
		int timeSeriesIndex = 0;
		for (double[] destSpikeTimes : vectorOfDestinationSpikeTimes) {
			double[] sourceSpikeTimes = sourceIterator.next();
			double[][] conditionalSpikeTimes = null;
			if (vectorOfConditionalSpikeTimes.size() > 0) {
				conditionalSpikeTimes = conditionalIterator.next();
			} else {
				conditionalSpikeTimes = new double[][] {};
			}
			processEventsFromSpikingTimeSeries(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes, conditioningEmbeddingsFromSpikes,
							   jointEmbeddingsFromSpikes, conditioningEmbeddingsFromSamples, jointEmbeddingsFromSamples,
							   numSamplesMultiplier);
		}

		// Convert the vectors to arrays so that they can be put in the trees
		double[][] arrayedTargetEmbeddingsFromSpikes = new double[conditioningEmbeddingsFromSpikes.size()][k];
		double[][] arrayedJointEmbeddingsFromSpikes = new double[conditioningEmbeddingsFromSpikes.size()][k + l];
		for (int i = 0; i < conditioningEmbeddingsFromSpikes.size(); i++) {
			arrayedTargetEmbeddingsFromSpikes[i] = conditioningEmbeddingsFromSpikes.elementAt(i);
			arrayedJointEmbeddingsFromSpikes[i] = jointEmbeddingsFromSpikes.elementAt(i);
		}
		double[][] arrayedTargetEmbeddingsFromSamples = new double[conditioningEmbeddingsFromSamples.size()][k];
		double[][] arrayedJointEmbeddingsFromSamples = new double[conditioningEmbeddingsFromSamples.size()][k + l];
		for (int i = 0; i < conditioningEmbeddingsFromSamples.size(); i++) {
			arrayedTargetEmbeddingsFromSamples[i] = conditioningEmbeddingsFromSamples.elementAt(i);
			arrayedJointEmbeddingsFromSamples[i] = jointEmbeddingsFromSamples.elementAt(i);
		}

		kdTreeJointAtSpikes = new KdTree(arrayedJointEmbeddingsFromSpikes);
		kdTreeJointAtSamples = new KdTree(arrayedJointEmbeddingsFromSamples);
		kdTreeConditioningAtSpikes = new KdTree(arrayedTargetEmbeddingsFromSpikes);
		kdTreeConditioningAtSamples = new KdTree(arrayedTargetEmbeddingsFromSamples);

		kdTreeJointAtSpikes.setNormType(normType);
		kdTreeJointAtSamples.setNormType(normType);
		kdTreeConditioningAtSpikes.setNormType(normType);
		kdTreeConditioningAtSamples.setNormType(normType);
	}

	protected void makeEmbeddingsAtPoints(double[] pointsAtWhichToMakeEmbeddings, int indexOfFirstPointToUse,
					      double[] sourceSpikeTimes, double[] destSpikeTimes,
					      double[][] conditionalSpikeTimes,
					      Vector<double[]> conditioningEmbeddings, 
					      Vector<double[]> jointEmbeddings) {

		Random random = new Random();

		int embeddingPointIndex = indexOfFirstPointToUse;
		int mostRecentDestIndex = k;
		int mostRecentSourceIndex = l;
		int[] mostRecentConditioningIndices = Arrays.copyOf(condEmbedDims, condEmbedDims.length);
		int totalLengthOfConditioningEmbeddings = 0;
		for (int i = 0; i < condEmbedDims.length; i++) {
			totalLengthOfConditioningEmbeddings += condEmbedDims[i];
		}

		// Loop through the points at which embeddings need to be made
		for (; embeddingPointIndex < pointsAtWhichToMakeEmbeddings.length; embeddingPointIndex++) {

			// Advance the tracker of the most recent dest index
			while (mostRecentDestIndex < (destSpikeTimes.length - 1)) {
				if (destSpikeTimes[mostRecentDestIndex + 1] < pointsAtWhichToMakeEmbeddings[embeddingPointIndex]) {
					mostRecentDestIndex++;
				} else {
					break;
				}
			}
			// Do the same for the most recent source index
			while (mostRecentSourceIndex < (sourceSpikeTimes.length - 1)) {
				if (sourceSpikeTimes[mostRecentSourceIndex
						+ 1] < pointsAtWhichToMakeEmbeddings[embeddingPointIndex]) {
					mostRecentSourceIndex++;
				} else {
					break;
				}
			}
			// Now advance the trackers for the most recent conditioning indices
			for (int j = 0; j < mostRecentConditioningIndices.length; j++) {
				while (mostRecentConditioningIndices[j] < (conditionalSpikeTimes[j].length - 1)) {
					if (conditionalSpikeTimes[j][mostRecentConditioningIndices[j] + 1] < pointsAtWhichToMakeEmbeddings[embeddingPointIndex]) {
						mostRecentConditioningIndices[j]++;
					} else {
						break;
					}
				}
			}

			
			double[] conditioningPast = new double[k + totalLengthOfConditioningEmbeddings];
			double[] jointPast = new double[k + totalLengthOfConditioningEmbeddings + l];

			// Add the embedding intervals from the target process
			conditioningPast[0] = pointsAtWhichToMakeEmbeddings[embeddingPointIndex] - destSpikeTimes[mostRecentDestIndex];
			jointPast[0] = pointsAtWhichToMakeEmbeddings[embeddingPointIndex]
					- destSpikeTimes[mostRecentDestIndex];
			for (int i = 1; i < k; i++) {
				conditioningPast[i] = destSpikeTimes[mostRecentDestIndex - i + 1]
						- destSpikeTimes[mostRecentDestIndex - i];
				jointPast[i] = destSpikeTimes[mostRecentDestIndex - i + 1]
						- destSpikeTimes[mostRecentDestIndex - i];
			}

			// Add the embeding intervals from the conditional processes
			int indexOfNextEmbeddingInterval = k;
			for (int i = 0; i < condEmbedDims.length; i++) {
				conditioningPast[indexOfNextEmbeddingInterval] =
					pointsAtWhichToMakeEmbeddings[embeddingPointIndex] - conditionalSpikeTimes[i][mostRecentConditioningIndices[i]]; 
				jointPast[indexOfNextEmbeddingInterval] = 
					pointsAtWhichToMakeEmbeddings[embeddingPointIndex] - conditionalSpikeTimes[i][mostRecentConditioningIndices[i]];
				indexOfNextEmbeddingInterval += 1;
				for (int j = 1; j < condEmbedDims[i]; j++) {
					conditioningPast[indexOfNextEmbeddingInterval] =
						conditionalSpikeTimes[i][mostRecentConditioningIndices[i] - j + 1] -
						conditionalSpikeTimes[i][mostRecentConditioningIndices[i] - j]; 
					jointPast[indexOfNextEmbeddingInterval] =
						conditionalSpikeTimes[i][mostRecentConditioningIndices[i] - j + 1] -
						conditionalSpikeTimes[i][mostRecentConditioningIndices[i] - j]; 
					indexOfNextEmbeddingInterval += 1;
				}
			}

			// Add the embedding intervals from the source process (this only gets added to the joint embeddings)
			jointPast[k + totalLengthOfConditioningEmbeddings] = pointsAtWhichToMakeEmbeddings[embeddingPointIndex]
					- sourceSpikeTimes[mostRecentSourceIndex];
			for (int i = 1; i < l; i++) {
				jointPast[k + totalLengthOfConditioningEmbeddings + i] = sourceSpikeTimes[mostRecentSourceIndex - i + 1]
						- sourceSpikeTimes[mostRecentSourceIndex - i];
			}

			// Add Gaussian noise, if necessary
			if (addNoise) {
				for (int i = 0; i < conditioningPast.length; i++) {
					conditioningPast[i] += random.nextGaussian() * noiseLevel;
				}
				for (int i = 0; i < jointPast.length; i++) {
					jointPast[i] += random.nextGaussian() * noiseLevel;
				}
			}

			conditioningEmbeddings.add(conditioningPast);
			jointEmbeddings.add(jointPast);
		}
	}



	protected int getFirstDestIndex(double[] sourceSpikeTimes, double[] destSpikeTimes, double[][] conditionalSpikeTimes, Boolean setProcessTimeLengths)
		throws Exception{

		// First sort the spike times in case they were not properly in ascending order:
		Arrays.sort(sourceSpikeTimes);
		Arrays.sort(destSpikeTimes);
		
		int firstTargetIndexOfEMbedding = k;
		while (destSpikeTimes[firstTargetIndexOfEMbedding] <= sourceSpikeTimes[l - 1]) {
			firstTargetIndexOfEMbedding++;
		}
		if (conditionalSpikeTimes.length != condEmbedDims.length) {
			throw new Exception("Number of conditional embedding lengths does not match the number of conditional processes");
		}
		for (int i = 0; i < conditionalSpikeTimes.length; i++) {
			while (destSpikeTimes[firstTargetIndexOfEMbedding] <= conditionalSpikeTimes[i][condEmbedDims[i]]) {
				firstTargetIndexOfEMbedding++;
			}
		}

		// We don't want to reset these lengths when resampling for surrogates
		if (setProcessTimeLengths) {
			processTimeLengths.add(destSpikeTimes[destSpikeTimes.length - 1] - destSpikeTimes[firstTargetIndexOfEMbedding]);
		}
		
		return firstTargetIndexOfEMbedding;
	}

	protected double[] generateRandomSampleTimes(double[] sourceSpikeTimes, double[] destSpikeTimes, double[][] conditionalSpikeTimes,
						   double actualNumSamplesMultiplier, int firstTargetIndexOfEMbedding) {

		double sampleLowerBound = destSpikeTimes[firstTargetIndexOfEMbedding];
		double sampleUpperBound = destSpikeTimes[destSpikeTimes.length - 1];
		int num_samples = (int) Math.round(actualNumSamplesMultiplier * (destSpikeTimes.length - firstTargetIndexOfEMbedding + 1));
		double[] randomSampleTimes = new double[num_samples];
		Random rand = new Random();
		for (int i = 0; i < randomSampleTimes.length; i++) {
			randomSampleTimes[i] = sampleLowerBound + rand.nextDouble() * (sampleUpperBound - sampleLowerBound);
		}
		Arrays.sort(randomSampleTimes);

		return randomSampleTimes;
	}
	
	protected void processEventsFromSpikingTimeSeries(double[] sourceSpikeTimes, double[] destSpikeTimes, double[][] conditionalSpikeTimes,
							  Vector<double[]> conditioningEmbeddingsFromSpikes, Vector<double[]> jointEmbeddingsFromSpikes,
							  Vector<double[]> conditioningEmbeddingsFromSamples, Vector<double[]> jointEmbeddingsFromSamples,
							  double actualNumSamplesMultiplier)
			throws Exception {

		int firstTargetIndexOfEMbedding = getFirstDestIndex(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes, true);
		double[] randomSampleTimes = generateRandomSampleTimes(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes, actualNumSamplesMultiplier, firstTargetIndexOfEMbedding);

		makeEmbeddingsAtPoints(destSpikeTimes, firstTargetIndexOfEMbedding, sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
				       conditioningEmbeddingsFromSpikes, jointEmbeddingsFromSpikes);
		makeEmbeddingsAtPoints(randomSampleTimes, 0, sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
				       conditioningEmbeddingsFromSamples, jointEmbeddingsFromSamples);
	}

	/*
	 * Method to do the 
	 */
	protected void processEventsFromSpikingTimeSeries(double[] sourceSpikeTimes, double[] destSpikeTimes, double[][] conditionalSpikeTimes,
							  Vector<double[]> conditioningEmbeddingsFromSamples, Vector<double[]> jointEmbeddingsFromSamples,
							  double actualNumSamplesMultiplier)
			throws Exception {

		int firstTargetIndexOfEMbedding = getFirstDestIndex(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes, false);
		double[] randomSampleTimes = generateRandomSampleTimes(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes, actualNumSamplesMultiplier, firstTargetIndexOfEMbedding);

		makeEmbeddingsAtPoints(randomSampleTimes, 0, sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
				       conditioningEmbeddingsFromSamples, jointEmbeddingsFromSamples);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#
	 * getAddedMoreThanOneObservationSet()
	 */
	@Override
	public boolean getAddedMoreThanOneObservationSet() {
		return (vectorOfDestinationSpikeTimes != null) && (vectorOfDestinationSpikeTimes.size() > 1);
	}

	// Class to allow returning two values in the subsequent method
	private static class distanceAndNumPoints {
		public double distance;
		public int numPoints;

		public distanceAndNumPoints(double distance, int numPoints) {
			this.distance = distance;
			this.numPoints = numPoints;
		}
	}
	
	private distanceAndNumPoints findMaxDistanceAndNumPointsFromIndices(double[] point, int[] indices, Vector<double[]> setOfPoints) {
		double maxDistance = 0;
		int i = 0;
		for (; indices[i] != -1; i++) {
			double distance = KdTree.norm(point, setOfPoints.elementAt(indices[i]), normType);
			if (distance > maxDistance) {
				maxDistance = distance;
			}
		}
		return new distanceAndNumPoints(maxDistance, i);
	}

	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		return computeAverageLocalOfObservations(kdTreeJointAtSpikes, jointEmbeddingsFromSpikes);
	}
	

	/*
	 * We take the actual joint tree at spikes (along with the associated embeddings) as an argument, as we will need to swap these out when
	 * computing surrogates.
	 */
	public double computeAverageLocalOfObservations(KdTree actualKdTreeJointAtSpikes, Vector<double[]> actualJointEmbeddingsFromSpikes) throws Exception {

		double currentSum = 0;
		for (int i = 0; i < conditioningEmbeddingsFromSpikes.size(); i++) {

			double radiusJointSpikes = actualKdTreeJointAtSpikes.findKNearestNeighbours(Knns, i).poll().norms[0];
			double radiusJointSamples = kdTreeJointAtSamples.findKNearestNeighbours(Knns,
					new double[][] { actualJointEmbeddingsFromSpikes.elementAt(i) }).poll().norms[0];

			/*
			  The algorithm specified in box 1 of doi.org/10.1371/journal.pcbi.1008054 specifies finding the maximum of the two radii 
			  just calculated and then redoing the searches in both sets at this radius. In this implementation, however, we make use 
			  of the fact that one radius is equal to the maximum, and so only one search needs to be redone. 
			 */
			double eps = 0.0;
			// Need variables for the number of neighbours as this is now variable within the maximum radius
			int kJointSpikes = 0;
			int kJointSamples = 0;
			if (radiusJointSpikes >= radiusJointSamples) {
				/* 
				   The maximum was the radius in the set of embeddings at spikes, so redo search in the set of embeddings at randomly
				   sampled points, using this larger radius.
				*/
				kJointSpikes = Knns;
				int[] indicesWithinR = new int[jointEmbeddingsFromSamples.size()];
				boolean[] isWithinR = new boolean[jointEmbeddingsFromSamples.size()];
				kdTreeJointAtSamples.findPointsWithinR(radiusJointSpikes + eps,
								       new double[][] { actualJointEmbeddingsFromSpikes.elementAt(i) },
								       true,
								       isWithinR,
								       indicesWithinR);
				distanceAndNumPoints temp = findMaxDistanceAndNumPointsFromIndices(actualJointEmbeddingsFromSpikes.elementAt(i), indicesWithinR,
												   jointEmbeddingsFromSamples);
				kJointSamples = temp.numPoints;
				radiusJointSamples = temp.distance;
			} else {
				/* 
				   The maximum was the radius in the set of embeddings at randomly sampled points, so redo search in the set of embeddings 
				   at spikes, using this larger radius.
				*/
				kJointSamples = Knns;
				int[] indicesWithinR = new int[jointEmbeddingsFromSamples.size()];
				boolean[] isWithinR = new boolean[jointEmbeddingsFromSamples.size()];
				actualKdTreeJointAtSpikes.findPointsWithinR(radiusJointSamples + eps,
								       new double[][] { actualJointEmbeddingsFromSpikes.elementAt(i) },
								       true,
								       isWithinR,
								       indicesWithinR);
				distanceAndNumPoints temp = findMaxDistanceAndNumPointsFromIndices(actualJointEmbeddingsFromSpikes.elementAt(i), indicesWithinR,
												   actualJointEmbeddingsFromSpikes);
				// -1 due to the point itself being in the set
				kJointSpikes = temp.numPoints - 1;
				radiusJointSpikes = temp.distance;
			}

			// Repeat the above steps, but in the conditioning (rather than joint) space.
			double radiusConditioningSpikes = kdTreeConditioningAtSpikes.findKNearestNeighbours(Knns, i).poll().norms[0];
			double radiusConditioningSamples = kdTreeConditioningAtSamples.findKNearestNeighbours(Knns, 
					new double[][] { conditioningEmbeddingsFromSpikes.elementAt(i) }).poll().norms[0];
			int kConditioningSpikes = 0;
			int kConditioningSamples = 0;
			if (radiusConditioningSpikes >= radiusConditioningSamples) {
				kConditioningSpikes = Knns;
				int[] indicesWithinR = new int[conditioningEmbeddingsFromSamples.size()];
				boolean[] isWithinR = new boolean[conditioningEmbeddingsFromSamples.size()];
				kdTreeConditioningAtSamples.findPointsWithinR(radiusConditioningSpikes + eps,
								       new double[][] { conditioningEmbeddingsFromSpikes.elementAt(i) },
								       true,
								       isWithinR,
								       indicesWithinR);
				distanceAndNumPoints temp = findMaxDistanceAndNumPointsFromIndices(conditioningEmbeddingsFromSpikes.elementAt(i), indicesWithinR,
												   conditioningEmbeddingsFromSamples);
				kConditioningSamples = temp.numPoints;
				radiusConditioningSamples = temp.distance;
			} else {
				kConditioningSamples = Knns;
				int[] indicesWithinR = new int[conditioningEmbeddingsFromSamples.size()];
				boolean[] isWithinR = new boolean[conditioningEmbeddingsFromSamples.size()];
				kdTreeConditioningAtSpikes.findPointsWithinR(radiusConditioningSamples + eps,
								       new double[][] { conditioningEmbeddingsFromSpikes.elementAt(i) },
								       true,
								       isWithinR,
								       indicesWithinR);
				distanceAndNumPoints temp = findMaxDistanceAndNumPointsFromIndices(conditioningEmbeddingsFromSpikes.elementAt(i), indicesWithinR,
												   conditioningEmbeddingsFromSpikes);
				// -1 due to the point itself being in the set
				kConditioningSpikes = temp.numPoints - 1;
				radiusConditioningSpikes = temp.distance;
			}
			
			/*
			 * TODO
			 * The KdTree class defaults to the squared euclidean distance when the euclidean norm is specified. This is fine for Kraskov estimators
			 * (as the radii are never used, just the numbers of points within radii). It causes problems here though, as we do use the radii and the
			 * squared euclidean distance is not a distance metric. We get around this by just taking the square root here, but it might be better to 
			 * fix this in the KdTree class. 
			 */
			if (normType == EuclideanUtils.NORM_EUCLIDEAN) {
				radiusJointSpikes = Math.sqrt(radiusJointSpikes);
				radiusJointSamples = Math.sqrt(radiusJointSamples);
				radiusConditioningSpikes = Math.sqrt(radiusConditioningSpikes);
				radiusConditioningSamples = Math.sqrt(radiusConditioningSamples);
			}
			
			currentSum += (MathsUtils.digamma(kJointSpikes) - MathsUtils.digamma(kJointSamples) + 
				       ((k + l) * (-Math.log(radiusJointSpikes) + Math.log(radiusJointSamples))) -
				       MathsUtils.digamma(kConditioningSpikes) + MathsUtils.digamma(kConditioningSamples) + 
				       + (k * (Math.log(radiusConditioningSpikes) - Math.log(radiusConditioningSamples))));
			if (Double.isNaN(currentSum)) {
				throw new Exception("NaNs in TE clac");
			}
		}
		// Normalise by time
		double timeSum = 0;
		for (Double time : processTimeLengths) {
			timeSum += time;
		}
		currentSum /= timeSum;
		return currentSum;
	}


	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck, double estimatedValue) throws Exception {
		return computeSignificance(numPermutationsToCheck,
					   estimatedValue,
					   System.currentTimeMillis());
	}

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck,
								    double estimatedValue, long randomSeed) throws Exception{

		Random random = new Random(randomSeed);
		double[] surrogateTEValues = new double[numPermutationsToCheck];

		for (int permutationNumber = 0; permutationNumber < numPermutationsToCheck; permutationNumber++) {
			Vector<double[]> resampledConditioningEmbeddingsFromSamples = new Vector<double[]>();
			Vector<double[]> resampledJointEmbeddingsFromSamples = new Vector<double[]>();

			// Send all of the observations through:
			Iterator<double[]> sourceIterator = vectorOfSourceSpikeTimes.iterator();
			Iterator<double[][]> conditionalIterator = null;
			if (vectorOfConditionalSpikeTimes.size() > 0) {
				conditionalIterator = vectorOfConditionalSpikeTimes.iterator();
			}
			int timeSeriesIndex = 0;
			for (double[] destSpikeTimes : vectorOfDestinationSpikeTimes) {
				double[] sourceSpikeTimes = sourceIterator.next();
				double[][] conditionalSpikeTimes = null;
				if (vectorOfConditionalSpikeTimes.size() > 0) {
					conditionalSpikeTimes = conditionalIterator.next();
				} else {
					conditionalSpikeTimes = new double[][] {};
				}
				processEventsFromSpikingTimeSeries(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
								   resampledConditioningEmbeddingsFromSamples, resampledJointEmbeddingsFromSamples,
								   surrogateNumSamplesMultiplier);
			}

			// Convert the vectors to arrays so that they can be put in the trees
			double[][] arrayedResampledConditioningEmbeddingsFromSamples = new double[resampledConditioningEmbeddingsFromSamples.size()][k];
			for (int i = 0; i < resampledConditioningEmbeddingsFromSamples.size(); i++) {
				arrayedResampledConditioningEmbeddingsFromSamples[i] = resampledConditioningEmbeddingsFromSamples.elementAt(i);
			}

			KdTree resampledKdTreeConditioningAtSamples = new KdTree(arrayedResampledConditioningEmbeddingsFromSamples);
			resampledKdTreeConditioningAtSamples.setNormType(normType);

			Vector<double[]> conditionallyPermutedJointEmbeddingsFromSpikes = new Vector(jointEmbeddingsFromSpikes);

			Vector<Integer> usedIndices = new Vector<Integer>();
			for (int i = 0; i < conditionallyPermutedJointEmbeddingsFromSpikes.size(); i++) {
				PriorityQueue<NeighbourNodeData> neighbours =
					resampledKdTreeConditioningAtSamples.findKNearestNeighbours(kPerm,
												    new double[][] {conditioningEmbeddingsFromSpikes.elementAt(i)});
				ArrayList<Integer> foundIndices = new ArrayList<Integer>();
				for (int j = 0; j < kPerm; j++) {
					foundIndices.add(neighbours.poll().sampleIndex);
				}

				ArrayList<Integer> prunedIndices = new ArrayList<Integer>(foundIndices);
				prunedIndices.removeAll(usedIndices);
				int chosenIndex = 0;
				if (prunedIndices.size() > 0) {
					chosenIndex = prunedIndices.get(random.nextInt(prunedIndices.size()));
				} else {
					chosenIndex = foundIndices.get(random.nextInt(foundIndices.size()));
				}
				usedIndices.add(chosenIndex);
				int embeddingLength = conditionallyPermutedJointEmbeddingsFromSpikes.elementAt(i).length;
				for(int j = 0; j < l; j++) {
					conditionallyPermutedJointEmbeddingsFromSpikes.elementAt(i)[embeddingLength - l + j] =
						resampledJointEmbeddingsFromSamples.elementAt(chosenIndex)[embeddingLength - l + j];
				}
			
			}

			double[][] arrayedConditionallyPermutedJointEmbeddingsFromSpikes = new double[conditionallyPermutedJointEmbeddingsFromSpikes.size()][];
			for (int i = 0; i < conditionallyPermutedJointEmbeddingsFromSpikes.size(); i++) {
				arrayedConditionallyPermutedJointEmbeddingsFromSpikes[i] = conditionallyPermutedJointEmbeddingsFromSpikes.elementAt(i);
			}
			KdTree conditionallyPermutedKdTreeJointFromSpikes = new KdTree(arrayedConditionallyPermutedJointEmbeddingsFromSpikes);
			conditionallyPermutedKdTreeJointFromSpikes.setNormType(normType);
		
			surrogateTEValues[permutationNumber] = computeAverageLocalOfObservations(conditionallyPermutedKdTreeJointFromSpikes,
												 conditionallyPermutedJointEmbeddingsFromSpikes);
		}
		return new EmpiricalMeasurementDistribution(surrogateTEValues, estimatedValue);
	}
	

	/*
	 * (non-Javadoc)
	 * 
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#
	 * computeLocalOfPreviousObservations()
	 */
	@Override
	public SpikingLocalInformationValues computeLocalOfPreviousObservations() throws Exception {
		// TODO Auto-generated method stub
		return null;
	}


	/*
	 * (non-Javadoc)
	 * 
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#setDebug(
	 * boolean)
	 */
	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#getLastAverage
	 * ()
	 */
	@Override
	public double getLastAverage() {
		// TODO Auto-generated method stub
		return 0;
	}
}
