package infodynamics.measures.spiking.integration;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Vector;
import java.util.ArrayList;
import java.lang.Math;

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
 * KdTree class (or may be do-able with existing KdTree class if excluding based on nearest target spikes
 * rather than in terms of timing)
 */
public class TransferEntropyCalculatorSpikingIntegration implements TransferEntropyCalculatorSpiking {

	/**
	 * The past destination interspike intervals to consider 
	 * (and associated property name and convenience length variable)
	 * The code assumes that the first interval is numbered 1, the next is numbered 2, etc.
	 * It also assumes that the intervals are sorted. The setter method performs sorting to ensure this.
	 */
	public static final String DEST_PAST_INTERVALS_PROP_NAME = "DEST_PAST_INTERVALS";
	protected int[] destPastIntervals = new int[] {};
	protected int numDestPastIntervals = 0;
	/**
	 * As above but for the sources
	 */
	public static final String SOURCE_PAST_INTERVALS_PROP_NAME = "SOURCE_PAST_INTERVALS";
	protected int[] sourcePastIntervals = new int[] {};
	protected int numSourcePastIntervals = 0;

	/**
	 * As above but for the conditioning processes. There is no property name for this variable,
	 * due to there not currently being a method for converting strings to 2d arrays. Instead,
	 * a separate setter method is implemented.
	 */
	protected Vector<int[]> vectorOfCondPastIntervals = new Vector<int[]>();
	protected int numCondPastIntervals = 0;

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
	 * Whether to use the jittered sampling approach. Useful for bursty spike trains. Explained in methods section of
	 * doi.org/10.1101/2021.06.29.450432
	 */
        protected boolean jitteredSamplesForSurrogates = false;
        public static final String DO_JITTERED_SAMPLING_PROP_NAME = "DO_JITTERED_SAMPLING";
        protected double jitteredSamplingNoiseLevel = 1;
        public static final String JITTERED_SAMPLING_NOISE_LEVEL = "JITTERED_SAMPLING_NOISE_LEVEL";

	/**
	 * Stores whether we are in debug mode
	 */
	protected boolean debug = false;

	/**
	 * Property name for the number of random sample points to use as a multiple
	 * of the number of target spikes.
	 */
	public static final String PROP_SAMPLE_MULTIPLIER = "NUM_SAMPLES_MULTIPLIER";
	protected double numSamplesMultiplier = 2.0;
	/**
	 * Property name for the number of random sample points to use in the construction of the surrogates as a multiple
	 * of the number of target spikes.
	 */
	public static final String PROP_SURROGATE_SAMPLE_MULTIPLIER = "SURROGATE_NUM_SAMPLES_MULTIPLIER";
	protected double surrogateNumSamplesMultiplier = 2.0;

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
		initialise(0, 0);
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
		initialise(0, 0);
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
		if (propertyName.equalsIgnoreCase(DEST_PAST_INTERVALS_PROP_NAME)) {
			if (propertyValue.length() == 0) {
				destPastIntervals = new int[] {};
			} else {
				int[] destPastIntervalsTemp = ParsedProperties.parseStringArrayOfInts(propertyValue);
				for (int interval : destPastIntervalsTemp) {
					if (interval < 1) {
						throw new Exception ("Invalid interval number less than 1.");
					}
				}
				destPastIntervals = destPastIntervalsTemp;
				Arrays.sort(destPastIntervals);
			}
		} else if (propertyName.equalsIgnoreCase(SOURCE_PAST_INTERVALS_PROP_NAME)) {
			if (propertyValue.length() == 0) {
				sourcePastIntervals = new int[] {};
			} else {
				int[] sourcePastIntervalsTemp = ParsedProperties.parseStringArrayOfInts(propertyValue);
				for (int interval : sourcePastIntervalsTemp) {
					if (interval < 1) {
						throw new Exception ("Invalid interval number less than 1.");
					}
				}
				sourcePastIntervals = sourcePastIntervalsTemp;
				Arrays.sort(sourcePastIntervals);
			}
		} else if (propertyName.equalsIgnoreCase(KNNS_PROP_NAME)) {
			Knns = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(DO_JITTERED_SAMPLING_PROP_NAME)) {
		        jitteredSamplesForSurrogates = Boolean.parseBoolean(propertyValue);
		} else if (propertyName.equalsIgnoreCase(JITTERED_SAMPLING_NOISE_LEVEL)) {
		        jitteredSamplingNoiseLevel = Double.parseDouble(propertyValue);
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
		if (propertyName.equalsIgnoreCase(KNNS_PROP_NAME)) {
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


	public void appendConditionalIntervals(int[] intervals) throws Exception{
		for (int interval : intervals) {
			if (interval < 1) {
				throw new Exception ("Invalid interval number less than 1.");
			}
		}
		Arrays.sort(intervals);
		vectorOfCondPastIntervals.add(intervals);
	}

	public void clearConditionalIntervals() {
		vectorOfCondPastIntervals = new Vector<int[]>();
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

		// Set these conveniance variables as they are used quite a bit later on
		numDestPastIntervals = destPastIntervals.length;
		numSourcePastIntervals = sourcePastIntervals.length;
		numCondPastIntervals = 0;
		for (int[] intervals : vectorOfCondPastIntervals) {
			numCondPastIntervals += intervals.length;
		}

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
							   numSamplesMultiplier, false);
		}



		// Convert the vectors to arrays so that they can be put in the trees
		double[][] arrayedTargetEmbeddingsFromSpikes = new double[conditioningEmbeddingsFromSpikes.size()][numDestPastIntervals + numCondPastIntervals];
		double[][] arrayedJointEmbeddingsFromSpikes = new double[conditioningEmbeddingsFromSpikes.size()][numDestPastIntervals +
														  numCondPastIntervals + numSourcePastIntervals];
		for (int i = 0; i < conditioningEmbeddingsFromSpikes.size(); i++) {
			arrayedTargetEmbeddingsFromSpikes[i] = conditioningEmbeddingsFromSpikes.elementAt(i);
			arrayedJointEmbeddingsFromSpikes[i] = jointEmbeddingsFromSpikes.elementAt(i);
		}
		double[][] arrayedTargetEmbeddingsFromSamples = new double[conditioningEmbeddingsFromSamples.size()][numDestPastIntervals + numCondPastIntervals];
		double[][] arrayedJointEmbeddingsFromSamples = new double[conditioningEmbeddingsFromSamples.size()][numDestPastIntervals +
														  numCondPastIntervals + numSourcePastIntervals];
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

		// Initialise the starting points of all the tracking variables
		int embeddingPointIndex = indexOfFirstPointToUse;
		int mostRecentDestIndex = destPastIntervals[destPastIntervals.length - 1];
		int mostRecentSourceIndex = sourcePastIntervals[sourcePastIntervals.length - 1] - 1;
		int[] mostRecentConditioningIndices = new int[vectorOfCondPastIntervals.size()];
		for (int i = 0; i < vectorOfCondPastIntervals.size(); i++) {
			mostRecentConditioningIndices[i] = vectorOfCondPastIntervals.elementAt(i)[vectorOfCondPastIntervals.elementAt(i).length - 1] - 1;
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
				if (sourceSpikeTimes[mostRecentSourceIndex + 1] < pointsAtWhichToMakeEmbeddings[embeddingPointIndex]) {
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

			
			double[] conditioningPast = new double[numDestPastIntervals + numCondPastIntervals];
			double[] jointPast = new double[numDestPastIntervals + numCondPastIntervals + numSourcePastIntervals];

			// Add the embedding intervals from the target process
			for (int i = 0; i < destPastIntervals.length; i++) {
				// Case where we are inserting an interval from an observation point back to the most recent event in the target process
				if (destPastIntervals[i] == 1) {
					conditioningPast[i] = pointsAtWhichToMakeEmbeddings[embeddingPointIndex] - destSpikeTimes[mostRecentDestIndex];
					jointPast[i] = pointsAtWhichToMakeEmbeddings[embeddingPointIndex]
						- destSpikeTimes[mostRecentDestIndex];
				// Case where we are inserting an inter-event intervvl from the target process
				} else {
					conditioningPast[i] = destSpikeTimes[mostRecentDestIndex - destPastIntervals[i] + 2]
						- destSpikeTimes[mostRecentDestIndex - destPastIntervals[i] + 1];
					jointPast[i] = destSpikeTimes[mostRecentDestIndex - destPastIntervals[i] + 2]
						- destSpikeTimes[mostRecentDestIndex - destPastIntervals[i] + 1];
				}
			}

			// Add the embeding intervals from the conditional processes
			int indexOfNextEmbeddingInterval = numDestPastIntervals;
			for (int i = 0; i < vectorOfCondPastIntervals.size(); i++) {
				for (int j = 0; j < vectorOfCondPastIntervals.elementAt(i).length; j++) {
					// Case where we are inserting an interval from an observation point back to the most recent event in the conditioning process
					if (vectorOfCondPastIntervals.elementAt(i)[j] == 1) {
						conditioningPast[indexOfNextEmbeddingInterval] = pointsAtWhichToMakeEmbeddings[embeddingPointIndex]
							- conditionalSpikeTimes[i][mostRecentConditioningIndices[i]];
						jointPast[indexOfNextEmbeddingInterval] = pointsAtWhichToMakeEmbeddings[embeddingPointIndex]
							- conditionalSpikeTimes[i][mostRecentConditioningIndices[i]];
					// Case where we are inserting an inter-event interval from the conditioning process
					} else {
						// Convenience variable
						int intervalNumber = vectorOfCondPastIntervals.elementAt(i)[j];
						conditioningPast[indexOfNextEmbeddingInterval] = conditionalSpikeTimes[i][mostRecentConditioningIndices[i] - intervalNumber + 2]
							- conditionalSpikeTimes[i][mostRecentConditioningIndices[i] - intervalNumber + 1];
						jointPast[indexOfNextEmbeddingInterval] = conditionalSpikeTimes[i][mostRecentConditioningIndices[i] - intervalNumber + 2]
							- conditionalSpikeTimes[i][mostRecentConditioningIndices[i] - intervalNumber + 1];
					}
					indexOfNextEmbeddingInterval++;
				}
			}

			// Add the embedding intervals from the source process (this only gets added to the joint embeddings)
			for (int i = 0; i < sourcePastIntervals.length; i++) {
				// Case where we are inserting an interval from an observation point back to the most recent event in the source process
				if (sourcePastIntervals[i] == 1) {
					jointPast[indexOfNextEmbeddingInterval] = pointsAtWhichToMakeEmbeddings[embeddingPointIndex]
						- sourceSpikeTimes[mostRecentSourceIndex];
			        // Case where we are inserting an inter-event interval from the source process					
				} else {
					jointPast[indexOfNextEmbeddingInterval] = sourceSpikeTimes[mostRecentSourceIndex - sourcePastIntervals[i] + 2]
						- sourceSpikeTimes[mostRecentSourceIndex - sourcePastIntervals[i] + 1];
				}
				indexOfNextEmbeddingInterval++;
			}

			// Add Gaussian noise, if necessary
			if (addNoise) {
				for (int i = 0; i < conditioningPast.length; i++) {
					conditioningPast[i] = Math.log(conditioningPast[i] + 1.1);
					conditioningPast[i] += random.nextGaussian() * noiseLevel;
				}
				for (int i = 0; i < jointPast.length; i++) {
					if (jointPast[i] < 0) {
						System.out.println("NEGATIVE");
					}
					jointPast[i] = Math.log(jointPast[i] + 1.1);
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
		
		int firstTargetIndexOfEmbedding = destPastIntervals[destPastIntervals.length - 1];

		int furthestInterval = sourcePastIntervals[sourcePastIntervals.length - 1];
		while (destSpikeTimes[firstTargetIndexOfEmbedding] < sourceSpikeTimes[furthestInterval - 1]) {
			firstTargetIndexOfEmbedding++;
		}
		if (conditionalSpikeTimes.length != vectorOfCondPastIntervals.size()) {
			throw new Exception("Number of conditional embedding lengths does not match the number of conditional processes");
		}		
		for (int i = 0; i < conditionalSpikeTimes.length; i++) {
			furthestInterval = vectorOfCondPastIntervals.elementAt(i)[vectorOfCondPastIntervals.elementAt(i).length - 1];
			while (destSpikeTimes[firstTargetIndexOfEmbedding] <= conditionalSpikeTimes[i][furthestInterval - 1]) {
				firstTargetIndexOfEmbedding++;
			}
		}

		// We don't want to reset these lengths when resampling for surrogates
		if (setProcessTimeLengths) {
			processTimeLengths.add(destSpikeTimes[destSpikeTimes.length - 1] - destSpikeTimes[firstTargetIndexOfEmbedding]);
		}
		
		return firstTargetIndexOfEmbedding;
	}

	protected double[] generateRandomSampleTimes(double[] sourceSpikeTimes, double[] destSpikeTimes, double[][] conditionalSpikeTimes,
						     double actualNumSamplesMultiplier, int firstTargetIndexOfEmbedding, boolean doJitteredSampling) {

		double sampleLowerBound = destSpikeTimes[firstTargetIndexOfEmbedding];
		double sampleUpperBound = destSpikeTimes[destSpikeTimes.length - 1];
		int num_samples = (int) Math.round(actualNumSamplesMultiplier * (destSpikeTimes.length - firstTargetIndexOfEmbedding + 1));
		double[] randomSampleTimes = new double[num_samples];
		Random rand = new Random();
		if (doJitteredSampling) {
		    //System.out.println("jittering " + jitteredSamplingNoiseLevel);
			for (int i = 0; i < randomSampleTimes.length; i++) {
				randomSampleTimes[i] = destSpikeTimes[firstTargetIndexOfEmbedding + (i % (destSpikeTimes.length - firstTargetIndexOfEmbedding - 1))]
					+ jitteredSamplingNoiseLevel * (rand.nextDouble() - 0.5);
				//randomSampleTimes[i] = -1.0;
				if ((randomSampleTimes[i] > sampleUpperBound) || (randomSampleTimes[i] < sampleLowerBound)) {
					randomSampleTimes[i] = sampleLowerBound + rand.nextDouble() * (sampleUpperBound - sampleLowerBound);
				}
			}
		} else {
			for (int i = 0; i < randomSampleTimes.length; i++) {
				randomSampleTimes[i] = sampleLowerBound + rand.nextDouble() * (sampleUpperBound - sampleLowerBound);
			}
		}
		Arrays.sort(randomSampleTimes);
		/*System.out.println(sampleLowerBound + " " + sampleUpperBound);
		for (int i = 0; i < randomSampleTimes.length; i += 100) {
			System.out.print(randomSampleTimes[i] + " ");			
		}
		System.out.println("\n\n\n");*/
		return randomSampleTimes;
	}
	
	protected void processEventsFromSpikingTimeSeries(double[] sourceSpikeTimes, double[] destSpikeTimes, double[][] conditionalSpikeTimes,
							  Vector<double[]> conditioningEmbeddingsFromSpikes, Vector<double[]> jointEmbeddingsFromSpikes,
							  Vector<double[]> conditioningEmbeddingsFromSamples, Vector<double[]> jointEmbeddingsFromSamples,
							  double actualNumSamplesMultiplier, boolean doJitteredSampling)
			throws Exception {

		int firstTargetIndexOfEmbedding = getFirstDestIndex(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes, true);
		double[] randomSampleTimes = generateRandomSampleTimes(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
								       actualNumSamplesMultiplier, firstTargetIndexOfEmbedding,
								       doJitteredSampling);

		makeEmbeddingsAtPoints(destSpikeTimes, firstTargetIndexOfEmbedding, sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
				       conditioningEmbeddingsFromSpikes, jointEmbeddingsFromSpikes);
		makeEmbeddingsAtPoints(randomSampleTimes, 0, sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
				       conditioningEmbeddingsFromSamples, jointEmbeddingsFromSamples);
	}

	protected void processEventsFromSpikingTimeSeries(double[] sourceSpikeTimes, double[] destSpikeTimes, double[][] conditionalSpikeTimes,
							  Vector<double[]> conditioningEmbeddingsFromSamples, Vector<double[]> jointEmbeddingsFromSamples,
							  double actualNumSamplesMultiplier, boolean doJitteredSampling)
			throws Exception {

		int firstTargetIndexOfEmbedding = getFirstDestIndex(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes, false);
		double[] randomSampleTimes = generateRandomSampleTimes(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
								       actualNumSamplesMultiplier, firstTargetIndexOfEmbedding,
								       doJitteredSampling);
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
			double tempRadiusJointSamples = radiusJointSamples;
			if (normType == EuclideanUtils.NORM_EUCLIDEAN) {
				radiusJointSpikes = Math.sqrt(radiusJointSpikes);
				radiusJointSamples = Math.sqrt(radiusJointSamples);
				radiusConditioningSpikes = Math.sqrt(radiusConditioningSpikes);
				radiusConditioningSamples = Math.sqrt(radiusConditioningSamples);
			}
			
			currentSum += (MathsUtils.digamma(kJointSpikes) - MathsUtils.digamma(kJointSamples) + 
				       ((numDestPastIntervals + numCondPastIntervals + numSourcePastIntervals) * (-Math.log(radiusJointSpikes) + Math.log(radiusJointSamples))) -
				       MathsUtils.digamma(kConditioningSpikes) + MathsUtils.digamma(kConditioningSamples) + 
				       + ((numDestPastIntervals + numCondPastIntervals) * (Math.log(radiusConditioningSpikes) - Math.log(radiusConditioningSamples))));
			if (Double.isNaN(currentSum)) {
				for (double[] embed : jointEmbeddingsFromSamples) {
					System.out.println(Arrays.toString(embed));
				}
				throw new Exception("NaNs in TE clac " + radiusJointSpikes + " " + radiusJointSamples + " " + tempRadiusJointSamples);
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
				if (jitteredSamplesForSurrogates) {
				    processEventsFromSpikingTimeSeries(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
								       resampledConditioningEmbeddingsFromSamples, resampledJointEmbeddingsFromSamples,
								       surrogateNumSamplesMultiplier, true);
				} else {
				    processEventsFromSpikingTimeSeries(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
								       resampledConditioningEmbeddingsFromSamples, resampledJointEmbeddingsFromSamples,
								       surrogateNumSamplesMultiplier, false);
				}
			}
			// Convert the vectors to arrays so that they can be put in the trees
			double[][] arrayedResampledConditioningEmbeddingsFromSamples = new double[resampledConditioningEmbeddingsFromSamples.size()][numDestPastIntervals];
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
				for(int j = 0; j < numSourcePastIntervals; j++) {
					conditionallyPermutedJointEmbeddingsFromSpikes.elementAt(i)[embeddingLength - numSourcePastIntervals + j] =
						resampledJointEmbeddingsFromSamples.elementAt(chosenIndex)[embeddingLength - numSourcePastIntervals + j];
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
