package infodynamics.measures.spiking.integration;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Vector;

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

	Vector<double[]> ConditioningEmbeddingsFromSpikes = null;
	Vector<double[]> jointEmbeddingsFromSpikes = null;
	Vector<double[]> ConditioningEmbeddingsFromSamples = null;
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
	protected double num_samples_multiplier = 1.0;
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
			int k_temp = Integer.parseInt(propertyValue);
			if (k_temp < 1) {
				throw new Exception ("Invalid k value less than 1.");
			} else {
				k = k_temp;
			}
		} else if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			int l_temp = Integer.parseInt(propertyValue);
			if (l_temp < 1) {
				throw new Exception ("Invalid l value less than 1.");
			} else {
				l = l_temp;
			}
		} else if (propertyName.equalsIgnoreCase(COND_EMBED_LENGTHS_PROP_NAME)) {
			int[] condEmbedDims_temp = ParsedProperties.parseStringArrayOfInts(propertyValue);
			for (int dim : condEmbedDims_temp) {
				if (dim < 1) {
					throw new Exception ("Invalid conditional embedding value less than 1.");
				}
			}
			condEmbedDims = condEmbedDims_temp;
		} else if (propertyName.equalsIgnoreCase(KNNS_PROP_NAME)) {
			Knns = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
			if (propertyValue.equals("0") || propertyValue.equalsIgnoreCase("false")) {
				addNoise = false;
				noiseLevel = 0;
			} else {
				addNoise = true;
				noiseLevel = Double.parseDouble(propertyValue);
			}
			
		} else if (propertyName.equalsIgnoreCase(PROP_SAMPLE_MULTIPLIER)) {
			double temp_num_samples_multiplier = Double.parseDouble(propertyValue);
			if (temp_num_samples_multiplier <= 0) {
				throw new Exception ("Num samples multiplier must be greater than 0.");
			} else {
				num_samples_multiplier = temp_num_samples_multiplier;
			}
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normType = KdTree.validateNormType(propertyValue);
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
			return Double.toString(num_samples_multiplier);
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

		ConditioningEmbeddingsFromSpikes = new Vector<double[]>();
		jointEmbeddingsFromSpikes = new Vector<double[]>();
		ConditioningEmbeddingsFromSamples = new Vector<double[]>();
		jointEmbeddingsFromSamples = new Vector<double[]>();
		processTimeLengths = new Vector<Double>();

		// Send all of the observations through:
		Iterator<double[]> sourceIterator = vectorOfSourceSpikeTimes.iterator();
		int timeSeriesIndex = 0;
		if (vectorOfConditionalSpikeTimes.size() > 0) {
			Iterator<double[][]> conditionalIterator = vectorOfConditionalSpikeTimes.iterator();
			for (double[] destSpikeTimes : vectorOfDestinationSpikeTimes) {
				double[] sourceSpikeTimes = sourceIterator.next();
				double[][] conditionalSpikeTimes = conditionalIterator.next();
				processEventsFromSpikingTimeSeries(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes, ConditioningEmbeddingsFromSpikes,
								   jointEmbeddingsFromSpikes, ConditioningEmbeddingsFromSamples, jointEmbeddingsFromSamples,
								   processTimeLengths);
			}
		} else {
			for (double[] destSpikeTimes : vectorOfDestinationSpikeTimes) {
				double[] sourceSpikeTimes = sourceIterator.next();
				double[][] conditionalSpikeTimes = new double[][] {};
				processEventsFromSpikingTimeSeries(sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes, ConditioningEmbeddingsFromSpikes,
								   jointEmbeddingsFromSpikes, ConditioningEmbeddingsFromSamples, jointEmbeddingsFromSamples,
								   processTimeLengths);
			}
		}

		// Convert the vectors to arrays so that they can be put in the trees
		double[][] arrayedTargetEmbeddingsFromSpikes = new double[ConditioningEmbeddingsFromSpikes.size()][k];
		double[][] arrayedJointEmbeddingsFromSpikes = new double[ConditioningEmbeddingsFromSpikes.size()][k + l];
		for (int i = 0; i < ConditioningEmbeddingsFromSpikes.size(); i++) {
			arrayedTargetEmbeddingsFromSpikes[i] = ConditioningEmbeddingsFromSpikes.elementAt(i);
			arrayedJointEmbeddingsFromSpikes[i] = jointEmbeddingsFromSpikes.elementAt(i);
		}
		double[][] arrayedTargetEmbeddingsFromSamples = new double[ConditioningEmbeddingsFromSamples.size()][k];
		double[][] arrayedJointEmbeddingsFromSamples = new double[ConditioningEmbeddingsFromSamples.size()][k + l];
		for (int i = 0; i < ConditioningEmbeddingsFromSamples.size(); i++) {
			arrayedTargetEmbeddingsFromSamples[i] = ConditioningEmbeddingsFromSamples.elementAt(i);
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

	protected void makeEmbeddingsAtPoints(double[] pointsAtWhichToMakeEmbeddings, int index_of_first_point_to_use,
					      double[] sourceSpikeTimes, double[] destSpikeTimes,
					      double[][] conditionalSpikeTimes,
					      Vector<double[]> ConditioningEmbeddings, 
					      Vector<double[]> jointEmbeddings) {

		Random random = new Random();

		int embedding_point_index = index_of_first_point_to_use;
		int most_recent_dest_index = k;
		int most_recent_source_index = l;
		int[] most_recent_conditioning_indices = Arrays.copyOf(condEmbedDims, condEmbedDims.length);
		int total_length_of_conditioning_embeddings = 0;
		for (int i = 0; i < condEmbedDims.length; i++) {
			total_length_of_conditioning_embeddings += condEmbedDims[i];
		}

		// Loop through the points at which embeddings need to be made
		for (; embedding_point_index < pointsAtWhichToMakeEmbeddings.length; embedding_point_index++) {

			// Advance the tracker of the most recent dest index
			while (most_recent_dest_index < (destSpikeTimes.length - 1)) {
				if (destSpikeTimes[most_recent_dest_index + 1] < pointsAtWhichToMakeEmbeddings[embedding_point_index]) {
					most_recent_dest_index++;
				} else {
					break;
				}
			}
			// Do the same for the most recent source index
			while (most_recent_source_index < (sourceSpikeTimes.length - 1)) {
				if (sourceSpikeTimes[most_recent_source_index
						+ 1] < pointsAtWhichToMakeEmbeddings[embedding_point_index]) {
					most_recent_source_index++;
				} else {
					break;
				}
			}
			// Now advance the trackers for the most recent conditioning indices
			for (int j = 0; j < most_recent_conditioning_indices.length; j++) {
				while (most_recent_conditioning_indices[j] < (conditionalSpikeTimes[j].length - 1)) {
					if (conditionalSpikeTimes[j][most_recent_conditioning_indices[j] + 1] < pointsAtWhichToMakeEmbeddings[embedding_point_index]) {
						most_recent_conditioning_indices[j]++;
					} else {
						break;
					}
				}
			}

			
			double[] conditioningPast = new double[k + total_length_of_conditioning_embeddings];
			double[] jointPast = new double[k + total_length_of_conditioning_embeddings + l];

			// Add the embedding intervals from the target process
			conditioningPast[0] = pointsAtWhichToMakeEmbeddings[embedding_point_index] - destSpikeTimes[most_recent_dest_index];
			jointPast[0] = pointsAtWhichToMakeEmbeddings[embedding_point_index]
					- destSpikeTimes[most_recent_dest_index];
			for (int i = 1; i < k; i++) {
				conditioningPast[i] = destSpikeTimes[most_recent_dest_index - i + 1]
						- destSpikeTimes[most_recent_dest_index - i];
				jointPast[i] = destSpikeTimes[most_recent_dest_index - i + 1]
						- destSpikeTimes[most_recent_dest_index - i];
			}

			// Add the embeding intervals from the conditional processes
			int index_of_next_embedding_interval = k;
			for (int i = 0; i < condEmbedDims.length; i++) {
				conditioningPast[index_of_next_embedding_interval] =
					pointsAtWhichToMakeEmbeddings[embedding_point_index] - conditionalSpikeTimes[i][most_recent_conditioning_indices[i]]; 
				jointPast[index_of_next_embedding_interval] = 
					pointsAtWhichToMakeEmbeddings[embedding_point_index] - conditionalSpikeTimes[i][most_recent_conditioning_indices[i]];
				index_of_next_embedding_interval += 1;
				for (int j = 1; j < condEmbedDims[i]; j++) {
					conditioningPast[index_of_next_embedding_interval] =
						conditionalSpikeTimes[i][most_recent_conditioning_indices[i] - j + 1] -
						conditionalSpikeTimes[i][most_recent_conditioning_indices[i] - j]; 
					jointPast[index_of_next_embedding_interval] =
						conditionalSpikeTimes[i][most_recent_conditioning_indices[i] - j + 1] -
						conditionalSpikeTimes[i][most_recent_conditioning_indices[i] - j]; 
					index_of_next_embedding_interval += 1;
				}
			}

			// Add the embedding intervals from the source process (this only gets added to the joint embeddings)
			jointPast[k + total_length_of_conditioning_embeddings] = pointsAtWhichToMakeEmbeddings[embedding_point_index]
					- sourceSpikeTimes[most_recent_source_index];
			for (int i = 1; i < l; i++) {
				jointPast[k + total_length_of_conditioning_embeddings + i] = sourceSpikeTimes[most_recent_source_index - i + 1]
						- sourceSpikeTimes[most_recent_source_index - i];
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

			ConditioningEmbeddings.add(conditioningPast);
			jointEmbeddings.add(jointPast);
		}
	}

	protected void processEventsFromSpikingTimeSeries(double[] sourceSpikeTimes, double[] destSpikeTimes, double[][] conditionalSpikeTimes,
			Vector<double[]> ConditioningEmbeddingsFromSpikes, Vector<double[]> jointEmbeddingsFromSpikes,
							  Vector<double[]> ConditioningEmbeddingsFromSamples, Vector<double[]> jointEmbeddingsFromSamples,
							  Vector<Double> processTimeLengths)
			throws Exception {

		// First sort the spike times in case they were not properly in ascending order:
		Arrays.sort(sourceSpikeTimes);
		Arrays.sort(destSpikeTimes);
		
		int first_target_index_of_embedding = k;
		while (destSpikeTimes[first_target_index_of_embedding] <= sourceSpikeTimes[l - 1]) {
			first_target_index_of_embedding++;
		}
		if (conditionalSpikeTimes.length != condEmbedDims.length) {
			throw new Exception("Number of conditional embedding lengths does not match the number of conditional processes");
		}
		for (int i = 0; i < conditionalSpikeTimes.length; i++) {
			while (destSpikeTimes[first_target_index_of_embedding] <= conditionalSpikeTimes[i][condEmbedDims[i]]) {
				first_target_index_of_embedding++;
			}
		}

		//processTimeLengths.add(destSpikeTimes[sourceSpikeTimes.length - 1] - destSpikeTimes[first_target_index_of_embedding]);
		processTimeLengths.add(destSpikeTimes[destSpikeTimes.length - 1] - destSpikeTimes[first_target_index_of_embedding]);
		
		double sample_lower_bound = destSpikeTimes[first_target_index_of_embedding];
		double sample_upper_bound = destSpikeTimes[destSpikeTimes.length - 1];
		int num_samples = (int) Math.round(num_samples_multiplier * (destSpikeTimes.length - first_target_index_of_embedding + 1));
		double[] randomSampleTimes = new double[num_samples];
		Random rand = new Random();
		for (int i = 0; i < randomSampleTimes.length; i++) {
			randomSampleTimes[i] = sample_lower_bound + rand.nextDouble() * (sample_upper_bound - sample_lower_bound);
		}
		Arrays.sort(randomSampleTimes);

		makeEmbeddingsAtPoints(destSpikeTimes, first_target_index_of_embedding, sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
				       ConditioningEmbeddingsFromSpikes, jointEmbeddingsFromSpikes);
		makeEmbeddingsAtPoints(randomSampleTimes, 0, sourceSpikeTimes, destSpikeTimes, conditionalSpikeTimes,
				       ConditioningEmbeddingsFromSamples, jointEmbeddingsFromSamples);
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#
	 * computeAverageLocalOfObservations()
	 */
	@Override
	public double computeAverageLocalOfObservations() throws Exception {

		double currentSum = 0;
		for (int i = 0; i < ConditioningEmbeddingsFromSpikes.size(); i++) {

			double radiusJointSpikes = kdTreeJointAtSpikes.findKNearestNeighbours(Knns, i).poll().norms[0];
			double radiusJointSamples = kdTreeJointAtSamples.findKNearestNeighbours(Knns,
					new double[][] { jointEmbeddingsFromSpikes.elementAt(i) }).poll().norms[0];

			/*
			  The algorithm specified in box 1 of doi.org/10.1371/journal.pcbi.1008054 specifies finding the maximum of the two radii 
			  just calculated and then redoing the searches in both sets at this radius. In this implementation, however, we make use 
			  of the fact that one radius is equal to the maximum, and so only one search needs to be redone. 
			 */
			double eps = 0.01;
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
								       new double[][] { jointEmbeddingsFromSpikes.elementAt(i) },
								       true,
								       isWithinR,
								       indicesWithinR);
				distanceAndNumPoints temp = findMaxDistanceAndNumPointsFromIndices(jointEmbeddingsFromSpikes.elementAt(i), indicesWithinR,
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
				kdTreeJointAtSpikes.findPointsWithinR(radiusJointSamples + eps,
								       new double[][] { jointEmbeddingsFromSpikes.elementAt(i) },
								       true,
								       isWithinR,
								       indicesWithinR);
				distanceAndNumPoints temp = findMaxDistanceAndNumPointsFromIndices(jointEmbeddingsFromSpikes.elementAt(i), indicesWithinR,
												   jointEmbeddingsFromSpikes);
				// -1 due to the point itself being in the set
				kJointSpikes = temp.numPoints - 1;
				radiusJointSpikes = temp.distance;
			}

			// Repeat the above steps, but in the conditioning (rather than joint) space.
			double radiusConditioningSpikes = kdTreeConditioningAtSpikes.findKNearestNeighbours(Knns, i).poll().norms[0];
			double radiusConditioningSamples = kdTreeConditioningAtSamples.findKNearestNeighbours(Knns, 
					new double[][] { ConditioningEmbeddingsFromSpikes.elementAt(i) }).poll().norms[0];
			int kConditioningSpikes = 0;
			int kConditioningSamples = 0;
			if (radiusConditioningSpikes >= radiusConditioningSamples) {
				kConditioningSpikes = Knns;
				int[] indicesWithinR = new int[ConditioningEmbeddingsFromSamples.size()];
				boolean[] isWithinR = new boolean[ConditioningEmbeddingsFromSamples.size()];
				kdTreeConditioningAtSamples.findPointsWithinR(radiusConditioningSpikes + eps,
								       new double[][] { ConditioningEmbeddingsFromSpikes.elementAt(i) },
								       true,
								       isWithinR,
								       indicesWithinR);
				distanceAndNumPoints temp = findMaxDistanceAndNumPointsFromIndices(ConditioningEmbeddingsFromSpikes.elementAt(i), indicesWithinR,
												   ConditioningEmbeddingsFromSamples);
				kConditioningSamples = temp.numPoints;
				radiusConditioningSamples = temp.distance;
			} else {
				kConditioningSamples = Knns;
				int[] indicesWithinR = new int[ConditioningEmbeddingsFromSamples.size()];
				boolean[] isWithinR = new boolean[ConditioningEmbeddingsFromSamples.size()];
				kdTreeConditioningAtSpikes.findPointsWithinR(radiusConditioningSamples + eps,
								       new double[][] { ConditioningEmbeddingsFromSpikes.elementAt(i) },
								       true,
								       isWithinR,
								       indicesWithinR);
				distanceAndNumPoints temp = findMaxDistanceAndNumPointsFromIndices(ConditioningEmbeddingsFromSpikes.elementAt(i), indicesWithinR,
												   ConditioningEmbeddingsFromSpikes);
				// -1 due to the point itself being in the set
				kConditioningSpikes = temp.numPoints - 1;
				radiusConditioningSpikes = temp.distance;
			}
			
			

			currentSum += (MathsUtils.digamma(kJointSpikes) - MathsUtils.digamma(kJointSamples) + 
				       ((k + l) * (-Math.log(radiusJointSpikes) + Math.log(radiusJointSamples))) -
				       MathsUtils.digamma(kConditioningSpikes) + MathsUtils.digamma(kConditioningSamples) + 
				       + (k * (Math.log(radiusConditioningSpikes) - Math.log(radiusConditioningSamples))));
			if (Double.isNaN(currentSum)) {
				throw new Exception(kJointSpikes + " " + kJointSamples + " " + kConditioningSpikes + " " + kConditioningSamples + "\n" +
						    radiusJointSpikes + " " + radiusJointSamples + " " + radiusConditioningSpikes + " " + radiusConditioningSamples);
			}
		}
		// Normalise by time
		double time_sum = 0;
		for (Double time : processTimeLengths) {
			time_sum += time;
		}
		currentSum /= time_sum;
		return currentSum;
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
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#
	 * computeSignificance(int)
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#
	 * computeSignificance(int[][])
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception {
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
