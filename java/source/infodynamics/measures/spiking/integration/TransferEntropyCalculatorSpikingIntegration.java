package infodynamics.measures.spiking.integration;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Vector;

import infodynamics.measures.spiking.TransferEntropyCalculatorSpiking;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.KdTree;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.NeighbourNodeData;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.UnivariateNearestNeighbourSearcher;

/**
 * Computes the transfer entropy between a pair of spike trains,
 *  using an integration-based measure in order to match the theoretical 
 *  form of TE between such spike trains.
 * 
 * <p>Usage paradigm is as per the interface {@link TransferEntropyCalculatorSpiking} </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class TransferEntropyCalculatorSpikingIntegration implements
		TransferEntropyCalculatorSpiking {

        protected final static boolean USE_POINT_ITSELF  = false;
        protected final static boolean TRIM_RADII  = false;
        protected final static boolean USE_SAME_RADII  = false;

	/**
	 * Number of past destination spikes to consider (akin to embedding length)
	 */
	protected int k = 1;
	/**
	 * Number of past source spikes to consider (akin to embedding length)
	 */
	protected int l = 1;

	/**
	 * Number of nearest neighbours to search for in the full joint space
	 */
	protected int Knns = 4;
	
	/**
	 * Storage for source observations supplied via {@link #addObservations(double[], double[])} etc.
	 */
	protected Vector<double[]> vectorOfSourceSpikeTimes = null;

	/**
	 * Storage for destination observations supplied via {@link #addObservations(double[], double[])} etc.
	 */
	protected Vector<double[]> vectorOfDestinationSpikeTimes = null;
		
	// constants for indexing our data storage
	protected final static int NEXT_DEST = 0;
	protected final static int NEXT_SOURCE = 1;
	protected final static int NEXT_POSSIBILITIES = 2;

	/**
	 * Cache of the timing data for each new observed spiking event in both the source
	 *  and destination
	 */
	Vector<double[][]>[] eventTimings = null;
	/**
	 * Cache of the timing data for each new observed spiking event for the
	 *  destination only
	 */

	Vector<double[]> targetEmbeddingsFromSpikes = null;
	Vector<double[]> jointEmbeddingsFromSpikes = null;
	Vector<double[]> targetEmbeddingsFromSamples = null;
	Vector<double[]> jointEmbeddingsFromSamples = null;

	protected KdTree kdTreeJointAtSpikes = null;
	protected KdTree kdTreeJointAtSamples = null;
	protected KdTree kdTreeConditioningAtSpikes = null;
	protected KdTree kdTreeConditioningAtSamples = null;
	
	Vector<double[][]> destPastAndNextTimings = null;
	/**
	 * Cache of the type of event for each new observed spiking event in both the source
	 *  and destination (i.e. which spiked next)
	 */	
	Vector<Integer> eventTypeLocator = null;
	/**
	 * Cache for each new observed spiking event of which index it has in the vector
	 *  of spiking events of the same type
	 */
	Vector<Integer> eventIndexLocator = null;
	/**
	 * Cache for each time-series of observed spiking events of how many
	 *  observations were in that set.  
	 */
	Vector<Integer> numEventsPerObservationSet = null;

	/**
	 * KdTrees for searching the joint past spaces and time to next spike,
	 *  for each possibility of which spiked next
	 */
	protected KdTree[] kdTreesJoint = null;

	/**
	 * KdTrees for searching the joint past spaces,
	 *  for each possibility of which spiked next
	 */
	protected KdTree[] kdTreesSourceDestHistories = null;

	/**
	 * KdTrees for searching the past destination space and time to next spike
	 */
	protected KdTree kdTreeDestNext = null;

	/**
	 * KdTrees for searching the past destination space
	 */
	protected KdTree kdTreeDestHistory = null;

	/**
	 * NN searcher for the time to next spike space only, if required
	 */
	protected UnivariateNearestNeighbourSearcher nnSearcherDestTimeToNextSpike = null;

	/**
	 * Property name for the number of nearest neighbours to search
	 */
	public static final String KNNS_PROP_NAME = "Knns";

	/**
	 * Property name for adjusting the search radius for the next spike such that 
	 *  it does not cover negative times (with respect to the previous spike, being either
	 *  source or destination spike)
	 */
	public static final String TRIM_TO_POS_PROP_NAME = "TRIM_RANGE_TO_POS_TIMES";
	/**
	 * Property name for an amount of random Gaussian noise to be
	 *  added to the data (default is 1e-8, matching the MILCA toolkit).
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

	protected boolean trimToPosNextSpikeTimes = false;
	
	/**
	 * Stores whether we are in debug mode
	 */
	protected boolean debug = false;
	
	public TransferEntropyCalculatorSpikingIntegration() {
		super();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#initialise(int)
	 */
	@Override
	public void initialise() throws Exception {
		initialise(k,l);
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#initialise(int)
	 */
	@Override
	public void initialise(int k) throws Exception {
		initialise(k,this.l);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#initialise(int, int)
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

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#setProperty(java.lang.String, java.lang.String)
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			l = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(KNNS_PROP_NAME)) {
			Knns = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(TRIM_TO_POS_PROP_NAME)) {
			trimToPosNextSpikeTimes = Boolean.parseBoolean(propertyValue);
	    } else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
	        if (propertyValue.equals("0") ||
	            propertyValue.equalsIgnoreCase("false")) {
	          addNoise = false;
	          noiseLevel = 0;
	        } else {
	          addNoise = true;
	          noiseLevel = Double.parseDouble(propertyValue);
	        }
		} else {
			// No property was set on this class
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#getProperty(java.lang.String)
	 */
	@Override
	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			return Integer.toString(k);
		} else if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			return Integer.toString(l);
		} else if (propertyName.equalsIgnoreCase(KNNS_PROP_NAME)) {
			return Integer.toString(Knns);
		} else if (propertyName.equalsIgnoreCase(TRIM_TO_POS_PROP_NAME)) {
			return Boolean.toString(trimToPosNextSpikeTimes);
	    } else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
	        return Double.toString(noiseLevel);
		} else {
			// No property matches for this class
			return null;
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#setObservations(double[], double[])
	 */
	@Override
	public void setObservations(double[] source, double[] destination)
			throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#startAddObservations()
	 */
	@Override
	public void startAddObservations() {
		vectorOfSourceSpikeTimes = new Vector<double[]>();
		vectorOfDestinationSpikeTimes = new Vector<double[]>();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#addObservations(double[], double[])
	 */
	@Override
	public void addObservations(double[] source, double[] destination)
			throws Exception {
		// Store these observations in our vector for now
		vectorOfSourceSpikeTimes.add(source);
		vectorOfDestinationSpikeTimes.add(destination);
	}


	
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {

		targetEmbeddingsFromSpikes = new Vector<double[]>();
		jointEmbeddingsFromSpikes = new Vector<double[]>();
		targetEmbeddingsFromSamples = new Vector<double[]>();
		jointEmbeddingsFromSamples = new Vector<double[]>();
		
		// Send all of the observations through:
		Iterator<double[]> sourceIterator = vectorOfSourceSpikeTimes.iterator();
		int timeSeriesIndex = 0;
		for (double[] destSpikeTimes : vectorOfDestinationSpikeTimes) {
			double[] sourceSpikeTimes = sourceIterator.next();
			timeSeriesIndex++;
			
			processEventsFromSpikingTimeSeries(sourceSpikeTimes, destSpikeTimes,
							   timeSeriesIndex, eventTimings, destPastAndNextTimings, 
							   eventTypeLocator, eventIndexLocator, numEventsPerObservationSet,
							   targetEmbeddingsFromSpikes, jointEmbeddingsFromSpikes,
							   targetEmbeddingsFromSamples, jointEmbeddingsFromSamples);
		}

		double[][] arrayedTargetEmbeddingsFromSpikes = new double[targetEmbeddingsFromSpikes.size()][k];
		double[][] arrayedJointEmbeddingsFromSpikes = new double[targetEmbeddingsFromSpikes.size()][k + l];
		for (int i = 0; i < targetEmbeddingsFromSpikes.size(); i++) {
			arrayedTargetEmbeddingsFromSpikes[i] = targetEmbeddingsFromSpikes.elementAt(i);
			arrayedJointEmbeddingsFromSpikes[i] = jointEmbeddingsFromSpikes.elementAt(i);
		}
		double[][] arrayedTargetEmbeddingsFromSamples = new double[targetEmbeddingsFromSamples.size()][k];
		double[][] arrayedJointEmbeddingsFromSamples = new double[targetEmbeddingsFromSamples.size()][k + l];
		for (int i = 0; i < targetEmbeddingsFromSamples.size(); i++) {
			arrayedTargetEmbeddingsFromSamples[i] = targetEmbeddingsFromSamples.elementAt(i);
			arrayedJointEmbeddingsFromSamples[i] = jointEmbeddingsFromSamples.elementAt(i);
		}

		kdTreeJointAtSpikes = new KdTree(
						 new int[] {k + l},
						 new double[][][] {arrayedJointEmbeddingsFromSpikes});
		kdTreeJointAtSamples = new KdTree(
						 new int[] {k + l},
						 new double[][][] {arrayedJointEmbeddingsFromSamples});
		kdTreeConditioningAtSpikes = new KdTree(
						 new int[] {k},
						 new double[][][] {arrayedTargetEmbeddingsFromSpikes});
		kdTreeConditioningAtSamples = new KdTree(
						 new int[] {k},
						 new double[][][] {arrayedTargetEmbeddingsFromSamples});

		/*kdTreeJointAtSpikes.setNormType("EUCLIDEAN");
		kdTreeJointAtSamples.setNormType("EUCLIDEAN");
		kdTreeConditioningAtSpikes.setNormType("EUCLIDEAN");
		kdTreeConditioningAtSamples.setNormType("EUCLIDEAN");*/
	}

	protected void makeEmbeddingsAtPoints(double[] pointsAtWhichToMakeEmbeddings, double[] sourceSpikeTimes, double[] destSpikeTimes,
					      Vector<double[]> targetEmbeddings, Vector<double[]> jointEmbeddings) {
		//System.out.println("foo");

		Random random = new Random();

		int embedding_point_index = 0;
		int most_recent_dest_index = k;
		int most_recent_source_index = l;

		// Make sure that the first point at which an embedding is made has enough preceding spikes in both source and
		// target for embeddings to be made.
		while (pointsAtWhichToMakeEmbeddings[embedding_point_index] <= destSpikeTimes[most_recent_dest_index] |
		       pointsAtWhichToMakeEmbeddings[embedding_point_index] <= sourceSpikeTimes[most_recent_source_index]) {
			embedding_point_index++;
		}

		// Loop through the points at which embeddings need to be made
		for (;embedding_point_index < pointsAtWhichToMakeEmbeddings.length; embedding_point_index++) {

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
				if (sourceSpikeTimes[most_recent_source_index + 1] < pointsAtWhichToMakeEmbeddings[embedding_point_index]) {
					most_recent_source_index++;
				} else {
					break;
				}
			}

			double[] destPast = new double[k];
			double[] jointPast = new double[k + l];
			destPast[0] = pointsAtWhichToMakeEmbeddings[embedding_point_index] -
					destSpikeTimes[most_recent_dest_index];
			jointPast[0] = pointsAtWhichToMakeEmbeddings[embedding_point_index] -
					destSpikeTimes[most_recent_dest_index];
			jointPast[k] = pointsAtWhichToMakeEmbeddings[embedding_point_index] -
					sourceSpikeTimes[most_recent_source_index];
			if (addNoise) {
				destPast[0] += random.nextGaussian()*noiseLevel;
				jointPast[0] += random.nextGaussian()*noiseLevel;
				jointPast[k] += random.nextGaussian()*noiseLevel;
			}
			for (int i = 1; i < k; i++) {
				destPast[i] = destSpikeTimes[most_recent_dest_index - i + 1] -
						destSpikeTimes[most_recent_dest_index - i];
				jointPast[i] = destSpikeTimes[most_recent_dest_index - i + 1] -
						destSpikeTimes[most_recent_dest_index - i];
				if (addNoise) {
					destPast[i] += random.nextGaussian()*noiseLevel;
					jointPast[i] += random.nextGaussian()*noiseLevel;
				}
			}
			for (int i = 1; i < l; i++) {
				jointPast[k + i] = sourceSpikeTimes[most_recent_source_index - i + 1] -
						sourceSpikeTimes[most_recent_source_index - i];
				if (addNoise) {
					jointPast[k + i] += random.nextGaussian()*noiseLevel;
				}
			}

			targetEmbeddings.add(destPast);
			jointEmbeddings.add(jointPast);
		}
	}

	protected void processEventsFromSpikingTimeSeries(double[] sourceSpikeTimes, double[] destSpikeTimes,
							  int timeSeriesIndex, Vector<double[][]>[] eventTimings, 
							  Vector<double[][]> destPastAndNextTimings, Vector<Integer> eventTypeLocator,
							  Vector<Integer> eventIndexLocator, Vector<Integer> numEventsPerObservationSet,
							  Vector<double[]> targetEmbeddingsFromSpikes, Vector<double[]> jointEmbeddingsFromSpikes,
							  Vector<double[]> targetEmbeddingsFromSamples, Vector<double[]> jointEmbeddingsFromSamples)
		throws Exception {
		// addObservationsAfterParamsDetermined(sourceSpikeTimes, destSpikeTimes);
		
		// First sort the spike times in case they were not properly in ascending order:
		Arrays.sort(sourceSpikeTimes);
		Arrays.sort(destSpikeTimes);


		// New
		int NUM_SAMPLES = sourceSpikeTimes.length;
		double sample_lower_bound = Arrays.stream(sourceSpikeTimes).min().getAsDouble();
		double sample_upper_bound = Arrays.stream(sourceSpikeTimes).max().getAsDouble();
		double[] randomSampleTimes = new double[NUM_SAMPLES];
		Random rand = new Random();
		for (int i = 0; i < randomSampleTimes.length; i++) {
			randomSampleTimes[i] = sample_lower_bound + rand.nextDouble() * (sample_upper_bound - sample_lower_bound);
		}
		Arrays.sort(randomSampleTimes);
		// End New

		makeEmbeddingsAtPoints(destSpikeTimes, sourceSpikeTimes, destSpikeTimes, targetEmbeddingsFromSpikes, jointEmbeddingsFromSpikes);
		makeEmbeddingsAtPoints(randomSampleTimes, sourceSpikeTimes, destSpikeTimes, targetEmbeddingsFromSamples, jointEmbeddingsFromSamples);
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#getAddedMoreThanOneObservationSet()
	 */
	@Override
	public boolean getAddedMoreThanOneObservationSet() {
		return (vectorOfDestinationSpikeTimes != null) &&
				(vectorOfDestinationSpikeTimes.size() > 1);
	}

	private double max_neighbour_distance(PriorityQueue<NeighbourNodeData> nnPQ) {
		double max_val = -1e9;
		while (nnPQ.peek() != null) {
			NeighbourNodeData nnData = nnPQ.poll();
			if (nnData.norms[0] > max_val) {
				max_val = nnData.norms[0];
			}
		}
		return max_val;
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#computeAverageLocalOfObservations()
	 */
	@Override
	public double computeAverageLocalOfObservations() throws Exception {

		double currentSum = 0;
		for (int i = 0; i < targetEmbeddingsFromSpikes.size(); i++) {

			PriorityQueue<NeighbourNodeData> nnPQJointSpikes =
				kdTreeJointAtSpikes.findKNearestNeighbours(Knns + 1, new double[][] {jointEmbeddingsFromSpikes.elementAt(i)});
			PriorityQueue<NeighbourNodeData> nnPQJointSamples =
				 kdTreeJointAtSamples.findKNearestNeighbours(Knns, new double[][] {jointEmbeddingsFromSpikes.elementAt(i)});
			PriorityQueue<NeighbourNodeData> nnPQConditioningSpikes =
				 kdTreeConditioningAtSpikes.findKNearestNeighbours(Knns + 1, new double[][] {targetEmbeddingsFromSpikes.elementAt(i)});
			PriorityQueue<NeighbourNodeData> nnPQConditioningSamples =
				 kdTreeConditioningAtSamples.findKNearestNeighbours(Knns, new double[][] {targetEmbeddingsFromSpikes.elementAt(i)});

			double radiusJointSpikes = max_neighbour_distance(nnPQJointSpikes);
			double radiusJointSamples = max_neighbour_distance(nnPQJointSamples);
			double radiusConditioningSpikes = max_neighbour_distance(nnPQConditioningSpikes);
			double radiusConditioningSamples = max_neighbour_distance(nnPQConditioningSamples);
			
			currentSum += ((k + l) * (- Math.log(radiusJointSpikes) + Math.log(radiusJointSamples))
				       + k * (Math.log(radiusConditioningSpikes) - Math.log(radiusConditioningSamples)));
		}
		currentSum /= (vectorOfDestinationSpikeTimes.elementAt(0)[vectorOfDestinationSpikeTimes.elementAt(0).length - 1]
			       - vectorOfDestinationSpikeTimes.elementAt(0)[0]);
		
		return currentSum;
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#computeLocalOfPreviousObservations()
	 */
	@Override
	public SpikingLocalInformationValues computeLocalOfPreviousObservations()
			throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#computeSignificance(int)
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#computeSignificance(int[][])
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#setDebug(boolean)
	 */
	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#getLastAverage()
	 */
	@Override
	public double getLastAverage() {
		// TODO Auto-generated method stub
		return 0;
	}

}
