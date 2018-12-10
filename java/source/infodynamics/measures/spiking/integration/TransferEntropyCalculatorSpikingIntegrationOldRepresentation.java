package infodynamics.measures.spiking.integration;

import java.util.Arrays;
import java.util.Iterator;
import java.util.PriorityQueue;
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
public class TransferEntropyCalculatorSpikingIntegrationOldRepresentation implements
		TransferEntropyCalculatorSpiking {

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
	protected final static int PREV_DEST = 0;
	protected final static int PREV_SOURCE = 1;
	protected final static int PREV_POSSIBILITIES = 2;
	protected final static int NEXT_DEST = 0;
	protected final static int NEXT_SOURCE = 1;
	protected final static int NEXT_POSSIBILITIES = 2;

	/**
	 * Cache of the timing data for each new observed spiking event in both the source
	 *  and destination
	 */
	Vector<double[][]>[][] eventTimings = null;
	/**
	 * Cache of the timing data for each new observed spiking event for the
	 *  destination only
	 */
	Vector<double[][]> destPastAndNextTimings = null;
	/**
	 * Cache of the type of event for each new observed spiking event in both the source
	 *  and destination (i.e. which spiked previously, which spiked next
	 */	
	Vector<int[]> eventTypeLocator = null;
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
	 *  for each combination of which spiked previously and next
	 */
	protected KdTree[][] kdTreesJoint = null;

	/**
	 * KdTrees for searching the joint past spaces,
	 *  for each combination of which spiked previously and next
	 */
	protected KdTree[][] kdTreesSourceDestHistories = null;

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
	 * Stores whether we are in debug mode
	 */
	protected boolean debug = false;
	
	public TransferEntropyCalculatorSpikingIntegrationOldRepresentation() {
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
		// TODO Auto embed if required
		// preFinaliseAddObservations();

		// Run through each spiking time series set and pull out the observation
		//  tuples we'll store.
		// Initialise our data stores:
		eventTimings = new Vector[PREV_POSSIBILITIES][NEXT_POSSIBILITIES];
		for (int prev = 0; prev < PREV_POSSIBILITIES; prev++) {
			for (int next = 0; next < NEXT_POSSIBILITIES; next++) {
				eventTimings[prev][next] = new Vector<double[][]>();
			}
		}
		destPastAndNextTimings = new Vector<double[][]>();		
		eventTypeLocator = new Vector<int[]>();
		eventIndexLocator = new Vector<Integer>();
		numEventsPerObservationSet = new Vector<Integer>();
		
		// Send all of the observations through:
		Iterator<double[]> sourceIterator = vectorOfSourceSpikeTimes.iterator();
		int timeSeriesIndex = 0;
		for (double[] destSpikeTimes : vectorOfDestinationSpikeTimes) {
			double[] sourceSpikeTimes = sourceIterator.next();
			timeSeriesIndex++;
			
			processEventsFromSpikingTimeSeries(sourceSpikeTimes, destSpikeTimes,
					timeSeriesIndex, eventTimings, destPastAndNextTimings, 
					eventTypeLocator, eventIndexLocator, numEventsPerObservationSet);
		}
		
		// Now we have collected all the events.
		// Load up the search structures:
		// 1. Full joint space:
		// 2. Histories of source and dest only:
		kdTreesJoint = new KdTree[PREV_POSSIBILITIES][NEXT_POSSIBILITIES];
		kdTreesSourceDestHistories = new KdTree[PREV_POSSIBILITIES][NEXT_POSSIBILITIES];
		for (int prev = 0; prev < PREV_POSSIBILITIES; prev++) {
			for (int next = 0; next < NEXT_POSSIBILITIES; next++) {
				// This line does not work:
				// double[][][] jointEventTimings = (double[][][]) eventTimings[prev][next].toArray();
				// So we'll do it manually:
				double[][] sourcePastTimings = new double[eventTimings[prev][next].size()][];
				double[][] destPastTimings = new double[eventTimings[prev][next].size()][];
				double[][] nextTimings = new double[eventTimings[prev][next].size()][];
				int i = 0;
				for (double[][] timing : eventTimings[prev][next]) {
					sourcePastTimings[i] = timing[0];
					destPastTimings[i] = timing[1];
					nextTimings[i] = timing[2];
					i++;
				}
				// TODO Should we normalise before we supply to the KdTree?
				//  Think about this later. I'm not convinced it's the best
				//  approach in this particular case.
				kdTreesJoint[prev][next] = new KdTree(
						new int[] {prev == PREV_DEST ? l : l - 1,
								   prev == PREV_DEST ? k - 1 : k,
								   1},
						new double[][][] {sourcePastTimings, destPastTimings, nextTimings});
				kdTreesSourceDestHistories[prev][next] = new KdTree(
						new int[] {prev == PREV_DEST ? l : l - 1,
								   prev == PREV_DEST ? k - 1 : k},
						new double[][][] {sourcePastTimings, destPastTimings});
			}
		}
		// 3. For the dest past and time to next spike
		// 4. For the dest past only
		double[][] destPastOnlyTimings = new double[destPastAndNextTimings.size()][];
		double[][] nextTimingsForDestPastOnly = new double[destPastAndNextTimings.size()][];
		int i = 0;
		for (double[][] timing : destPastAndNextTimings) {
			destPastOnlyTimings[i] = timing[0];
			nextTimingsForDestPastOnly[i] = timing[1];
			i++;
		}
		kdTreeDestNext = new KdTree(
				new int[] {k - 1, 1},
				new double[][][] {destPastOnlyTimings, nextTimingsForDestPastOnly});
		
		if (k == 1) {
			// We need an NN searcher for the time to next spike (dest only)
			nnSearcherDestTimeToNextSpike = new UnivariateNearestNeighbourSearcher(nextTimingsForDestPastOnly);
		} else {
			kdTreeDestHistory = new KdTree(destPastOnlyTimings);
		}
	}

	protected void processEventsFromSpikingTimeSeries(double[] sourceSpikeTimes, double[] destSpikeTimes,
			int timeSeriesIndex, Vector<double[][]>[][] eventTimings, 
			Vector<double[][]> destPastAndNextTimings, Vector<int[]> eventTypeLocator,
			Vector<Integer> eventIndexLocator, Vector<Integer> numEventsPerObservationSet) throws Exception {
		// addObservationsAfterParamsDetermined(sourceSpikeTimes, destSpikeTimes);
		
		// First sort the spike times in case they were not properly in ascending order:
		Arrays.sort(sourceSpikeTimes);
		Arrays.sort(destSpikeTimes);
		
		// Scan to find the indices by which we have k and l spikes for dest and source
		//  respectively		
		int dest_index = k - 1;
		int source_index = l - 1;
		boolean previousIsDest = false;
		double[] spikeTimesForPreviousSpiker = sourceSpikeTimes;
		if (sourceSpikeTimes[source_index] > destSpikeTimes[dest_index]) {
			// Minimum required Source spikes are later than the dest.
			previousIsDest = false;
			spikeTimesForPreviousSpiker = sourceSpikeTimes;
			// Need to advance dest_index until it's the most recent before source_index
			for(;dest_index < destSpikeTimes.length; dest_index++) {
				if (destSpikeTimes[dest_index] > sourceSpikeTimes[source_index]) {
					// We've gone past the set of source spikes we have, we
					//  can move back one in the dest series
					dest_index--;
					break;
				}
			}
			if (dest_index == destSpikeTimes.length) {
				// We didn't have enough spikes in this series to generate any observations
				// TODO work out how to handle this later -- I think this is ok
				numEventsPerObservationSet.add(0);
				return;
				// throw new Exception("Dest spikes stop before enough source spikes in time-series " + timeSeriesIndex);
			}
		} else {
			// Minimum required Dest spikes are later than the source.
			previousIsDest = true;
			spikeTimesForPreviousSpiker = destSpikeTimes;
			// Need to advance source_index until it's the most recent before dest_index
			for(;source_index < sourceSpikeTimes.length; source_index++) {
				if (sourceSpikeTimes[source_index] > destSpikeTimes[dest_index]) {
					// We've gone past the set of dest spikes we have, we
					//  can move back one in the source series
					source_index--;
					break;
				}
			}
			if (source_index == sourceSpikeTimes.length) {
				// We didn't have enough spikes in this series to generate any observations
				// TODO work out how to handle this later -- I think this is ok
				numEventsPerObservationSet.add(0);
				return;
				// throw new Exception("Source spikes stop before enough dest spikes in time-series " + timeSeriesIndex);
			}
		}
		// Post-condition: dest_index and source_index are set correctly for the first set of pasts
		int indexForPreviousSpiker = previousIsDest ? dest_index : source_index;
		
		double timeToNextSpike;
		boolean nextIsDest = false;
		double[] spikeTimesForNextSpiker = sourceSpikeTimes;
		double timeOfPrevSpike = spikeTimesForPreviousSpiker[indexForPreviousSpiker];
		int numEvents = 0;
		while(true) {
			// 0. Check whether we're finished
			if ((source_index == sourceSpikeTimes.length - 1) &&
				(dest_index == destSpikeTimes.length - 1)) {
				// We have no next spike so we can't take an observation here
				//  and we're done
				break;
			}
			// Otherwise:
			// 1. Determine which of source / dest fires next
			if (source_index == sourceSpikeTimes.length - 1) {
				nextIsDest = true;
			} else if (dest_index == destSpikeTimes.length - 1) {
				nextIsDest = false;
			} else if (sourceSpikeTimes[source_index+1] < destSpikeTimes[dest_index+1]) {
				nextIsDest = false;
			} else {
				nextIsDest = true;
			}
			spikeTimesForNextSpiker = nextIsDest ? destSpikeTimes : sourceSpikeTimes;
			int indexForNextSpiker = nextIsDest ? dest_index : source_index;
			timeToNextSpike = spikeTimesForNextSpiker[indexForNextSpiker+1] - timeOfPrevSpike;
			// 2. Embed the past spikes
			double[] sourcePast = new double[previousIsDest ? l : l - 1];
			double[] destPast = new double[previousIsDest ? k - 1 : k];
			/* if (debug) {
				System.out.println("previousIsDest = " + previousIsDest + " and nextIsDest = " + nextIsDest);
			}*/
			if (previousIsDest) {
				sourcePast[0] = timeOfPrevSpike -
					sourceSpikeTimes[source_index];
			} else {
				destPast[0] = timeOfPrevSpike -
						destSpikeTimes[dest_index];
			}
			for (int i = 1; i < k; i++) {
				destPast[previousIsDest ? i - 1 : i] = destSpikeTimes[dest_index - i + 1] -
						destSpikeTimes[dest_index - i];
			}
			for (int i = 1; i < l; i++) {
				sourcePast[previousIsDest ? i : i - 1] = sourceSpikeTimes[source_index - i + 1] -
						sourceSpikeTimes[source_index - i];
			}
			// 3. Store these embedded observations
			double[][] observations = new double[][]{sourcePast, destPast,
										new double[] {timeToNextSpike}};
			if (debug) {
				System.out.printf("Adding event %d with: timeToNextSpike=%.4f, sourceSpikeTimes=", numEvents, timeToNextSpike);
				MatrixUtils.printArray(System.out, sourcePast, 3);
				System.out.printf(", destSpikeTimes=");
				MatrixUtils.printArray(System.out, destPast, 3);
				System.out.println();
			}
			// Add the index locator first so it gets the index correct before
			//  we add the new event in:
			eventIndexLocator.add(eventTimings[previousIsDest ? PREV_DEST : PREV_SOURCE][nextIsDest ? NEXT_DEST : NEXT_SOURCE].size());
			eventTimings[previousIsDest ? PREV_DEST : PREV_SOURCE][nextIsDest ? NEXT_DEST : NEXT_SOURCE].add(observations);
			eventTypeLocator.add(new int[] {previousIsDest ? PREV_DEST : PREV_SOURCE,
											nextIsDest ? NEXT_DEST : NEXT_SOURCE});
			// And finally store the observations for the dest only 
			//  search structure if required:
			if (nextIsDest) {
				double[][] destOnlyObservations;
				if (previousIsDest) {
					destOnlyObservations = new double[][] {
							destPast,
							new double[] {timeToNextSpike}
					};
				} else {
					// previous is source:
					//  We can take a copy of destPast, removing the first entry
					//  (since this only signals time the dest last fired before the source)
					//  and add that entry to the timeToNextSpike (which was back to the
					//  source firing).
					double[] destPastOnly = Arrays.copyOfRange(destPast, 1, destPast.length);
					double timeToNextSpikeSincePreviousDestSpike =
							destPast[0] + timeToNextSpike;
					destOnlyObservations = new double[][] {
							destPastOnly,
							new double[] {timeToNextSpikeSincePreviousDestSpike}
					};
				}
				destPastAndNextTimings.add(destOnlyObservations);
			}
			// 4. Reset prev as next ...
			previousIsDest = nextIsDest;
			if (previousIsDest) {
				dest_index++;
			} else {
				source_index++;
			}
			spikeTimesForPreviousSpiker = previousIsDest ? destSpikeTimes : sourceSpikeTimes;
			indexForPreviousSpiker = previousIsDest ? dest_index : source_index;
			timeOfPrevSpike = spikeTimesForPreviousSpiker[indexForPreviousSpiker];
			numEvents++;
		}
		numEventsPerObservationSet.add(numEvents);
		if (debug) {
			System.out.printf("Finished processing %d source-target events for observation set %d\n", numEvents, timeSeriesIndex);
		}
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#getAddedMoreThanOneObservationSet()
	 */
	@Override
	public boolean getAddedMoreThanOneObservationSet() {
		return (vectorOfDestinationSpikeTimes != null) &&
				(vectorOfDestinationSpikeTimes.size() > 1);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.spiking.TransferEntropyCalculatorSpiking#computeAverageLocalOfObservations()
	 */
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		
		int numberOfEvents = eventTypeLocator.size();
		
		double te = 0;
		double contributionFromSpikes = 0;
		double totalTimeLength = 0;
		
		double digammaK = MathsUtils.digamma(Knns);
		double inverseKTerm = 2.0 / (double) k;
		
		// Create temporary storage for arrays used in the neighbour counting:
		boolean[] isWithinR = new boolean[numberOfEvents]; // dummy, we don't really use this
		int[] indicesWithinR = new int[numberOfEvents];
		
		// Iterate over all the spiking events:
		Iterator<Integer> eventIndexIterator = eventIndexLocator.iterator();
		int eventIndex = -1;
		int indexForNextIsDest = -1;
		for (int[] eventType : eventTypeLocator) {
			eventIndex++;
			int eventIndexWithinType = eventIndexIterator.next().intValue();
			double[][] thisEventTimings = eventTimings[eventType[0]][eventType[1]].elementAt(eventIndexWithinType);
			totalTimeLength += thisEventTimings[2][0];
			// Find the Knns nearest neighbour matches to this event,
			//  with the same previous spiker and the next.
			// TODO Add dynamic exclusion time later
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreesJoint[eventType[0]][eventType[1]].findKNearestNeighbours(
							Knns, eventIndexWithinType);
			
			// Find eps_{x,y,z} as the maximum x, y and z norms amongst this set:
			double radius_sourcePast = 0.0;
			double radius_destPast = 0.0;
			double radius_destNext = 0.0;
			int radius_destNext_sampleIndex = -1;
			for (int j = 0; j < k; j++) {
				// Take the furthest remaining of the nearest neighbours from the PQ:
				NeighbourNodeData nnData = nnPQ.poll();
				if (nnData.norms[0] > radius_sourcePast) {
					radius_sourcePast = nnData.norms[0];
				}
				if (nnData.norms[1] > radius_destPast) {
					radius_destPast = nnData.norms[1];
				}
				if (nnData.norms[2] > radius_destNext) {
					radius_destNext = nnData.norms[2];
					radius_destNext_sampleIndex = nnData.sampleIndex;
				}
			}
			// TODO Do we need to correct radius_destNext to have a different value below, chopping it
			//  where it pushes into negative times (i.e. *before* the previous spike)?

			if (debug && (eventIndex < 10000)) {
				// Pull out the data for this observation:
				System.out.print("index = " + eventIndex + ", " +
						eventIndexWithinType + " for " +
						(eventType[0] == PREV_DEST ? "dst" : "src") +
						"->" +
						(eventType[1] == NEXT_DEST ? "dst" : "src") +
						", timings: src: ");
				MatrixUtils.printArray(System.out, thisEventTimings[0], 3);
				System.out.print(", dest: ");
				MatrixUtils.printArray(System.out, thisEventTimings[1], 3);
				System.out.print(", time to next: ");
				MatrixUtils.printArray(System.out, thisEventTimings[2], 3);
				System.out.printf("index=%d: K=%d NNs at next_range %.5f (point %d)", eventIndexWithinType, Knns, radius_destNext, radius_destNext_sampleIndex);
			}

			// Select only events where the destination spiked next:
			if (eventType[1] == NEXT_DEST) {
				indexForNextIsDest++;
				
				// Now find the matching samples in each sub-space;
				//  first match dest history and source history, with a next spike in dest:
				kdTreesSourceDestHistories[eventType[0]][NEXT_DEST].
					findPointsWithinRs(eventIndexWithinType,
							new double[] {radius_sourcePast, radius_destPast}, 0,
								true, isWithinR, indicesWithinR);
				// And check which of these samples had spike time in dest after ours:
				int countOfDestNextAndGreater = 0;
				for (int nIndex = 0; indicesWithinR[nIndex] != -1; nIndex++) {
					// Pull out this matching event from the full joint space
					double[][] matchedHistoryEventTimings = eventTimings[eventType[0]][NEXT_DEST].elementAt(indicesWithinR[nIndex]);
					if (matchedHistoryEventTimings[2][0] > thisEventTimings[2][0] + radius_destNext) {
						// This sample had a matched history and next spike was a destination
						//  spike with a longer interval than the current sample
						countOfDestNextAndGreater++;
					}
					// Reset the isWithinR array while we're here
					isWithinR[indicesWithinR[nIndex]] = false;
				}
				// And count how many samples with the matching history actually had a
				//  *source* spike next, after ours.
				// Note that we now must go to the other kdTree for next source spike
				kdTreesSourceDestHistories[eventType[0]][NEXT_SOURCE].
					findPointsWithinRs(
							new double[] {radius_sourcePast, radius_destPast}, thisEventTimings,
								true, isWithinR, indicesWithinR);
				// And check which of these samples had spike time in source at or after ours:
				int countOfSourceNextAndGreater = 0;
				for (int nIndex = 0; indicesWithinR[nIndex] != -1; nIndex++) {
					// Pull out this matching event from the full joint space
					double[][] matchedHistoryEventTimings = eventTimings[eventType[0]][NEXT_SOURCE].elementAt(indicesWithinR[nIndex]);
					if (matchedHistoryEventTimings[2][0] >= thisEventTimings[2][0] - radius_destNext) {
						// This sample had a matched history and next spike was a source
						//  spike with an interval longer than or considered equal to the current sample.
						// (The "equal to" is why we look for matches within radius_destNext here as well.)
						countOfSourceNextAndGreater++;
					}
					// Reset the isWithinR array while we're here
					isWithinR[indicesWithinR[nIndex]] = false;
				}
				
				if (debug && (eventIndex < 10000)) {
					System.out.printf(" of %d + %d + %d points with matching S-D history",
							Knns, countOfSourceNextAndGreater, countOfDestNextAndGreater);
				}
				
				// Now find the matching samples in the dest history and
				//  with a next spike timing.
				// Construct the appropriate timings to compare to here:
				double timeToNextSpikeSincePreviousDestSpike;
				if (eventType[0] == PREV_DEST) {
					timeToNextSpikeSincePreviousDestSpike = thisEventTimings[2][0];
				} else {
					// previous is source:
					timeToNextSpikeSincePreviousDestSpike =
							thisEventTimings[1][0] + thisEventTimings[2][0];
				}
				int countOfDestNextAndGreaterMatchedDest = 0;
				int countOfDestNextMatched = 0;
				if (k > 1) {
					// Search only the space of dest past -- no point
					//  searching dest past and next, since we need to run through
					//  all matches of dest past to count those with greater next spike
					//  times we might as well count those with matching spike times
					//  while we're at it.					
					kdTreeDestHistory.findPointsWithinR(indexForNextIsDest, radius_destPast,
							true, isWithinR, indicesWithinR);
					// And check which of these samples had next spike time after ours:
					for (int nIndex = 0; indicesWithinR[nIndex] != -1; nIndex++) {
						// Pull out this matching event from the dest history space
						double[][] matchedHistoryEventTimings = destPastAndNextTimings.elementAt(indicesWithinR[nIndex]);
						if (matchedHistoryEventTimings[1][0] >= timeToNextSpikeSincePreviousDestSpike - radius_destNext) {
							// This sample had a matched history and next spike was a
							//  spike with an interval longer than or considered equal to the current sample.
							// (The "equal to" is why we look for matches within kthNnData.distance here as well.)
							countOfDestNextAndGreaterMatchedDest++;
							if (matchedHistoryEventTimings[1][0] <= timeToNextSpikeSincePreviousDestSpike + radius_destNext) {
								// Then we also have a match on the next spike itself
								countOfDestNextMatched++;
							}
						}
						// Reset the isWithinR array while we're here
						isWithinR[indicesWithinR[nIndex]] = false;
					}
				} else {
					// We don't take any past dest spike ISIs into account, so we just need to look at the proportion of next
					//  spike times that match.
					countOfDestNextMatched = nnSearcherDestTimeToNextSpike.countPointsWithinOrOnR(indexForNextIsDest, radius_destNext);
					countOfDestNextAndGreaterMatchedDest = countOfDestNextMatched + 
							nnSearcherDestTimeToNextSpike.countPointsWithinROrLarger(indexForNextIsDest, radius_destNext, true);
				}
				
				if (debug && (eventIndex < 10000)) {
					System.out.printf(", and %d of %d points for D history only; ",
							countOfDestNextMatched, countOfDestNextAndGreaterMatchedDest);
				}

				// With these neighbours counted, we're ready to compute the probability of the spike given the past
				//  of source and dest.
				double logPGivenSourceAndDest = digammaK - inverseKTerm
						- MathsUtils.digamma(Knns + countOfSourceNextAndGreater + countOfDestNextAndGreater)
						+ 1.0 / (double) (Knns + countOfSourceNextAndGreater + countOfDestNextAndGreater);
				double logPGivenDest = MathsUtils.digamma(countOfDestNextMatched)
						- MathsUtils.digamma(countOfDestNextAndGreaterMatchedDest)
						+ 1.0 / ((double) countOfDestNextAndGreaterMatchedDest);
				if (debug && (eventIndex < 10000)) {
					System.out.printf(" te ~~ log (%d/%d)/(%d/%d) = %.4f -> %.4f  (inferred rates %.4f vs %.4f)\n", Knns,
							Knns + countOfSourceNextAndGreater + countOfDestNextAndGreater,
							countOfDestNextMatched, countOfDestNextAndGreaterMatchedDest,
							Math.log(((double) Knns / (double) (Knns + countOfSourceNextAndGreater + countOfDestNextAndGreater)) /
									 ((double) (countOfDestNextMatched) / (double) (countOfDestNextAndGreaterMatchedDest))), 
							logPGivenSourceAndDest - logPGivenDest,
							(double) Knns / (double) (Knns + countOfSourceNextAndGreater + countOfDestNextAndGreater) / (2.0*radius_destNext),
							(double) (countOfDestNextMatched) / (double) (countOfDestNextAndGreaterMatchedDest) / (2.0*radius_destNext));
				}
				contributionFromSpikes += logPGivenSourceAndDest - logPGivenDest;
			} else {
				if (debug) {
					System.out.println();
				}
			}
			
		}
		contributionFromSpikes /= totalTimeLength;
		te = contributionFromSpikes;
		return te;
	}

	public double computeAverageLocalOfObservationsAlg1() throws Exception {
		
		int numberOfEvents = eventTypeLocator.size();
		
		double te = 0;
		double contributionFromSpikes = 0;
		double contributionFromNonSpikes = 0;
		double contributionFromNonSpikes_destOnly = 0;
		double contributionFromNonSpikes_destAndSource = 0;
		double totalTimeLength = 0;
		
		// Create temporary storage for arrays used in the neighbour counting:
		boolean[] isWithinR = new boolean[numberOfEvents]; // dummy, we don't really use this
		int[] indicesWithinR = new int[numberOfEvents];
		
		// Iterate over all the spiking events:
		Iterator<Integer> eventIndexIterator = eventIndexLocator.iterator();
		int eventIndex = -1;
		int indexForNextIsDest = -1;
		for (int[] eventType : eventTypeLocator) {
			eventIndex++;
			int eventIndexWithinType = eventIndexIterator.next().intValue();
			double[][] thisEventTimings = eventTimings[eventType[0]][eventType[1]].elementAt(eventIndexWithinType);
			totalTimeLength += thisEventTimings[2][0];
			// Find the Knns nearest neighbour matches to this event,
			//  with the same previous spiker and the next.
			// TODO Add dynamic exclusion time later
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreesJoint[eventType[0]][eventType[1]].findKNearestNeighbours(
							Knns, eventIndexWithinType);
			// First element in the PQ is the kth NN,
			//  and epsilon = kthNnData.distance
			NeighbourNodeData kthNnData = nnPQ.poll();
			double radiusToKnn = kthNnData.distance;
			if (debug && (eventIndex < 10000)) {
				// Pull out the data for this observation:
				System.out.print("index = " + eventIndex + ", " +
						eventIndexWithinType + " for " +
						(eventType[0] == PREV_DEST ? "dst" : "src") +
						"->" +
						(eventType[1] == NEXT_DEST ? "dst" : "src") +
						", timings: src: ");
				MatrixUtils.printArray(System.out, thisEventTimings[0], 3);
				System.out.print(", dest: ");
				MatrixUtils.printArray(System.out, thisEventTimings[1], 3);
				System.out.print(", time to next: ");
				MatrixUtils.printArray(System.out, thisEventTimings[2], 3);
				System.out.printf("index=%d: K=%d NNs at range %.5f (point %d)", eventIndexWithinType, Knns, radiusToKnn, kthNnData.sampleIndex);
			}
			
			// Select only events where the destination spiked next:
			if (eventType[1] == NEXT_DEST) {
				indexForNextIsDest++;
				
				// Now find the matching samples in each sub-space;
				//  first match dest history and source history, with a next spike in dest:
				kdTreesSourceDestHistories[eventType[0]][NEXT_DEST].
					findPointsWithinR(eventIndexWithinType, radiusToKnn, 0,
						false, isWithinR, indicesWithinR);
				// And check which of these samples had spike time in dest after ours:
				int countOfDestNextAndGreater = 0;
				for (int nIndex = 0; indicesWithinR[nIndex] != -1; nIndex++) {
					// Pull out this matching event from the full joint space
					double[][] matchedHistoryEventTimings = eventTimings[eventType[0]][NEXT_DEST].elementAt(indicesWithinR[nIndex]);
					if (matchedHistoryEventTimings[2][0] >= thisEventTimings[2][0] + radiusToKnn) {
						// This sample had a matched history and next spike was a destination
						//  spike with a longer interval than the current sample
						countOfDestNextAndGreater++;
					}
					// Reset the isWithinR array while we're here
					isWithinR[indicesWithinR[nIndex]] = false;
				}
				// And count how many samples with the matching history actually had a
				//  *source* spike next, after ours.
				// Note that we now must go to the other kdTree for next source spike
				kdTreesSourceDestHistories[eventType[0]][NEXT_SOURCE].
					findPointsWithinR(radiusToKnn, thisEventTimings,
							false, isWithinR, indicesWithinR);
				// And check which of these samples had spike time in source after ours:
				int countOfSourceNextAndGreater = 0;
				for (int nIndex = 0; indicesWithinR[nIndex] != -1; nIndex++) {
					// Pull out this matching event from the full joint space
					double[][] matchedHistoryEventTimings = eventTimings[eventType[0]][NEXT_SOURCE].elementAt(indicesWithinR[nIndex]);
					if (matchedHistoryEventTimings[2][0] > thisEventTimings[2][0] - radiusToKnn) {
						// This sample had a matched history and next spike was a source
						//  spike with an interval longer than or considered equal to the current sample.
						// (The "equal to" is why we look for matches within kthNnData.distance here as well.)
						countOfSourceNextAndGreater++;
					}
					// Reset the isWithinR array while we're here
					isWithinR[indicesWithinR[nIndex]] = false;
				}
				
				if (debug && (eventIndex < 10000)) {
					System.out.printf(" of %d + %d + %d points with matching S-D history",
							Knns, countOfSourceNextAndGreater, countOfDestNextAndGreater);
				}
				
				// Now find the matching samples in the dest history and
				//  with a next spike timing.
				// Construct the appropriate timings to compare to here:
				double[][] destOnlyObservations;
				double[][] destPastOnlyObservations;
				double timeToNextSpikeSincePreviousDestSpike;
				if (eventType[0] == PREV_DEST) {
					destOnlyObservations = new double[][] {
							thisEventTimings[1], // timing of past dest spikes
							thisEventTimings[2]  // time to next spike
					};
					timeToNextSpikeSincePreviousDestSpike = thisEventTimings[2][0];
					destPastOnlyObservations = new double[][] {
							thisEventTimings[1] // timing of past dest spikes
					};
				} else {
					// previous is source:
					//  We can take a copy of the dest past timings, removing the first entry
					//  (since this only signals time the dest last fired before the source)
					//  and add that entry to the timeToNextSpike (which was back to the
					//  source firing).
					double[] destPastOnly = Arrays.copyOfRange(
							thisEventTimings[1], 1, thisEventTimings[1].length);
					timeToNextSpikeSincePreviousDestSpike =
							thisEventTimings[1][0] + thisEventTimings[2][0];
					destOnlyObservations = new double[][] {
							destPastOnly,
							new double[] {timeToNextSpikeSincePreviousDestSpike}
					};
					destPastOnlyObservations = new double[][] {
							destPastOnly
					};
				}
				int countOfDestNextAndGreaterMatchedDest = 0;
				int countOfDestNextMatched = 0;
				if (k > 1) {
					// Search only the space of dest past -- no point
					//  searching dest past and next, since we need to run through
					//  all matches of dest past to count those with greater next spike
					//  times we might as well count those with matching spike times
					//  while we're at it.
					
					// OLD WAY:
					// NO NO NO -- Can't search for it this way, because it's biased -- 
					//  should search for it by giving the index of this dest past-next 
					//  observation, so that it doesn't match to this observation.
					//  Should be able to use indexForNextIsDest here
					//kdTreeDestHistory.findPointsWithinR(radiusToKnn, destPastOnlyObservations,
					//		false, isWithinR, indicesWithinR);
					// Proper way:
					kdTreeDestHistory.findPointsWithinR(indexForNextIsDest, radiusToKnn,
							false, isWithinR, indicesWithinR);
					// And check which of these samples had next spike time after ours:
					for (int nIndex = 0; indicesWithinR[nIndex] != -1; nIndex++) {
						// Pull out this matching event from the dest history space
						double[][] matchedHistoryEventTimings = destPastAndNextTimings.elementAt(indicesWithinR[nIndex]);
						if (matchedHistoryEventTimings[1][0] > timeToNextSpikeSincePreviousDestSpike - radiusToKnn) {
							// This sample had a matched history and next spike was a
							//  spike with an interval longer than or considered equal to the current sample.
							// (The "equal to" is why we look for matches within kthNnData.distance here as well.)
							countOfDestNextAndGreaterMatchedDest++;
							if (matchedHistoryEventTimings[1][0] < timeToNextSpikeSincePreviousDestSpike + radiusToKnn) {
								// Then we also have a match on the next spike itself
								countOfDestNextMatched++;
							}
						}
						// Reset the isWithinR array while we're here
						isWithinR[indicesWithinR[nIndex]] = false;
					}
				} else {
					// We don't take any past dest spike times into account, so we just need to look at the proportion of next
					//  spike times that match.
					countOfDestNextMatched = nnSearcherDestTimeToNextSpike.countPointsStrictlyWithinR(indexForNextIsDest, radiusToKnn);
					countOfDestNextAndGreaterMatchedDest = countOfDestNextMatched + 
							nnSearcherDestTimeToNextSpike.countPointsWithinROrLarger(indexForNextIsDest, radiusToKnn, false);
				}
				
				if (debug && (eventIndex < 10000)) {
					System.out.printf(", and %d of %d points for D history only; ",
							countOfDestNextMatched, countOfDestNextAndGreaterMatchedDest);
				}

				// With these neighbours counted, we're ready to compute the probability of the spike given the past
				//  of source and dest.
				// Digammas for algorithm 1 include the extra "+1" on all terms except
				//  for the full joint space
				double logPGivenSourceAndDest = MathsUtils.digamma(Knns) -
						MathsUtils.digamma(Knns + countOfSourceNextAndGreater + countOfDestNextAndGreater + 1);
				double logPGivenDest = MathsUtils.digamma(countOfDestNextMatched + 1) -
						MathsUtils.digamma(countOfDestNextAndGreaterMatchedDest + 1);
				if (debug && (eventIndex < 10000)) {
					System.out.printf(" te ~~ log (%d/%d)/(%d/%d) = %.4f -> %.4f  (inferred rates %.4f vs %.4f)\n", Knns,
							Knns + countOfSourceNextAndGreater + countOfDestNextAndGreater + 1,
							countOfDestNextMatched + 1, countOfDestNextAndGreaterMatchedDest + 1,
							Math.log(((double) Knns / (double) (Knns + countOfSourceNextAndGreater + countOfDestNextAndGreater + 1)) /
									 ((double) (countOfDestNextMatched + 1) / (double) (countOfDestNextAndGreaterMatchedDest + 1))), 
							logPGivenSourceAndDest - logPGivenDest,
							(double) Knns / (double) (Knns + countOfSourceNextAndGreater + countOfDestNextAndGreater + 1) / (2.0*radiusToKnn),
							(double) (countOfDestNextMatched + 1) / (double) (countOfDestNextAndGreaterMatchedDest + 1) / (2.0*radiusToKnn));
				}
				contributionFromSpikes += logPGivenSourceAndDest - logPGivenDest;
			} else {
				if (debug) {
					System.out.println();
				}
			}
			
			// Regardless of which type of event it was, we need to integrate
			//  the spiking rates up until the next spiking event
			// Our first attempt at a solution uses the search width defined
			//  using the history and the next spike.
			
			// Consider first the destination process only.
			// Match dest history
			double[][] destPastOnlyObservations;
			double timeToNextSpikeSincePreviousDestSpike;
			if (eventType[0] == PREV_DEST) {
				timeToNextSpikeSincePreviousDestSpike = thisEventTimings[2][0];
				destPastOnlyObservations = new double[][] {
						thisEventTimings[1] // timing of past dest spikes
				};
			} else {
				// previous is source:
				//  We can take a copy of the dest past timings, removing the first entry
				//  (since this only signals time the dest last fired before the source)
				//  and add that entry to the timeToNextSpike (which was back to the
				//  source firing).
				double[] destPastOnly = Arrays.copyOfRange(
						thisEventTimings[1], 1, thisEventTimings[1].length);
				timeToNextSpikeSincePreviousDestSpike =
						thisEventTimings[1][0] + thisEventTimings[2][0];
				destPastOnlyObservations = new double[][] {
						destPastOnly
				};
			}
			int countOfDestNextEarlier = 0;
			int countOfDestMatches = 0;
			if (k > 1) {
				kdTreeDestHistory.findPointsWithinR(radiusToKnn, destPastOnlyObservations,
					false, isWithinR, indicesWithinR);
				// And check which of these samples had next spike time before ours:
				for (int nIndex = 0; indicesWithinR[nIndex] != -1; nIndex++) {
					// Pull out this matching event from the dest history space
					double[][] matchedHistoryEventTimings = destPastAndNextTimings.elementAt(indicesWithinR[nIndex]);
					if (matchedHistoryEventTimings[1][0] < timeToNextSpikeSincePreviousDestSpike) {
						// This sample had a matched history and next spike was a
						//  spike with an interval shorter than the current sample.
						countOfDestNextEarlier++;
					}
					// Reset the isWithinR array while we're here
					isWithinR[indicesWithinR[nIndex]] = false;
					countOfDestMatches++;
				}
			} else {
				// We're not using the past, so we match on everything up to the last spike
				countOfDestMatches = nnSearcherDestTimeToNextSpike.getNumObservations();
				countOfDestNextEarlier = nnSearcherDestTimeToNextSpike.countPointsSmallerAndOutsideR(
						// indexForNextIsDest must point to the next event (possibly this one) where the dest spikes next.
						(eventType[1] == NEXT_DEST) ? indexForNextIsDest : indexForNextIsDest + 1,
						radiusToKnn, false);

			}
			// And include the contribution for each of these
			double integralForDestHistorySpace = 0;
			for (int hi = 0; hi < countOfDestNextEarlier; hi++) {
				integralForDestHistorySpace += (double) 1 /
							(double) (countOfDestMatches - hi);
			}
			
			// First match dest history and source history, with a next spike in dest:
			kdTreesSourceDestHistories[eventType[0]][NEXT_DEST].
				findPointsWithinR(radiusToKnn, thisEventTimings,
						false, isWithinR, indicesWithinR);
			// And store which of these samples had spike time in dest before ours.
			//  Store them in a vector of double arrays, with each array holding the
			//  spike time then -1 for a next dest spike and +1 for a source spike
			Vector<double[]> spikesBeforeOurs = new Vector<double[]>();
			int countOfSpikesAfterAndIncludingOurs = 0;
			for (int nIndex = 0; indicesWithinR[nIndex] != -1; nIndex++) {
				// Pull out this matching event from the full joint space
				double[][] matchedHistoryEventTimings = eventTimings[eventType[0]][NEXT_DEST].elementAt(indicesWithinR[nIndex]);
				if (matchedHistoryEventTimings[2][0] < thisEventTimings[2][0]) {
					// This sample had a matched history and next spike was a destination
					//  spike with a shorted interval than the current sample
					spikesBeforeOurs.add(new double[] {
							matchedHistoryEventTimings[2][0], -1});
				} else {
					countOfSpikesAfterAndIncludingOurs++;
				}
				// Reset the isWithinR array while we're here
				isWithinR[indicesWithinR[nIndex]] = false;
			}
			// And store which of these samples had spike time in source before ours.
			//  Store them in a vector of double arrays, with each array holding the
			//  spike time then -1 for a next dest spike and +1 for a source spike
			// Note that we now must go to the other kdTree for next source spike
			kdTreesSourceDestHistories[eventType[0]][NEXT_SOURCE].
				findPointsWithinR(radiusToKnn, thisEventTimings,
						false, isWithinR, indicesWithinR);
			for (int nIndex = 0; indicesWithinR[nIndex] != -1; nIndex++) {
				// Pull out this matching event from the full joint space
				double[][] matchedHistoryEventTimings = eventTimings[eventType[0]][NEXT_SOURCE].elementAt(indicesWithinR[nIndex]);
				if (matchedHistoryEventTimings[2][0] < thisEventTimings[2][0]) {
					// This sample had a matched history and next spike was a source
					//  spike with a shorted interval than the current sample
					spikesBeforeOurs.add(new double[] {
							matchedHistoryEventTimings[2][0], +1});
				} else {
					countOfSpikesAfterAndIncludingOurs++;
				}
				// Reset the isWithinR array while we're here
				isWithinR[indicesWithinR[nIndex]] = false;
			}
			// Now we can sort the spikes which occur before ours and process
			//  them in order:
			// Next line doesn't work, so replaced with clunkier code:
			// double[][] nextSpikeTimesAndType = (double[][]) spikesBeforeOurs.toArray();
			double[][] nextSpikeTimesAndType = new double[spikesBeforeOurs.size()][];
			for (int si = 0; si < nextSpikeTimesAndType.length; si++) {
				nextSpikeTimesAndType[si] = spikesBeforeOurs.elementAt(si);
			}
			Arrays.sort(nextSpikeTimesAndType, FirstIndexComparatorDouble.getInstance());
			double integralForJointSpace = 0;
			for (int si = 0; si < nextSpikeTimesAndType.length; si++) {
				if (nextSpikeTimesAndType[si][1] < 0) {
					// We have a next spike from the dest, which is 
					//  earlier than our spike.
					// Integrated Prob for getting a spike here is 1 / N, where
					//  N is the number of properly matched histories (i.e.
					//  which don't have a next spike before this one)
					double intPNext = (double) 1 / (double) 
							(nextSpikeTimesAndType.length - si + countOfSpikesAfterAndIncludingOurs);
					integralForJointSpace += intPNext;
				}
				// Ignore next spikes on the source, they simply get removed
				//  from the matched histories count
			}
			
			// We now have the integral of spike rates given the dest and
			// joint histories, so subtract this out:
			contributionFromNonSpikes += integralForDestHistorySpace - integralForJointSpace;
			contributionFromNonSpikes_destAndSource += integralForJointSpace;
			contributionFromNonSpikes_destOnly += integralForDestHistorySpace;
		}
		contributionFromSpikes /= totalTimeLength;
		contributionFromNonSpikes /= totalTimeLength;
		contributionFromNonSpikes_destAndSource /= totalTimeLength;
		contributionFromNonSpikes_destOnly /= totalTimeLength;
		te = contributionFromSpikes + contributionFromNonSpikes;
		System.out.printf("TE = %.4f (spikes) + %.4f (non-spikes: d:%.4f - s-d:%.4f) = %.4f\n",
				contributionFromSpikes, contributionFromNonSpikes,
				contributionFromNonSpikes_destOnly, contributionFromNonSpikes_destAndSource, te);
		return te;
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
