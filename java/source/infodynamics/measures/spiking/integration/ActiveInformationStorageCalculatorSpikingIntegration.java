package infodynamics.measures.spiking.integration;

import java.util.Arrays;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Vector;
import java.util.ArrayList;
import java.util.Collections;

import infodynamics.measures.spiking.ActiveInformationStorageCalculatorSpiking;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.KdTree;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.NeighbourNodeData;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.ParsedProperties;

/**
 * Computes the active information storage for jump process based on the theoretical form of AIS
 * for spike trains.
 * 
 * <p>
 * Usage paradigm is as per the interface
 * {@link ActiveInformationStorageCalculatorSpiking}
 * </p>
 * 
 * @author Michael Fang
 */
public class ActiveInformationStorageCalculatorSpikingIntegration implements ActiveInformationStorageCalculatorSpiking {

    /**
     * The past target interspike intervals will be used for AIS past state
     * (and associated property name and convenience length variable)
     * The code assumes that the first interval is numbered 1, the next is numbered 2, etc.
     * It also assumes that the intervals are sorted. The setter method performs sorting to ensure this.
     */
    public static final String PAST_INTERVALS_PROP_NAME = "PAST_INTERVALS";
    protected int[] pastIntervals = new int[] {};
    protected int numPastIntervals = 0;

    /**
     * Number of nearest neighbours to search for in the full joint space
     */
    protected int Knns = 4;

    /**
     * Storage for target observations supplied via
     * {@link #addObservations(double[])} etc.
     */
    protected Vector<double[]> vectorOfSpikeTimes = null;

    // Data structures for embeddings
    Vector<double[]> pastEmbeddingsFromSpikes = null;
    Vector<double[]> pastEmbeddingsFromSamples = null;
    Vector<Double> processTimeLengths = null;

    // KdTrees for nearest neighbor searches
    protected KdTree kdTreePastAtSpikes = null;
    protected KdTree kdTreePastAtSamples = null;

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
     * Whether to use the jittered sampling approach. Useful for bursty spike trains.
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
     * for each marginal variable -- Options are defined by 
     * {@link KdTree#setNormType(String)} and the
     * default is {@link EuclideanUtils#NORM_EUCLIDEAN}.
     */
    public final static String PROP_NORM_TYPE = "NORM_TYPE";
    protected int normType = EuclideanUtils.NORM_EUCLIDEAN;
    
    // Last computed AIS value
    protected double lastComputedAIS = 0.0;
    
    // Whether more than one observation set was added
    protected boolean addedMoreThanOneObservationSet = false;

    public ActiveInformationStorageCalculatorSpikingIntegration() {
        super();
    }

    @Override
    public void initialise() throws Exception {
        initialise(0);
    }

    @Override
    public void initialise(int k) throws Exception {
        vectorOfSpikeTimes = null;
    }

    @Override
    public void setProperty(String propertyName, String propertyValue) throws Exception {
        boolean propertySet = true;
        if (propertyName.equalsIgnoreCase(PAST_INTERVALS_PROP_NAME)) {
            if (propertyValue.length() == 0) {
                pastIntervals = new int[] {};
            } else {
                int[] pastIntervalsTemp = ParsedProperties.parseStringArrayOfInts(propertyValue);
                for (int interval : pastIntervalsTemp) {
                    if (interval < 1) {
                        throw new Exception("Invalid interval number less than 1.");
                    }
                }
                pastIntervals = pastIntervalsTemp;
                Arrays.sort(pastIntervals);
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
                throw new Exception("Num samples multiplier must be greater than 0.");
            } else {
                numSamplesMultiplier = tempNumSamplesMultiplier;
            }
        } else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
            normType = KdTree.validateNormType(propertyValue);
        } else if (propertyName.equalsIgnoreCase(PROP_SURROGATE_SAMPLE_MULTIPLIER)) {
            double tempSurrogateNumSamplesMultiplier = Double.parseDouble(propertyValue);
            if (tempSurrogateNumSamplesMultiplier <= 0) {
                throw new Exception("Surrogate Num samples multiplier must be greater than 0.");
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

    @Override
    public void setObservations(double[] observations) throws Exception {
        startAddObservations();
        addObservations(observations);
        finaliseAddObservations();
    }

    @Override
    public void startAddObservations() {
        vectorOfSpikeTimes = new Vector<double[]>();
    }

    @Override
    public void addObservations(double[] observations) throws Exception {
        vectorOfSpikeTimes.add(observations);
    }

    @Override
    public void finaliseAddObservations() throws Exception {
        // Set convenience variables
        numPastIntervals = pastIntervals.length;

        // Initialize data structures
        pastEmbeddingsFromSpikes = new Vector<double[]>();
        pastEmbeddingsFromSamples = new Vector<double[]>();
        processTimeLengths = new Vector<Double>();

        for (double[] spikeTimes : vectorOfSpikeTimes) {
            
            // Process the spike train to generate embeddings
            processEventsFromSpikingTimeSeries(spikeTimes,
                                              pastEmbeddingsFromSpikes, 
                                              pastEmbeddingsFromSamples,
                                              numSamplesMultiplier, false);
        }

        // Convert vectors to arrays for KdTree
        double[][] arrayedPastEmbeddingsFromSpikes = new double[pastEmbeddingsFromSpikes.size()][numPastIntervals];
        for (int i = 0; i < pastEmbeddingsFromSpikes.size(); i++) {
            arrayedPastEmbeddingsFromSpikes[i] = pastEmbeddingsFromSpikes.elementAt(i);
        }

        double[][] arrayedPastEmbeddingsFromSamples = new double[pastEmbeddingsFromSamples.size()][numPastIntervals];
        for (int i = 0; i < pastEmbeddingsFromSamples.size(); i++) {
            arrayedPastEmbeddingsFromSamples[i] = pastEmbeddingsFromSamples.elementAt(i);
        }

        // Create KdTrees
        kdTreePastAtSpikes = new KdTree(arrayedPastEmbeddingsFromSpikes);
        kdTreePastAtSamples = new KdTree(arrayedPastEmbeddingsFromSamples);

        // Set norm types
        kdTreePastAtSpikes.setNormType(normType);
        kdTreePastAtSamples.setNormType(normType);
    }
    
    protected void makeEmbeddingsAtPoints(double[] pointsAtWhichToMakeEmbeddings, int indexOfFirstPointToUse,
                      double[] spikeTimes,
                      Vector<double[]> pastEmbeddings) {

        Random random = new Random();

        // Initialise the starting points of all the tracking variables
        int embeddingPointIndex = indexOfFirstPointToUse;
        // Start from the spike immediately *before* the first embedding point
        int mostRecentIndex = pastIntervals[pastIntervals.length - 1];

        // Loop through the points at which embeddings need to be made
        for (; embeddingPointIndex < pointsAtWhichToMakeEmbeddings.length; embeddingPointIndex++) {

            // Advance the tracker of the most recent spike index
            while (mostRecentIndex < (spikeTimes.length - 1)) {
                if (spikeTimes[mostRecentIndex + 1] < pointsAtWhichToMakeEmbeddings[embeddingPointIndex]) {
                    mostRecentIndex++;
                } else {
                    break;
                }
            }
            

            // Create past embedding vector
            double[] pastEmbedding = new double[numPastIntervals];

            // Add the embedding intervals from the target process
            for (int i = 0; i < pastIntervals.length; i++) {
                // Case where we are inserting an interval from an observation point back to the most recent event
                if (pastIntervals[i] == 1) {
                    pastEmbedding[i] = pointsAtWhichToMakeEmbeddings[embeddingPointIndex] - spikeTimes[mostRecentIndex];
                // Case where we are inserting an inter-event interval
                } else {
                    pastEmbedding[i] = spikeTimes[mostRecentIndex - pastIntervals[i] + 2]
                        - spikeTimes[mostRecentIndex - pastIntervals[i] + 1];
                }
            }
    
            // Add Gaussian noise if needed
            if (addNoise) {
                for (int i = 0; i < pastEmbedding.length; i++) {
                    pastEmbedding[i] = Math.log(pastEmbedding[i] + 1.1);
                    pastEmbedding[i] += random.nextGaussian() * noiseLevel;
                }
            }

            pastEmbeddings.add(pastEmbedding);
        }
    }

    /**
     * Find the first spike that has enough history to create embeddings
     */
    protected int getFirstIndex(double[] spikeTimes, boolean setProcessTimeLengths)
        throws Exception {

        // First sort the spike times in case they were not properly in ascending order
        Arrays.sort(spikeTimes);
        
        // Start from the furthest past interval needed
        int firstIndex = pastIntervals[pastIntervals.length - 1];
        
        
        // Store the time length for normalization if requested
        if (setProcessTimeLengths) {
            if (firstIndex < spikeTimes.length) {
                processTimeLengths.add(spikeTimes[spikeTimes.length - 1] - spikeTimes[firstIndex]);
            } else {
                // If no valid index was found, add a default value to avoid errors
                processTimeLengths.add(0.0);
            }
        }
        
        return firstIndex;
    }

     /**
     * Generate random sample times for estimating unconditional distributions
     */
    protected double[] generateRandomSampleTimes(double[] spikeTimes,
                     double actualNumSamplesMultiplier, int firstEmbeddingIndex, boolean doJitteredSampling) {
        
        // Define bounds for random sampling
        double sampleLowerBound = spikeTimes[firstEmbeddingIndex];
        double sampleUpperBound = spikeTimes[spikeTimes.length - 1];
        // Calculate number of samples based on the multiplier and available spikes
        int numSamples = (int) Math.round(actualNumSamplesMultiplier * (spikeTimes.length - firstEmbeddingIndex + 1));
        double[] randomSampleTimes = new double[numSamples];
        Random rand = new Random();
        
        if (doJitteredSampling) {
            // Jittered sampling: add noise to existing spike times
            for (int i = 0; i < randomSampleTimes.length; i++) {
                // Cycle through available spikes if we need more samples than spikes
                randomSampleTimes[i] = spikeTimes[firstEmbeddingIndex + 
                    (i % (spikeTimes.length - firstEmbeddingIndex - 1))] +
                    jitteredSamplingNoiseLevel * (rand.nextDouble() - 0.5);
                    
                // If the jittered time is outside our bounds, generate a uniform random time
                if ((randomSampleTimes[i] > sampleUpperBound) || (randomSampleTimes[i] < sampleLowerBound)) {
                    randomSampleTimes[i] = sampleLowerBound + rand.nextDouble() * (sampleUpperBound - sampleLowerBound);
                }
            }
        } else {
            // Uniform random sampling across the time range
            for (int i = 0; i < randomSampleTimes.length; i++) {
                randomSampleTimes[i] = sampleLowerBound + rand.nextDouble() * (sampleUpperBound - sampleLowerBound);
            }
        }
        
        // Sort the sample times in ascending order
        Arrays.sort(randomSampleTimes);
        return randomSampleTimes;
    }

    /**
     * Process spike trains to generate past embeddings
     */
    protected void processEventsFromSpikingTimeSeries(double[] spikeTimes,
                              Vector<double[]> pastEmbeddingsFromSpikes,
                              Vector<double[]> pastEmbeddingsFromSamples,
                              double actualNumSamplesMultiplier, boolean doJitteredSampling)
            throws Exception {

        // Find the first spike that has enough history to create embeddings
        int firstEmbeddingIndex = getFirstIndex(spikeTimes, true);
        
        // Generate random sample times for estimating unconditional distributions
        double[] randomSampleTimes = generateRandomSampleTimes(spikeTimes,
                               actualNumSamplesMultiplier, firstEmbeddingIndex,
                               doJitteredSampling);

        // Create embeddings at actual spike times and at random sample times
        makeEmbeddingsAtPoints(spikeTimes, firstEmbeddingIndex, spikeTimes,
                   pastEmbeddingsFromSpikes);
        makeEmbeddingsAtPoints(randomSampleTimes, 0, spikeTimes,
                   pastEmbeddingsFromSamples);
    }

     @Override
    public boolean getAddedMoreThanOneObservationSet() {
        return (vectorOfSpikeTimes != null) && (vectorOfSpikeTimes.size() > 1);
    }


    /**
     * Helper class for storing distance and point count results
     */
    protected class DistanceAndNumPoints {
        public double distance;
        public int numPoints;
        
        public DistanceAndNumPoints(double distance, int numPoints) {
            this.distance = distance;
            this.numPoints = numPoints;
        }
    }
    
    /**
     * Find maximum distance and number of points from indices
     */
    protected DistanceAndNumPoints findMaxDistanceAndNumPointsFromIndices(double[] point, int[] indices, Vector<double[]> embeddings) {
        double maxDistance = 0;
        int i = 0;
        for (; i < indices.length && indices[i] != -1; i++) {
            double distance = KdTree.norm(point, embeddings.elementAt(indices[i]), normType);
            if (distance > maxDistance) {
                maxDistance = distance;
            }
        }
        return new DistanceAndNumPoints(maxDistance, i);
    }

    @Override
    public double computeAverageLocalOfObservations() throws Exception {
        return computeAverageLocalOfObservations(kdTreePastAtSpikes, pastEmbeddingsFromSpikes);
    }
    
    /**
     * We take the actual past tree at spikes (along with the associated embeddings) as an argument, 
     * as we will need to swap these out when computing surrogates.
     */
    public double computeAverageLocalOfObservations(KdTree actualKdTreePastAtSpikes, Vector<double[]> actualPastEmbeddingsFromSpikes) throws Exception {
        double currentSum = 0;
        
        for (int i = 0; i < actualPastEmbeddingsFromSpikes.size(); i++) {
            // Find nearest neighbors in past space
            double radiusPastSpikes = actualKdTreePastAtSpikes.findKNearestNeighbours(Knns, i).poll().norms[0];
            double radiusPastSamples = kdTreePastAtSamples.findKNearestNeighbours(Knns, 
                    new double[][] { actualPastEmbeddingsFromSpikes.elementAt(i) }).poll().norms[0];
            int kPastSpikes = 0;
            int kPastSamples = 0;
            
            if (radiusPastSpikes >= radiusPastSamples) {
                kPastSpikes = Knns;
                int[] indicesWithinR = new int[pastEmbeddingsFromSamples.size() + 1];
                boolean[] isWithinR = new boolean[pastEmbeddingsFromSamples.size() + 1];
                kdTreePastAtSamples.findPointsWithinR(radiusPastSpikes,
                               new double[][] { actualPastEmbeddingsFromSpikes.elementAt(i) },
                               true,
                               isWithinR,
                               indicesWithinR);
                DistanceAndNumPoints temp = findMaxDistanceAndNumPointsFromIndices(
                                           actualPastEmbeddingsFromSpikes.elementAt(i), 
                                           indicesWithinR,
                                           pastEmbeddingsFromSamples);
                kPastSamples = temp.numPoints;
                radiusPastSamples = temp.distance;
            } else {
                kPastSamples = Knns;
                int[] indicesWithinR = new int[pastEmbeddingsFromSpikes.size() + 1];
                boolean[] isWithinR = new boolean[pastEmbeddingsFromSpikes.size() + 1];
                actualKdTreePastAtSpikes.findPointsWithinR(radiusPastSamples,
                               new double[][] { actualPastEmbeddingsFromSpikes.elementAt(i) },
                               true,
                               isWithinR,
                               indicesWithinR);
                DistanceAndNumPoints temp = findMaxDistanceAndNumPointsFromIndices(
                                           actualPastEmbeddingsFromSpikes.elementAt(i), 
                                           indicesWithinR,
                                           actualPastEmbeddingsFromSpikes);

                // -1 due to the point itself being in the set
                kPastSpikes = temp.numPoints - 1;
                radiusPastSpikes = temp.distance;
            }
            
            if (normType == EuclideanUtils.NORM_EUCLIDEAN) {
                radiusPastSpikes = Math.sqrt(radiusPastSpikes);
                radiusPastSamples = Math.sqrt(radiusPastSamples);
            }
            
            // Calculate AIS using only past space
            int pastDimension = numPastIntervals;
            double localAIS = (
                - MathsUtils.digamma(kPastSpikes) + MathsUtils.digamma(kPastSamples) + 
                pastDimension * (Math.log(radiusPastSpikes) - Math.log(radiusPastSamples))
            );

            currentSum += localAIS;
        }
        
        // Add correction factor 
        currentSum += actualPastEmbeddingsFromSpikes.size() * (
            Math.log(actualPastEmbeddingsFromSpikes.size() - 1) -
            Math.log(pastEmbeddingsFromSamples.size())
        );
        
        // Normalize by time
        double timeSum = 0;
        for (Double time : processTimeLengths) {
            timeSum += time;
        }
        currentSum /= timeSum;
        
        // Store the result
        lastComputedAIS = -currentSum;
        
        return -currentSum;
    }

    /**
     * Generate fully randomized surrogate spike sequence, breaking temporal dependencies while maintaining statistical properties
     * @param originalSpikeTrain The original spike train
     * @param random Random number generator
     * @return Randomized surrogate spike train
     */
    protected double[] generateSurrogateSpikeTrain(double[] originalSpikeTrain, Random random) {
        // Check boundary cases
        if (originalSpikeTrain == null || originalSpikeTrain.length <= 1) {
            return originalSpikeTrain == null ? null : originalSpikeTrain.clone();
        }
        
        // Maintain spike count but completely randomize positions
        double[] surrogateSpikeTrain = new double[originalSpikeTrain.length];
        double minTime = originalSpikeTrain[0];
        double maxTime = originalSpikeTrain[originalSpikeTrain.length - 1];
        
        // Generate random spike times
        for (int i = 0; i < surrogateSpikeTrain.length; i++) {
            surrogateSpikeTrain[i] = minTime + random.nextDouble() * (maxTime - minTime);
        }
        
        // Sort to ensure temporal order
        Arrays.sort(surrogateSpikeTrain);
        return surrogateSpikeTrain;
    }

    @Override
    public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck, double estimatedValue) throws Exception {
        return computeSignificance(numPermutationsToCheck, estimatedValue, System.currentTimeMillis());
    }

    @Override
    public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck, double estimatedValue, long randomSeed) throws Exception {
        Random random = new Random(randomSeed);
        double[] surrogateAISValues = new double[numPermutationsToCheck];
        
        // Save original data for later restoration
        Vector<double[]> originalPastEmbeddingsFromSpikes = pastEmbeddingsFromSpikes;
        Vector<double[]> originalPastEmbeddingsFromSamples = pastEmbeddingsFromSamples;
        Vector<Double> originalProcessTimeLengths = processTimeLengths;
        KdTree originalKdTreePastAtSpikes = kdTreePastAtSpikes;
        KdTree originalKdTreePastAtSamples = kdTreePastAtSamples;
    
        for (int permutationNumber = 0; permutationNumber < numPermutationsToCheck; permutationNumber++) {
            try {
                // Reset data structures
                pastEmbeddingsFromSpikes = new Vector<double[]>();
                pastEmbeddingsFromSamples = new Vector<double[]>();
                processTimeLengths = new Vector<Double>();
                
                // Process each original spike train
                for (double[] originalSpikeTrain : vectorOfSpikeTimes) {
                    // Create randomized surrogate spike sequence
                    double[] surrogateSpikeTrain = generateSurrogateSpikeTrain(originalSpikeTrain, random);
                    
                    // Generate embedding vectors using surrogate sequence
                    processEventsFromSpikingTimeSeries(surrogateSpikeTrain,
                                          pastEmbeddingsFromSpikes, 
                                          pastEmbeddingsFromSamples,
                                          surrogateNumSamplesMultiplier, 
                                          jitteredSamplesForSurrogates);
                }
                
                
                // Create new KdTrees
                double[][] surrogatePastArray = new double[pastEmbeddingsFromSpikes.size()][numPastIntervals];
                for (int i = 0; i < pastEmbeddingsFromSpikes.size(); i++) {
                    surrogatePastArray[i] = pastEmbeddingsFromSpikes.elementAt(i);
                }
                
                double[][] surrogateSamplesArray = new double[pastEmbeddingsFromSamples.size()][numPastIntervals];
                for (int i = 0; i < pastEmbeddingsFromSamples.size(); i++) {
                    surrogateSamplesArray[i] = pastEmbeddingsFromSamples.elementAt(i);
                }
                
                kdTreePastAtSpikes = new KdTree(surrogatePastArray);
                kdTreePastAtSpikes.setNormType(normType);
                
                kdTreePastAtSamples = new KdTree(surrogateSamplesArray);
                kdTreePastAtSamples.setNormType(normType);
                
                // Calculate surrogate AIS value
                double surrogateValue = computeAverageLocalOfObservations();
                
                // Check if calculation result is valid
                if (Double.isNaN(surrogateValue) || Double.isInfinite(surrogateValue)) {
                    throw new Exception("Invalid surrogate AIS value: " + surrogateValue);
                }
                
                surrogateAISValues[permutationNumber] = surrogateValue;
            }
            catch (Exception e) {
                // Retry this permutation, using continue instead of break
                System.out.println("Surrogate generation failed: " + e.getMessage());
                permutationNumber--;
                continue;
            }
        }
        
        // Restore original data
        pastEmbeddingsFromSpikes = originalPastEmbeddingsFromSpikes;
        pastEmbeddingsFromSamples = originalPastEmbeddingsFromSamples;
        processTimeLengths = originalProcessTimeLengths;
        kdTreePastAtSpikes = originalKdTreePastAtSpikes;
        kdTreePastAtSamples = originalKdTreePastAtSamples;
        
        // Return empirical measurement distribution
        return new EmpiricalMeasurementDistribution(surrogateAISValues, estimatedValue);
    }

    @Override
    public SpikingLocalInformationValues computeLocalOfPreviousObservations() throws Exception {
        // Implementation for local AIS values
        // Would return a custom class implementing SpikingLocalInformationValues
        // containing the local AIS values for each spike
        return null; // Placeholder
    }

    @Override
    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    @Override
    public double getLastAverage() {
        return 0;
    }


}