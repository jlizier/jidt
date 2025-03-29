/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.measures.continuous;

import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;
import java.util.Vector;

/**
 * Implements {@link MutualInfoCalculatorMultiVariate} to provide a base
 * class with common functionality for child class implementations of
 * {@link MutualInfoCalculatorMultiVariate}
 * via various estimators. 
 * 
 * <p>These various estimators include: e.g. box-kernel estimation, KSG estimators, etc
 * (see the child classes linked above).
 * </p>
 * 
 * <p>Usage is as outlined in {@link MutualInfoCalculatorMultiVariate}.</p>
 * 
  * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class MutualInfoMultiVariateCommon implements
		MutualInfoCalculatorMultiVariate {

	/**
	 * Number of dimenions for our source multivariate data set
	 */
	protected int dimensionsSource = 1;
	/**
	 * Number of dimenions for our destination multivariate data set
	 */
	protected int dimensionsDest = 1;
	
	/**
	 * The set of source observations, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link addObservations(double[][], double[][])} functions.
	 */
	protected double[][] sourceObservations;
	/**
	 * Stored means of source variables before any normalising or noise addition
	 */
	protected double[] sourceMeansBeforeNorm;
	/**
	 * Stored standard deviations of source variables before any normalising or noise addition
	 */
	protected double[] sourceStdsBeforeNorm;

	/**
	 * The set of destination observations, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link addObservations(double[][], double[][])} functions.
	 */
	protected double[][] destObservations;
	/**
	 * Stored means of destination variables before any normalising or noise addition
	 */
	protected double[] destMeansBeforeNorm;
	/**
	 * Stored standard deviations of destination variablesbefore any normalising or noise addition
	 */
	protected double[] destStdsBeforeNorm;

	/**
	 * Track which observation set each sample came from
	 */
	protected int[] observationSetIndices;

	/**
	 * Track which sample index within an observation set that each sample came from
	 */
	protected int[] observationTimePoints;

	/**
	 * Total number of observations supplied.
	 * Only valid after {@link #finaliseAddObservations()} is called.
	 */
	protected int totalObservations = 0;

	/**
	 * Store the last computed average MI
	 */
	protected double lastAverage;

	/**
	 * Track whether we've computed the average for the supplied
	 *  observations yet
	 */
	protected boolean miComputed;

	/**
	 * Whether to report debug messages or not
	 */
	protected boolean debug;

	/**
	 * Time difference from the source to the destination observations
	 *  (ie destination lags the source by this time:
	 *  	we compute I(source_{n}; dest_{n+timeDiff}).
	 * (Note that our internal sourceObservations and destObservations
	 *  are adjusted so that there is no timeDiff between them).
	 */
	protected int timeDiff = 0;

	/**
	 * Storage for source observations supplied via {@link #addObservations(double[][], double[][])}
	 * type calls
	 */
	protected Vector<double[][]> vectorOfSourceObservations;
	/**
	 * Storage for destination observations supplied via {@link #addObservations(double[][], double[][])}
	 * type calls
	 */
	protected Vector<double[][]> vectorOfDestinationObservations;
	/**
	 * Tracks separate (time-series) observation sets
	 *  we are taking samples from
	 */
	protected int observationSetIndex = 0;
	/**
	 * Storage for which observation set each
	 *  block of samples comes from
	 */
	protected Vector<Integer> vectorOfObservationSetIndices;
	/**
	 * Storage for start time point for the observation
	 *  set within its block of samples
	 */
	protected Vector<Integer> vectorOfObservationStartTimePoints;

	/**
	 * Whether the user has supplied more than one (disjoint) set of samples
	 */
	protected boolean addedMoreThanOneObservationSet;
	/**
	 * Whether to normalise the incoming data 
	 */
	protected boolean normalise = true;
	/**
	 * Which strategy to choose to select the surrogate data 
	 */
	protected String surrogate_type = PROP_SHUFFLE;
	/**
	 * Whether to add an amount of random noise to the incoming data
	 */
	protected boolean addNoise = false;
	/**
	 * Amount of random Gaussian noise to add to the incoming data.
	 * 0 by default except for KSG estimators (where it is recommended
	 *  and 1e-8 is used to match MILCA toolkit)
	 */
	protected double noiseLevel = (double) 0.0;
	/**
	 * Whether we use dynamic correlation exclusion
	 */
	protected boolean dynCorrExcl = false;
	/**
	 * Size of dynamic correlation exclusion window.
	 */
	protected int dynCorrExclTime = 0;
	 * Has the user set a seed for the random noise
	 */
	protected boolean noiseSeedSet = false;
	/**
	 * Seed that the user set for the random noise
	 */
	protected long noiseSeed = 0;
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorCommon#initialise()
	 */
	public void initialise() throws Exception {
		initialise(dimensionsSource, dimensionsDest);
	}

	public void initialise(int sourceDimensions, int destDimensions) {
		dimensionsSource = sourceDimensions;
		dimensionsDest = destDimensions;
		lastAverage = 0.0;
		totalObservations = 0;
		miComputed = false;
		sourceObservations = null;
		destObservations = null;
		sourceMeansBeforeNorm = null;
		sourceStdsBeforeNorm = null;
		destMeansBeforeNorm = null;
		destStdsBeforeNorm = null;
		observationSetIndices = null;
		observationTimePoints = null;
		observationSetIndex = 0;
		vectorOfObservationSetIndices = null;
		vectorOfObservationStartTimePoints = null;
		addedMoreThanOneObservationSet = false;
	}

	/**
	 * Set properties for the calculator.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} --
	 *  		Time difference between source and destination (0 by default).
	 * 			Must be >= 0.
	 * 		</li>
	 *  <li>{@link #PROP_NORMALISE} -- whether to normalise the incoming individual
	 *      variables to mean 0 and standard deviation 1 (true by default, except for Gaussian calculator)</li>
	 *  <li>{@link #PROP_SURROGATE_TYPE} -- Which strategy to choose to select the surrogate data. Valid values are 
	 *		"SHUFFLE" - randomly suffling data, or "ROTATE", keeping data consecutive but shuffling orders. Default is shuffle.</li>
	 *  <li>{@link #PROP_ADD_NOISE} -- a standard deviation for an amount of
	 *    random Gaussian noise to add to
	 *      each variable, to avoid having neighbourhoods with artificially
	 *      large counts. (We also accept "false" to indicate "0".)
	 *      The amount is added in after any normalisation,
	 *      so can be considered as a number of standard deviations of the data.
	 *      Default is 0 for most estimators; this is strongly recommended by
	 *      by Kraskov for the KSG method though, so for that estimator we 
	 *      use 1e-8 to match the MILCA toolkit (though note it adds in
	 *      a random amount of noise in [0,noiseLevel) ).</li>
	 *  <li>{@link #PROP_NOISE_SEED} -- a long value seed for the random noise generator or
	 *      the string {@link MutualInfoCalculatorMultiVariate#NOISE_NO_SEED_VALUE} for no seed (default)</li>
	 * </ul>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_TIME_DIFF)) {
			timeDiff = Integer.parseInt(propertyValue);
			if (timeDiff < 0) {
				throw new Exception("Time difference must be >= 0. Flip data1 and data2 around if required.");
			}
		} else if (propertyName.equalsIgnoreCase(PROP_DYN_CORR_EXCL_TIME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);	
	    } else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
	        normalise = Boolean.parseBoolean(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_SURROGATE_TYPE)){
			if (Arrays.asList(VALID_SURROGATE_TYPES).contains(propertyValue.toUpperCase())) {
				surrogate_type = propertyValue.toUpperCase();
			} else {
				throw new IllegalArgumentException("Invalid value for surrogate_type. Allowed values are: " + VALID_SURROGATE_TYPES);
			}
	    } else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
	        if (propertyValue.equals("0") ||
	            propertyValue.equalsIgnoreCase("false")) {
	          addNoise = false;
	          noiseLevel = 0;
	        } else {
	          addNoise = true;
	          noiseLevel = Double.parseDouble(propertyValue);
	        }
	    } else if (propertyName.equalsIgnoreCase(PROP_NOISE_SEED)) {
	    	if (propertyValue.equals(NOISE_NO_SEED_VALUE)) {
	    		noiseSeedSet = false;
	    	} else {
	    		noiseSeedSet = true;
	    		noiseSeed = Long.parseLong(propertyValue);
	    	}
		} else {
			// No property was set here
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public String getProperty(String propertyName)
			throws Exception {
		
		if (propertyName.equalsIgnoreCase(PROP_TIME_DIFF)) {
			return Integer.toString(timeDiff);
		} else if (propertyName.equalsIgnoreCase(PROP_DYN_CORR_EXCL_TIME)) {
			return Integer.toString(dynCorrExclTime);
	    } else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
	        return Boolean.toString(normalise);
	    } else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
	        return Double.toString(noiseLevel);
  		} else if (propertyName.equalsIgnoreCase(PROP_SURROGATE_TYPE)){
			  return surrogate_type;
	    } else if (propertyName.equalsIgnoreCase(PROP_NOISE_SEED)) {
	    	if (noiseSeedSet) {
	    		return Long.toString(noiseSeed);
	    	} else {
	    		return NOISE_NO_SEED_VALUE;
	    	}
    } else {
			// No property was recognised here
			throw new Exception("Unknown property to get: " + propertyName);
		}
	}

	@Override
	public void setObservations(double[][] source, double[][] destination) throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}

	/**
	 * A non-overloaded method signature for setObservations with 2D arguments, as there have been
	 *  some problems calling overloaded versions of setObservations from jpype. 
	 * 
	 * @param source
	 * @param destination
	 * @throws Exception
	 */
	public void setObservations2D(double[][] source, double[][] destination) throws Exception {
		setObservations(source, destination);
	}

	@Override
	public void setObservations(double[] source, double[] destination) throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}

	/**
	 * A non-overloaded method signature for setObservations with 1D arguments, as there have been
	 *  some problems calling overloaded versions of setObservations from jpype. 
	 * 
	 * @param source
	 * @param destination
	 * @throws Exception
	 */
	public void setObservations1D(double[] source, double[] destination) throws Exception {
		setObservations(source, destination);
	}
	
	@Override
	public void setObservations(double[] source, double[][] destination) throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}

	@Override
	public void setObservations(double[][] source, double[] destination) throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}

	@Override
	public void setObservations(double[][] source, double[][] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		startAddObservations();
		addObservations(source, destination, sourceValid, destValid);
		finaliseAddObservations();
	}
	
	@Override
	public void setObservations(double[] source, double[] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {

		if ((dimensionsDest != 1) || (dimensionsSource != 1)) {
			throw new Exception("The number of source and dest dimensions (having been initialised to " +
					dimensionsSource + " and " + dimensionsDest + ") can only be 1 when " +
					"the univariate addObservations(double[],double[]) and " + 
					"setObservations(double[],double[]) methods are called");
		}
		setObservations(MatrixUtils.reshape(source, source.length, 1),
				MatrixUtils.reshape(destination, destination.length, 1),
				sourceValid, destValid);
	}
	
	@Override
	public void setObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception {

		boolean[] allSourceValid = MatrixUtils.andRows(sourceValid);
		boolean[] allDestValid = MatrixUtils.andRows(destValid);
		setObservations(source, destination, allSourceValid, allDestValid);
	}

	@Override
	public void startAddObservations() {
		vectorOfSourceObservations = new Vector<double[][]>();
		vectorOfDestinationObservations = new Vector<double[][]>();
		vectorOfObservationSetIndices = new Vector<Integer>();
		vectorOfObservationStartTimePoints = new Vector<Integer>();
	}
	
	@Override
	public void addObservations(double[][] source, double[][] destination) throws Exception {
		// Use the current observationSetIndex and increment for next use:
		addObservationsTrackObservationIDs(source, destination, observationSetIndex++, 0);
	}
	
	public void addObservationsTrackObservationIDs(double[][] source, double[][] destination,
			int observationSetIndexToUse, int startTimeIndex) throws Exception {
		if (vectorOfSourceObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length <= timeDiff) {
			// we won't be taking any observations here
			return;
		}
		if (source[0].length != dimensionsSource) {
			throw new Exception(String.format("Number of joint variables in source data %d " +
					"does not match the initialised value %d", source[0].length, dimensionsSource));
		}
		if (destination[0].length != dimensionsDest) {
			throw new Exception(String.format("Number of joint variables in destination data %d " +
					"does not match the initialised value %d", destination[0].length, dimensionsDest));
		}
		vectorOfSourceObservations.add(source);
		vectorOfDestinationObservations.add(destination);
		vectorOfObservationSetIndices.add(observationSetIndexToUse);
		vectorOfObservationStartTimePoints.add(startTimeIndex);
	}

	/**
	 * A non-overloaded method signature for addObservations with 2D arguments, as there have been
	 *  some problems calling overloaded versions of setObservations from jpype. 
	 * 
	 * @param source
	 * @param destination
	 * @throws Exception
	 */
	public void addObservations2D(double[][] source, double[][] destination) throws Exception {
		addObservations(source, destination);
	}

	@Override
	public void addObservations(double[] source, double[] destination) throws Exception {
		
		if ((dimensionsDest != 1) || (dimensionsSource != 1)) {
			throw new Exception("The number of source and dest dimensions (having been initialised to " +
					dimensionsSource + " and " + dimensionsDest + ") can only be 1 when " +
					"the univariate addObservations(double[],double[]) and " + 
					"setObservations(double[],double[]) methods are called");
		}
		addObservations(MatrixUtils.reshape(source, source.length, 1),
						MatrixUtils.reshape(destination, destination.length, 1));
	}

	/**
	 * A non-overloaded method signature for addObservations with 1D arguments, as there have been
	 *  some problems calling overloaded versions of setObservations from jpype. 
	 * 
	 * @param source
	 * @param destination
	 * @throws Exception
	 */
	public void addObservations1D(double[] source, double[] destination) throws Exception {
		addObservations(source, destination);
	}

	@Override
	public void addObservations(double[] source, double[][] destination) throws Exception {
		
		if (dimensionsSource != 1) {
			throw new Exception("The number of source dimensions (having been initialised to " +
					dimensionsSource + ") can only be 1 when " +
					"the partially univariate addObservations(double[],double[][]) and " + 
					"setObservations(double[],double[][]) methods are called");
		}
		addObservations(MatrixUtils.reshape(source, source.length, 1),
						destination);
	}
	
	@Override
	public void addObservations(double[][] source, double[] destination) throws Exception {
		
		if (dimensionsDest != 1) {
			throw new Exception("The number of dest dimensions (having been initialised to " +
					dimensionsDest + ") can only be 1 when " +
					"the partially univariate addObservations(double[][],double[]) and " + 
					"setObservations(double[][],double[]) methods are called");
		}
		addObservations(source,
				MatrixUtils.reshape(destination, destination.length, 1));
	}
	
	@Override
	public void addObservations(double[][] source, double[][] destination,
			int startTime, int numTimeSteps) throws Exception {
		addObservations(source, destination, startTime, numTimeSteps, observationSetIndex++);
	}
	
	protected void addObservations(double[][] source, double[][] destination,
			int startTime, int numTimeSteps, int observationSetIndexToUse) throws Exception {
		if (vectorOfSourceObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (numTimeSteps <= timeDiff) {
			// We won't be taking any observations here
			return;
		}
		double[][] sourceToAdd = new double[numTimeSteps][];
		System.arraycopy(source, startTime, sourceToAdd, 0, numTimeSteps);
		double[][] destToAdd = new double[numTimeSteps][];
		System.arraycopy(destination, startTime, destToAdd, 0, numTimeSteps);
		addObservationsTrackObservationIDs(sourceToAdd, destToAdd, observationSetIndexToUse, startTime);
	}

	@Override
	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception {
		
		if ((dimensionsDest != 1) || (dimensionsSource != 1)) {
			throw new Exception("The number of source and dest dimensions (having been initialised to " +
					dimensionsSource + " and " + dimensionsDest + ") can only be 1 when " +
					"the univariate addObservations(double[],double[]) and " + 
					"setObservations(double[],double[]) methods are called");
		}
		addObservations(MatrixUtils.reshape(source, source.length, 1),
						MatrixUtils.reshape(destination, destination.length, 1),
						startTime, numTimeSteps);
	}

	public void addObservations(double[] source, double[] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		addObservations(MatrixUtils.reshape(source, source.length, 1),
				MatrixUtils.reshape(destination, destination.length, 1),
				sourceValid, destValid);
	}
	
	public void addObservations(double[][] source, double[][] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		
		Vector<int[]> startAndEndTimePairs = computeStartAndEndTimePairs(sourceValid, destValid);
		
		// We've found the set of start and end times for this pair
		startAddObservations();
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservations(source, destination, startTime, endTime - startTime + 1, observationSetIndex);
		}
		observationSetIndex++;
		finaliseAddObservations();
	}

	public void addObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		
		boolean[] allSourceValid = MatrixUtils.andRows(sourceValid);
		boolean[] allDestValid = MatrixUtils.andRows(destValid);
		addObservations(source, destination, allSourceValid, allDestValid);
	}

	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
	 * 
	 * <p>This default implementation simply puts all of the observations into
	 *  the {@link #sourceObservations} and {@link #destObservations} arrays.
	 * Usually child implementations will override this, call this implementation
	 *  to perform the common processing, then perform their own processing.
	 * </p>
	 * 
	 * @throws Exception Allow child classes to throw an exception if there
	 *  is an issue detected specific to that calculator.
	 */
	public void finaliseAddObservations() throws Exception {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[][] destination : vectorOfDestinationObservations) {
			totalObservations += destination.length - timeDiff;
		}
		destObservations = new double[totalObservations][dimensionsDest];
		sourceObservations = new double[totalObservations][dimensionsSource];
		observationSetIndices = new int[totalObservations];
		observationTimePoints = new int[totalObservations];
		
		// Construct the joint vectors from the given observations
		//  (removing redundant data which is outside any timeDiff)
		int startObservation = 0;
		Iterator<double[][]> iterator = vectorOfDestinationObservations.iterator();
		Iterator<Integer> iteratorObsSetIndices = vectorOfObservationSetIndices.iterator();
		Iterator<Integer> iteratorObsStartTimePoints = vectorOfObservationStartTimePoints.iterator();
		for (double[][] source : vectorOfSourceObservations) {
			double[][] destination = iterator.next();
			// Copy the data from these given observations into our master 
			//  array, aligning them incorporating the timeDiff:
			MatrixUtils.arrayCopy(source, 0, 0,
					sourceObservations, startObservation, 0,
					source.length - timeDiff, dimensionsSource);
			MatrixUtils.arrayCopy(destination, timeDiff, 0,
					destObservations, startObservation, 0,
					destination.length - timeDiff, dimensionsDest);
			int numNewObservations = destination.length - timeDiff;
			// And update which observation set and time index each sample came from:
			Arrays.fill(observationSetIndices, startObservation, startObservation + numNewObservations, iteratorObsSetIndices.next());
			int firstTimeSampleId = iteratorObsStartTimePoints.next(); // This is the first sample of the source; storing the time index for destination below.
			for (int i = 0; i < numNewObservations; i++) {
				observationTimePoints[startObservation + i] = firstTimeSampleId + timeDiff + i;
			}
			startObservation += numNewObservations;
		}
		if (vectorOfSourceObservations.size() > 1) {
			addedMoreThanOneObservationSet = true;
		}
		// We don't need to keep the vectors of observation sets anymore:
		vectorOfSourceObservations = null;
		vectorOfDestinationObservations = null;
		
		// Normalise the data if required, and store means/stds any normalising
		sourceMeansBeforeNorm = MatrixUtils.means(sourceObservations);
		sourceStdsBeforeNorm = MatrixUtils.stdDevs(sourceObservations, sourceMeansBeforeNorm);
		destMeansBeforeNorm = MatrixUtils.means(destObservations);
		destStdsBeforeNorm = MatrixUtils.stdDevs(destObservations, destMeansBeforeNorm);
		if (normalise) {
			normaliseData();
		}

		// Add Gaussian noise of std dev noiseLevel to the data if required
		if (addNoise) {
			Random random = new Random();
	    	if (noiseSeedSet) {
	    		random.setSeed(noiseSeed);
	    	}
			for (int r = 0; r < sourceObservations.length; r++) {
				for (int c = 0; c < dimensionsSource; c++) {
					sourceObservations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
				for (int c = 0; c < dimensionsDest; c++) {
					destObservations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
			}
		}
	}
	
	/**
	 * Protected method to normalise the stored data samples for each variable.
	 * This method can be overriden by children if required to perform
	 * specific actions for their estimation methods.
	 * Assumes that the member variables for means and stds have already been computed.
	 */
	protected void normaliseData() {
		// We can overwrite these since they're already
		//  a copy of the users' data.
		MatrixUtils.normalise(sourceObservations, sourceMeansBeforeNorm, sourceStdsBeforeNorm);
		MatrixUtils.normalise(destObservations, destMeansBeforeNorm, destStdsBeforeNorm);
	}
	
	/**
	 * Generate a resampled distribution of what the MI would look like,
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination value.
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI.
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p>This method (in contrast to {@link #computeSignificance(int[][])})
	 * creates <i>random</i> shufflings of the next values for the surrogate MI
	 * calculations.</p>
	 * 
	 * @param numSurrogatesToCheck number of surrogate samples for permutations
	 *  to generate the distribution.
	 * @return the distribution of channel measure scores under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numSurrogatesToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = new int[numSurrogatesToCheck][sourceObservations.length];
		
		if (surrogate_type.equalsIgnoreCase(PROP_ROTATE)){
			newOrderings = rg.generateRotatedSurrogates(sourceObservations.length, numSurrogatesToCheck, dynCorrExclTime);
		} else {
			newOrderings = rg.generateRandomPerturbations(sourceObservations.length, numSurrogatesToCheck);
		}
		return computeSignificance(newOrderings);
	}

	/**
	 * Generate a resampled distribution of what the MI would look like,
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination value.
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI.
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p>This method (in contrast to {@link #computeSignificance(int)})
	 * allows the user to specify how to construct the surrogates,
	 * such that repeatable results may be obtained.</p>
	 * 
	 * <p>We provide a simple implementation which would be suitable for
	 *  any child class, though the child class may prefer to make its
	 *  own implementation to make class-specific optimisations.
	 *  Child classes must implement {@link java.lang.Cloneable}
	 *  for this method to be callable for them, and indeed implement
	 *  a <code>clone()</code> method in a way that protects their structure
	 *  from alteration by surrogate data being supplied to it.</p>
	 *  
	 * <p>We permute the source variable against the destination
	 *  to be consistent with the description in {@link ChannelCalculator#computeSignificance(int[][])}
	 *  (though in theory this doesn't matter for this function call).
	 * </p>
	 * 
	 * @param newOrderings a specification of how to shuffle the next values
	 *  to create the surrogates to generate the distribution with. The first
	 *  index is the permutation number (i.e. newOrderings.length is the number
	 *  of surrogate samples we use to bootstrap to generate the distribution here.)
	 *  Each array newOrderings[i] should be an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 * @return the distribution of channel measure scores under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception where the length of each permutation in newOrderings
	 *   is not equal to the number N samples that were previously supplied.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception {
		
		int numSurrogatesToCheck = newOrderings.length;
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		
		double[] surrogateMeasurements = new double[numSurrogatesToCheck];
		
		// Now compute the MI for each set of shuffled data:
		for (int i = 0; i < numSurrogatesToCheck; i++) {
			// Compute a new surrogate MI
			surrogateMeasurements[i] = computeAverageLocalOfObservations(newOrderings[i]);
			if (debug){
				System.out.println("New MI was " + surrogateMeasurements[i]);
			}
		}
		
		return new EmpiricalMeasurementDistribution(surrogateMeasurements, lastAverage);
	}

	/**
	 * <p>Compute the mutual information if the first (source) variable were
	 *  ordered as per the ordering specified in newOrdering.</p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of a shuffled source series here.</p>
	 * 
	 * <p>This method is primarily intended for use in {@link #computeSignificance(int[][])}
	 * however has been made public in case users wish to access it.
	 * </p>
	 * 
	 * <p>We provide a simple implementation which would be suitable for
	 *  any child class, though the child class may prefer to make its
	 *  own implementation to make class-specific optimisations.
	 *  Child classes must implement {@link java.lang.Cloneable}
	 *  for this method to be callable for them, and indeed implement
	 *  the <code>clone()</code> method in a way that protects their structure
	 *  from alteration by surrogate data being supplied to it.</p>
	 *  
	 * <p>Child implementations must only over-write lastAverage
	 * if <code>newOrdering</code> is null, which indicates that the original
	 * ordering should be used (making this equivalent to a call to
	 * {@link #computeAverageLocalOfObservations()}).</p>
	 *  
	 * @param newOrdering a specification of how to shuffle the source values
	 *  to create a surrogate source time series.
	 *  It is an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 * @return the surrogate MI score if the source values were shuffled as specified.
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[] newOrdering)
			throws Exception {
		
		if (newOrdering == null) {
			return computeAverageLocalOfObservations();
		}
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays)
		MutualInfoMultiVariateCommon miSurrogateCalculator =
				(MutualInfoMultiVariateCommon) this.clone();
		// Turn off normalisation and adding noise here since the data will already have been normalised 
		//  and noise added with the first run of the calculator. Normalising again can cause complication if
		//  the original data had no standard deviation (normalising again now would inflate 
		//  the small added noise values to the standard scale, and bring a range of CMI values
		//  instead of just the zeros that we should otherwise get).
		miSurrogateCalculator.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE, "false");
		miSurrogateCalculator.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0");

		// Generate a new re-ordered source data
		double[][] shuffledSourceData =
				MatrixUtils.extractSelectedTimePointsReusingArrays(
					sourceObservations, newOrdering);
		// Perform new initialisations
		miSurrogateCalculator.initialise(dimensionsSource, dimensionsDest);
		// Set new observations
		miSurrogateCalculator.setObservations(shuffledSourceData, destObservations);
		// Compute the MI
		return miSurrogateCalculator.computeAverageLocalOfObservations();
	}

	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations)
					throws Exception {
		if ((dimensionsDest != 1) || (dimensionsSource != 1)) {
			throw new Exception("The number of source and dest dimensions (having been initialised to " +
					dimensionsSource + " and " + dimensionsDest + ") can only be 1 when " +
					"the univariate addObservations(double[],double[]) and " + 
					"setObservations(double[],double[]) methods are called");
		}
		return computeLocalUsingPreviousObservations(
				MatrixUtils.reshape(newSourceObservations, newSourceObservations.length, 1),
				MatrixUtils.reshape(newDestObservations, newDestObservations.length, 1));
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	/**
	 * Return the MI last calculated in a call to {@link #computeAverageLocalOfObservations()}
	 * or {@link #computeLocalOfPreviousObservations()} after the previous
	 * {@link #initialise()} call.
	 * 
	 * @return the last computed average mutual information
	 */
	public double getLastAverage() {
		return lastAverage;
	}

	/**
	 * Get the number of samples to be used for the PDFs here 
	 * which have been supplied by calls to
	 * {@link #setObservations(double[][], double[][])},
	 * {@link #addObservations(double[][], double[][])}
	 * etc.
	 * 
	 * <p>Note that the number of samples may not be equal to the length of time-series
	 * supplied (i.e. where a {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}
	 * is set).
	 * </p>
	 * 
	 * @throws Exception if child class computes MI without explicit observations
	 *  (e.g. {@link infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian})
	 */
	public int getNumObservations() throws Exception {
		return totalObservations;
	}

	public boolean getAddedMoreThanOneObservationSet() {
		return addedMoreThanOneObservationSet;
	}

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid series of both source and destinations.
	 * 
	 * <p>Made public so it can be used if one wants to compute the number of
	 *  observations prior to setting the observations.</p>
	 * 
	 * @param sourceValid a time series (with indices the same as observations)
	 *  indicating whether the entry in observations at that index is valid for the source; 
	 * @param destValid as described for <code>sourceValid</code>
	 * @return a vector for start and end time pairs of valid series
	 *  of observations.
	 */
	public Vector<int[]> computeStartAndEndTimePairs(boolean[] sourceValid, boolean[] destValid) {
		// Scan along the data avoiding invalid values
		int startTime = 0;
		int endTime = 0;
		boolean lookingForStart = true;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();
		for (int t = 0; t < destValid.length; t++) {
			if (lookingForStart) {
				// Precondition: startTime holds a candidate start time
				//  (source value is at startTime == t - timeDiff)
				if (destValid[t] && sourceValid[t - timeDiff]) {
					// This point is OK at the source and destination
					// Set a candidate endTime
					endTime = t;
					lookingForStart = false;
					if (t == destValid.length - 1) {
						// we need to terminate now
						int[] timePair = new int[2];
						timePair[0] = startTime;
						timePair[1] = endTime;
						startAndEndTimePairs.add(timePair);
						// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
					}
				} else {
					// We need to keep looking.
					// Move the potential start time to the next point
					startTime++;
				}
			} else {
				// Precondition: startTime holds the start time for this set, 
				//  endTime holds a candidate end time
				// Check if we can include the current time step
				boolean terminateSequence = false;
				if (destValid[t] && sourceValid[t - timeDiff]) {
					// We can extend
					endTime = t;
				} else {
					terminateSequence = true;
				}
				if (t == destValid.length - 1) {
					// we need to terminate the sequence anyway
					terminateSequence = true;
				}
				if (terminateSequence) {
					// This section is done
					int[] timePair = new int[2];
					timePair[0] = startTime;
					timePair[1] = endTime;
					startAndEndTimePairs.add(timePair);
					// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
					lookingForStart = true;
					startTime = t + 1;
				}
			}
		}
		return startAndEndTimePairs;
	}	

	@Override
	public int[] getObservationSetIndices() {
		return observationSetIndices;
	}

	@Override
	public int[] getObservationTimePoints() {
		return observationTimePoints;
	}
}
