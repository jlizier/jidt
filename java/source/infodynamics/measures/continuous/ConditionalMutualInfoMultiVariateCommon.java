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

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;
import java.util.Vector;

/**
 * Implements {@link ConditionalMutualInfoCalculatorMultiVariate}
 * to provide a base
 * class with common functionality for child class implementations of
 * {@link ConditionalMutualInfoCalculatorMultiVariate}
 * via various estimators. 
 * 
 * <p>These various estimators include: e.g. box-kernel estimation, KSG estimators, etc
 * (see the child classes linked above).
 * </p>
 * 
 * <p>Usage is as outlined in {@link ConditionalMutualInfoCalculatorMultiVariate}.</p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class ConditionalMutualInfoMultiVariateCommon implements
		ConditionalMutualInfoCalculatorMultiVariate {

	/**
	 * Number of dimenions for variable 1
	 */
	protected int dimensionsVar1 = 1;
	/**
	 * Number of dimenions for variable 2
	 */
	protected int dimensionsVar2 = 1;
	/**
	 * Number of dimenions for the conditional variable
	 */
	protected int dimensionsCond = 1;
	
	/**
	 * The set of observations for var1, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link #addObservations(double[][], double[][], double[][])} functions.
	 */
	protected double[][] var1Observations;

	/**
	 * The set of observations for var2, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link #addObservations(double[][], double[][], double[][])} functions.
	 */
	protected double[][] var2Observations;

	/**
	 * The set of observations for the conditional, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link #addObservations(double[][], double[][], double[][])} functions.
	 */
	protected double[][] condObservations;

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
	 * Store the last computed average conditional MI
	 */
	protected double lastAverage;

	/**
	 * Track whether we've computed the average for the supplied
	 *  observations yet
	 */
	protected boolean condMiComputed;

	/**
	 * Whether to report debug messages or not
	 */
	protected boolean debug;

	/**
	 * Storage for var1 observations supplied via
	 * {@link #addObservations(double[][], double[][], double[][])} etc
	 */
	protected Vector<double[][]> vectorOfVar1Observations;
	/**
	 * Storage for var2 observations supplied via
	 * {@link #addObservations(double[][], double[][], double[][])} etc
	 */
	protected Vector<double[][]> vectorOfVar2Observations;
	/**
	 * Storage for conditional variable observations supplied via
	 * {@link #addObservations(double[][], double[][], double[][])} etc
	 */
	protected Vector<double[][]> vectorOfCondObservations;
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
	 * Whether the user has added more than one disjoint observation set
	 * via {@link #addObservations(double[][], double[][], double[][])} etc
	 */
	protected boolean addedMoreThanOneObservationSet;

	/**
	 * Member to track whether PROP_NORMALISE has been set
	 */
	protected boolean normalise = true;
	/**
	 * Whether to add an amount of random noise to the incoming data
	 */
	protected boolean addNoise = true;
	/**
	 * Amount of random Gaussian noise to add to the incoming data.
	 * 0 by default except for KSG estimators (where it is recommended
	 *  and 1e-8 is used to match MILCA toolkit)
	 */
	protected double noiseLevel = (double) 0;
	/**
	 * Has the user set a seed for the random noise
	 */
	protected boolean noiseSeedSet = false;
	/**
	 * Seed that the user set for the random noise
	 */
	protected long noiseSeed = 0;

	/**
	 * Cache for the means of each dimension in variable 1, in case we need to normalise 
	 *  new observations later
	 */
	protected double[] var1Means = null;
	/**
	 * Cache for the standard deviations of each dimension in variable 1, in case we need to normalise 
	 *  new observations later
	 */
	protected double[] var1Stds = null;
	/**
	 * Cache for the means of each dimension in variable 2, in case we need to normalise 
	 *  new observations later
	 */
	protected double[] var2Means = null;
	/**
	 * Cache for the standard deviations of each dimension in variable 2, in case we need to normalise 
	 *  new observations later
	 */
	protected double[] var2Stds = null;
	/**
	 * Cache for the means of each dimension in conditional variable, in case we need to normalise 
	 *  new observations later
	 */
	protected double[] condMeans = null;
	/**
	 * Cache for the standard deviations of each dimension in conditional variable, in case we need to normalise 
	 *  new observations later
	 */
	protected double[] condStds = null;

	@Override
	public void initialise() {
		initialise(dimensionsVar1, dimensionsVar2, dimensionsCond);
	}
	
	@Override
	public void initialise(int var1Dimensions, int var2Dimensions, int condDimensions) {
		dimensionsVar1 = var1Dimensions;
		dimensionsVar2 = var2Dimensions;
		dimensionsCond = condDimensions;
		lastAverage = 0.0;
		totalObservations = 0;
		condMiComputed = false;
		var1Observations = null;
		var2Observations = null;
		condObservations = null;
		observationSetIndices = null;
		observationTimePoints = null;
		observationSetIndex = 0;
		vectorOfObservationSetIndices = null;
		vectorOfObservationStartTimePoints = null;
		addedMoreThanOneObservationSet = false;
		var1Means = null;
		var1Stds = null;
		var2Means = null;
		var2Stds = null;
		condMeans = null;
		condStds = null;
	}

	/**
	 * Set properties for the calculator.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 *   <li>{@link #PROP_NORMALISE} - whether to normalise the individual
	 *      variables to mean 0, standard deviation 1
	 *      (true by default, except for child class
	 *      {@link infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian}
	 *      for which this property is false and cannot be altered)</li>
	 *   <li>{@link #PROP_ADD_NOISE} -- a standard deviation for an amount of
	 *  	random Gaussian noise to add to
	 *      each variable, to avoid having neighbourhoods with artificially
	 *      large counts. (We also accept "false" to indicate "0".)
	 *      The amount is added in after any normalisation,
	 *      so can be considered as a number of standard deviations of the data.
	 *      (Default is 0, except for KSG estimators where it is recommended by Kraskov
	 *      and so they use 1e-8 to match the MILCA toolkit, although that adds in
	 *      a random amount of noise in [0,noiseLevel) ).</li>
	 *  <li>{@link #PROP_NOISE_SEED} -- a long value seed for the random noise generator or
	 *      the string {@link ConditionalMutualInfoCalculatorMultiVariate#NOISE_NO_SEED_VALUE} for no seed (default)</li>
	 * </ul>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) {

		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
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
	public String getProperty(String propertyName) {
		if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			return Boolean.toString(normalise);
		} else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
			return Double.toString(noiseLevel);
	    } else if (propertyName.equalsIgnoreCase(PROP_NOISE_SEED)) {
	    	if (noiseSeedSet) {
	    		return Long.toString(noiseSeed);
	    	} else {
	    		return NOISE_NO_SEED_VALUE;
	    	}
		} else {
			// No property matches for this class
			return null;
		}
	}

	@Override
	public void setObservations(double[][] var1, double[][] var2, double[][] cond) throws Exception {
		startAddObservations();
		addObservations(var1, var2, cond);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	/**
	 * A non-overloaded method signature for setObservations with 2D arguments, as there have been
	 *  some problems calling overloaded versions of setObservations from python jpype.
	 *  Resolved if one follows the AutoAnalyser generated code, but left for back compatibility.
	 * 
	 * @param var1
	 * @param var2
	 * @param cond
	 * @throws Exception
	 */
	public void setObservations2D(double[][] var1, double[][] var, double[][] cond) throws Exception {
		setObservations(var1, var, cond);
	}

	@Override
	public void setObservations(double[] var1, double[] var2, double[] cond) throws Exception {
		startAddObservations();
		addObservations(var1, var2, cond);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	@Override
	public void setObservations(double[][] var1, double[] var2, double[] cond) throws Exception {
		startAddObservations();
		addObservations(var1, var2, cond);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	@Override
	public void setObservations(double[] var1, double[][] var2, double[] cond) throws Exception {
		startAddObservations();
		addObservations(var1, var2, cond);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	@Override
	public void setObservations(double[] var1, double[] var2, double[][] cond) throws Exception {
		startAddObservations();
		addObservations(var1, var2, cond);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	@Override
	public void setObservations(double[] var1, double[][] var2, double[][] cond) throws Exception {
		startAddObservations();
		addObservations(var1, var2, cond);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	@Override
	public void setObservations(double[][] var1, double[] var2, double[][] cond) throws Exception {
		startAddObservations();
		addObservations(var1, var2, cond);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	@Override
	public void setObservations(double[][] var1, double[][] var2, double[] cond) throws Exception {
		startAddObservations();
		addObservations(var1, var2, cond);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	/**
	 * A non-overloaded method signature for setObservations with 1D arguments, as there have been
	 *  some problems calling overloaded versions of setObservations from python jpype.
	 *  Resolved if one follows the AutoAnalyser generated code, but left for back compatibility.
	 * 
	 * @param var1
	 * @param var2
	 * @param cond
	 * @throws Exception
	 */
	public void setObservations1D(double[] var1, double[] var, double[] cond) throws Exception {
		setObservations(var1, var, cond);
	}
	
	@Override
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[] var1Valid, boolean[] var2Valid,
			boolean[] condValid) throws Exception {
		
		startAddObservations();
		addObservations(var1, var2, cond, var1Valid, var2Valid, condValid);
		finaliseAddObservations();
	}

	@Override
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[][] var1Valid, boolean[][] var2Valid,
			boolean[][] condValid) throws Exception {

		startAddObservations();
		addObservations(var1, var2, cond, var1Valid, var2Valid, condValid);
		finaliseAddObservations();
	}

	@Override
	public void startAddObservations() {
		vectorOfVar1Observations = new Vector<double[][]>();
		vectorOfVar2Observations = new Vector<double[][]>();
		vectorOfCondObservations = new Vector<double[][]>();
		vectorOfObservationSetIndices = new Vector<Integer>();
		vectorOfObservationStartTimePoints = new Vector<Integer>();
	}
	
	@Override
	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond) throws Exception {
		// Use the current observationSetIndex and increment for next use:
		addObservationsTrackObservationIDs(var1, var2, cond, observationSetIndex++, 0);
	}
	
	@Override
	public void addObservationsTrackObservationIDs(double[][] var1, double[][] var2,
			double[][] cond, int observationSetIndexToUse, int startTimeIndex) throws Exception {
		if (vectorOfVar1Observations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if ((var1.length != var2.length) ||
				((dimensionsCond != 0) && (var1.length != cond.length))) {
			throw new Exception(String.format("Observation vector lengths (%d, %d and %d) must match!",
					var1.length, var2.length, cond.length));
		}
		if (var1[0].length != dimensionsVar1) {
			throw new Exception("Number of joint variables in var1 data " +
					"does not match the initialised value");
		}
		if (var2[0].length != dimensionsVar2) {
			throw new Exception("Number of joint variables in var2 data " +
					"does not match the initialised value");
		}
		if ((dimensionsCond != 0) && (cond[0].length != dimensionsCond)) {
			throw new Exception("Number of joint variables in cond data " +
					"does not match the initialised value");
		}
		vectorOfVar1Observations.add(var1);
		vectorOfVar2Observations.add(var2);
		vectorOfCondObservations.add(cond);
		vectorOfObservationSetIndices.add(observationSetIndexToUse);
		vectorOfObservationStartTimePoints.add(startTimeIndex);
		if (vectorOfVar1Observations.size() > 1) {
			addedMoreThanOneObservationSet = true;
		}
	}

	/**
	 * A non-overloaded method signature for addObservations with 2D arguments, as there have been
	 *  some problems calling overloaded versions of setObservations from jpype. 
	 * 
	 * @param var1
	 * @param var2
	 * @param cond
	 * @throws Exception
	 */
	public void addObservations2D(double[][] var1, double[][] var, double[][] cond) throws Exception {
		addObservations(var1, var, cond);
	}

	@Override
	public void addObservations(double[] var1, double[] var2,
			double[] cond) throws Exception {
		if ((dimensionsVar1 != 1) || (dimensionsVar2 != 1) ||
				((dimensionsCond != 1) && (dimensionsCond != 0))) {
			throw new Exception("The number of dimensions for each variable (having been initialised to " +
					dimensionsVar1 + ", " + dimensionsVar2 + " & " +
					dimensionsCond + ") can only be 1 (or 0 for conditional) when " +
					"the univariate addObservations(double[],double[],double[]) and " + 
					"setObservations(double[],double[],double[]) methods are called");
		}
		double[][] reshapedConditional = null;
		if (dimensionsCond == 1) {
			// This won't execute if dimensionsCond == 0
			reshapedConditional = MatrixUtils.reshape(cond, cond.length, 1);
		}
		addObservations(MatrixUtils.reshape(var1, var1.length, 1),
						MatrixUtils.reshape(var2, var2.length, 1),
						reshapedConditional);

	}
	
	@Override
	public void addObservations(double[][] var1, double[] var2,
			double[] cond) throws Exception {
		if ((dimensionsVar2 != 1) ||
				((dimensionsCond != 1) && (dimensionsCond != 0))) {
			throw new Exception("The number of dimensions for variables var2 and cond (having been initialised to " +
					dimensionsVar2 + " & " +
					dimensionsCond + ") can only be 1 (or 0 for conditional) when " +
					"the addObservations(double[][],double[],double[]) and " + 
					"setObservations(double[][],double[],double[]) methods are called");
		}
		double[][] reshapedConditional = null;
		if (dimensionsCond == 1) {
			// This won't execute if dimensionsCond == 0
			reshapedConditional = MatrixUtils.reshape(cond, cond.length, 1);
		}
		addObservations(var1,
						MatrixUtils.reshape(var2, var2.length, 1),
						reshapedConditional);
	}

	@Override
	public void addObservations(double[] var1, double[][] var2,
			double[] cond) throws Exception {
		if ((dimensionsVar1 != 1) ||
				((dimensionsCond != 1) && (dimensionsCond != 0))) {
			throw new Exception("The number of dimensions for variables var1 and cond (having been initialised to " +
					dimensionsVar1 + " & " +
					dimensionsCond + ") can only be 1 (or 0 for conditional) when " +
					"the addObservations(double[],double[][],double[]) and " + 
					"setObservations(double[],double[][],double[]) methods are called");
		}
		double[][] reshapedConditional = null;
		if (dimensionsCond == 1) {
			// This won't execute if dimensionsCond == 0
			reshapedConditional = MatrixUtils.reshape(cond, cond.length, 1);
		}
		addObservations(MatrixUtils.reshape(var1, var1.length, 1),
						var2,
						reshapedConditional);
	}

	@Override
	public void addObservations(double[] var1, double[] var2,
			double[][] cond) throws Exception {
		if ((dimensionsVar1 != 1) || (dimensionsVar2 != 1)) {
			throw new Exception("The number of dimensions for variables var1 and var2 (having been initialised to " +
					dimensionsVar1 + " & " +
					dimensionsVar2 + ") can only be 1 when " +
					"the addObservations(double[],double[],double[][]) and " + 
					"setObservations(double[],double[],double[][]) methods are called");
		}
		addObservations(MatrixUtils.reshape(var1, var1.length, 1),
						MatrixUtils.reshape(var2, var2.length, 1),
						cond);
	}

	@Override
	public void addObservations(double[] var1, double[][] var2,
			double[][] cond) throws Exception {
		if (dimensionsVar1 != 1) {
			throw new Exception("The number of dimensions for variable var1 (having been initialised to " +
					dimensionsVar1 + ") can only be 1 when " +
					"the addObservations(double[],double[][],double[][]) and " + 
					"setObservations(double[],double[][],double[][]) methods are called");
		}
		addObservations(MatrixUtils.reshape(var1, var1.length, 1),
						var2,
						cond);
	}

	@Override
	public void addObservations(double[][] var1, double[] var2,
			double[][] cond) throws Exception {
		if (dimensionsVar2 != 1) {
			throw new Exception("The number of dimensions for variable var2 (having been initialised to " +
					dimensionsVar2 + ") can only be 1 when " +
					"the addObservations(double[][],double[],double[][]) and " + 
					"setObservations(double[][],double[],double[][]) methods are called");
		}
		addObservations(var1,
						MatrixUtils.reshape(var2, var2.length, 1),
						cond);
	}

	@Override
	public void addObservations(double[][] var1, double[][] var2,
			double[] cond) throws Exception {
		if ((dimensionsCond != 1) && (dimensionsCond != 0)) {
			throw new Exception("The number of dimensions for variable cond (having been initialised to " +
					dimensionsCond + ") can only be 1 or 0 when " +
					"the addObservations(double[][],double[][],double[]) and " + 
					"setObservations(double[][],double[][],double[]) methods are called");
		}
		double[][] reshapedConditional = null;
		if (dimensionsCond == 1) {
			// This won't execute if dimensionsCond == 0
			reshapedConditional = MatrixUtils.reshape(cond, cond.length, 1);
		}
		addObservations(var1,
						var2,
						reshapedConditional);
	}

	/**
	 * A non-overloaded method signature for addObservations with 1D arguments, as there have been
	 *  some problems calling overloaded versions of setObservations from jpype. 
	 * 
	 * @param var1
	 * @param var2
	 * @param cond
	 * @throws Exception
	 */
	public void addObservations1D(double[] var1, double[] var, double[] cond) throws Exception {
		addObservations(var1, var, cond);
	}
	
	@Override
	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond,
			int startTime, int numTimeSteps) throws Exception {
		addObservations(var1, var2, cond, startTime, numTimeSteps, observationSetIndex++);
	}
	
	protected void addObservations(double[][] var1, double[][] var2,
			double[][] cond,
			int startTime, int numTimeSteps, int observationSetIndexToUse) throws Exception {
		if (vectorOfVar1Observations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		double[][] var1ToAdd = new double[numTimeSteps][];
		System.arraycopy(var1, startTime, var1ToAdd, 0, numTimeSteps);
		double[][] var2ToAdd = new double[numTimeSteps][];
		System.arraycopy(var2, startTime, var2ToAdd, 0, numTimeSteps);
		double[][] condToAdd = null;
		if (dimensionsCond != 0) {
			condToAdd = new double[numTimeSteps][];
			System.arraycopy(cond, startTime, condToAdd, 0, numTimeSteps);
		}
		addObservationsTrackObservationIDs(var1ToAdd, var2ToAdd, condToAdd, observationSetIndexToUse, startTime);
	}

	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[] var1Valid, boolean[] var2Valid,
			boolean[] condValid) throws Exception {

		Vector<int[]> startAndEndTimePairs =
				computeStartAndEndTimePairs(var1Valid, var2Valid, condValid);
		
		// We've found the set of start and end times for this pair
		startAddObservations();
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservations(var1, var2, cond, startTime, endTime - startTime + 1, observationSetIndex);
		}
		observationSetIndex++;
		finaliseAddObservations();
	}
	
	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[][] var1Valid, boolean[][] var2Valid,
			boolean[][] condValid) throws Exception {
		boolean[] allVar1Valid = MatrixUtils.andRows(var1Valid);
		boolean[] allVar2Valid = MatrixUtils.andRows(var2Valid);
		boolean[] allCondValid = MatrixUtils.andRows(condValid);
		addObservations(var1, var2, cond, allVar1Valid, allVar2Valid, allCondValid);
	}
	
	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
	 * 
	 * This default implementation simply puts all of the observations into
	 *  the {@link #var1Observations}, {@link #var2Observations} 
	 *  and {@link #condObservations} arrays.
	 * Usually child implementations will override this, call this implementation
	 *  to perform the common processing, then perform their own processing.
	 *  
	 * @throws Exception Allow child classes to throw an exception if there
	 *  is an issue detected specific to that calculator.
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[][] var2 : vectorOfVar2Observations) {
			totalObservations += var2.length;
		}
		var1Observations = new double[totalObservations][dimensionsVar1];
		var2Observations = new double[totalObservations][dimensionsVar2];
		condObservations = new double[totalObservations][dimensionsCond];
		observationSetIndices = new int[totalObservations];
		observationTimePoints = new int[totalObservations];
		
		int startObservation = 0;
		Iterator<double[][]> iteratorVar2 = vectorOfVar2Observations.iterator();
		Iterator<double[][]> iteratorCond = vectorOfCondObservations.iterator();
		Iterator<Integer> iteratorObsSetIndices = vectorOfObservationSetIndices.iterator();
		Iterator<Integer> iteratorObsStartTimePoints = vectorOfObservationStartTimePoints.iterator();
		for (double[][] var1 : vectorOfVar1Observations) {
			double[][] var2 = iteratorVar2.next();
			double[][] cond = iteratorCond.next();
			// Copy the data from these given observations into our master 
			//  array
			MatrixUtils.arrayCopy(var1, 0, 0,
					var1Observations, startObservation, 0,
					var1.length, dimensionsVar1);
			MatrixUtils.arrayCopy(var2, 0, 0,
					var2Observations, startObservation, 0,
					var2.length, dimensionsVar2);
			if (dimensionsCond != 0) {
				MatrixUtils.arrayCopy(cond, 0, 0,
					condObservations, startObservation, 0,
					cond.length, dimensionsCond);
			} // else we can do nothing there
			// And update which observation set and time index each sample came from:
			Arrays.fill(observationSetIndices, startObservation, startObservation+var1.length, iteratorObsSetIndices.next());
			int firstTimeSampleId = iteratorObsStartTimePoints.next();
			for (int i = 0; i < var1.length; i++) {
				observationTimePoints[startObservation + i] = firstTimeSampleId + i;
			}
			startObservation += var2.length;
		}
		
		// Normalise the data if required
		var1Means = MatrixUtils.means(var1Observations);
		var1Stds = MatrixUtils.stdDevs(var1Observations, var1Means);
		var2Means = MatrixUtils.means(var2Observations);
		var2Stds = MatrixUtils.stdDevs(var2Observations, var2Means);
		if (dimensionsCond != 0) {
			condMeans = MatrixUtils.means(condObservations);
			condStds = MatrixUtils.stdDevs(condObservations, condMeans);
		}
		if (normalise) {
			normaliseData();
		}
		
		// We don't need to keep the vectors of observation sets anymore:
		vectorOfVar1Observations = null;
		vectorOfVar2Observations = null;
		vectorOfCondObservations = null;

		// Add Gaussian noise of std dev noiseLevel to the data if required
		if (addNoise) {
			Random random = new Random();
	    	if (noiseSeedSet) {
	    		random.setSeed(noiseSeed);
	    	}
			for (int r = 0; r < var1Observations.length; r++) {
				for (int c = 0; c < dimensionsVar1; c++) {
					var1Observations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
				for (int c = 0; c < dimensionsVar2; c++) {
					var2Observations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
				// This next loop will only execute if dimensionsCond > 0
				for (int c = 0; c < dimensionsCond; c++) {
					condObservations[r][c] +=
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
		MatrixUtils.normalise(var1Observations, var1Means, var1Stds);
		MatrixUtils.normalise(var2Observations, var2Means, var2Stds);
		if (dimensionsCond != 0) {
			MatrixUtils.normalise(condObservations, condMeans, condStds);
		}
	}
	
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		return computeSignificance(1, numPermutationsToCheck);
	}
	
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// Use var1 length (all variables have same length) even though
		//  we may be randomising the other variable:
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(
				totalObservations, numPermutationsToCheck);
		return computeSignificance(variableToReorder, newOrderings);
	}

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		return computeSignificance(1, newOrderings);
	}
	
	/**
	 * <p>As described in
	 * {@link ConditionalMutualInfoCalculatorMultiVariate#computeSignificance(int, int[][])}
	 * </p>
	 * 
	 * <p>Here we provide a simple implementation which would be suitable for
	 *  any child class, though the child class may prefer to make its
	 *  own implementation to make class-specific optimisations.
	 *  Child classes must implement {@link java.lang.Cloneable}
	 *  for this method to be callable for them, and indeed implement
	 *  the clone() method in a way that protects their structure
	 *  from alteration by surrogate data being supplied to it.</p>
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		if (!condMiComputed) {
			computeAverageLocalOfObservations();
		}
		
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays - child classes should override this)
		ConditionalMutualInfoMultiVariateCommon miSurrogateCalculator =
				(ConditionalMutualInfoMultiVariateCommon) this.clone();
		// Turn off normalisation and adding noise here since the data will already have been normalised 
		//  and noise added with the first run of the calculator. Normalising again can cause complication if
		//  the original data had no standard deviation (normalising again now would inflate 
		//  the small added noise values to the standard scale, and bring a range of CMI values
		//  instead of just the zeros that we should otherwise get).
		miSurrogateCalculator.setProperty(PROP_NORMALISE, "false");
		miSurrogateCalculator.setProperty(PROP_ADD_NOISE, "0");
		
		double[] surrogateMeasurements = new double[numPermutationsToCheck];
		
		// Now compute the MI for each set of shuffled data:
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Generate a new re-ordered source data
			double[][] shuffledData = 
					MatrixUtils.extractSelectedTimePointsReusingArrays(
							(variableToReorder == 1) ? var1Observations : var2Observations,
							newOrderings[i]);
			// Perform new initialisations
			miSurrogateCalculator.initialise(
					dimensionsVar1, dimensionsVar2, dimensionsCond);
			// Set new observations
			if (variableToReorder == 1) {
				miSurrogateCalculator.setObservations(shuffledData,
						var2Observations, condObservations);
			} else {
				miSurrogateCalculator.setObservations(var1Observations,
						shuffledData, condObservations);
			}
			// Compute the MI
			surrogateMeasurements[i] = miSurrogateCalculator.computeAverageLocalOfObservations();
			if (debug){
				System.out.println("New MI was " + surrogateMeasurements[i]);
			}
		}
		
		return new EmpiricalMeasurementDistribution(surrogateMeasurements, lastAverage);
	}

	/**
	 * <p>As described in
	 * {@link ConditionalMutualInfoCalculatorMultiVariate#computeAverageLocalOfObservations(int, int[])}
	 * </p>
	 * 
	 * <p>We provide a simple implementation which would be suitable for
	 *  any child class, though the child class may prefer to make its
	 *  own implementation to make class-specific optimisations.
	 *  Child classes must implement {@link java.lang.Cloneable}
	 *  for this method to be callable for them, and indeed implement
	 *  the clone() method in a way that protects their structure
	 *  from alteration by surrogate data being supplied to it.</p>
	 */
	@Override
	public double computeAverageLocalOfObservations(int variableToReorder, int[] newOrdering)
			throws Exception {
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays - child class should override this)
		ConditionalMutualInfoMultiVariateCommon miSurrogateCalculator =
				(ConditionalMutualInfoMultiVariateCommon) this.clone();
		
		// Generate a new re-ordered source data
		double[][] shuffledData =
				MatrixUtils.extractSelectedTimePointsReusingArrays(
					(variableToReorder == 1) ? var1Observations : var2Observations,
					newOrdering);
		// Perform new initialisations
		miSurrogateCalculator.initialise(
				dimensionsVar1, dimensionsVar2, dimensionsCond);
		// Set new observations
		if (variableToReorder == 1) {
			miSurrogateCalculator.setObservations(shuffledData,
					var2Observations, condObservations);
		} else {
			miSurrogateCalculator.setObservations(var1Observations,
					shuffledData, condObservations);
		}
		// Compute the MI
		return miSurrogateCalculator.computeAverageLocalOfObservations();
	}
	
	@Override
	public double[] computeLocalUsingPreviousObservations(double[] states1, double[] states2, double[][] condStates)
			throws Exception {
		if ((dimensionsVar1 != 1) || (dimensionsVar2 != 1)) {
			throw new Exception("The number of source and dest dimensions (having been initialised to " +
					dimensionsVar1 + " and " + dimensionsVar2 + ") can only be 1 when " +
					"the univariate computeLocalUsingPreviousObservations(double[],double[],double[][]) " + 
					"method is called");
		}
		return computeLocalUsingPreviousObservations(MatrixUtils.reshape(states1, states1.length, 1),
						MatrixUtils.reshape(states2, states2.length, 1),
						condStates);
	}

	@Override
	public double[] computeLocalUsingPreviousObservations(double[] states1, double[] states2, double[] condStates)
			throws Exception {
		if ((dimensionsVar1 != 1) || (dimensionsVar2 != 1) || (dimensionsCond != 1)) {
			throw new Exception("The number of source, dest and conditional dimensions (having been initialised to " +
					dimensionsVar1 + ", " + dimensionsVar2 + " and " + dimensionsCond + ") can only be 1 when " +
					"the univariate computeLocalUsingPreviousObservations(double[],double[],double[]) " + 
					"method is called");
		}
		return computeLocalUsingPreviousObservations(MatrixUtils.reshape(states1, states1.length, 1),
						MatrixUtils.reshape(states2, states2.length, 1),
						MatrixUtils.reshape(condStates, condStates.length, 1));
	}

	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return lastAverage;
	}

	public int getNumObservations() throws Exception {
		return totalObservations;
	}

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid series of all variables.
	 * 
	 * Made public so it can be used if one wants to compute the number of
	 *  observations prior to setting the observations.
	 * 
	 * @param var1Valid a series (indexed by observation number or time)
	 *  indicating whether the entry in observations at that index is valid for variable 1; 
	 * @param var2Valid as described for <code>var1Valid</code>
	 * @param condValid as described for <code>var1Valid</code>
	 * @return a vector for start and end time pairs of valid series
	 *  of observations.
	 */
	public Vector<int[]> computeStartAndEndTimePairs(
			boolean[] var1Valid, boolean[] var2Valid, boolean[] condValid) {
		// Scan along the data avoiding invalid values
		int startTime = 0;
		int endTime = 0;
		boolean lookingForStart = true;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();
		for (int t = 0; t < var2Valid.length; t++) {
			if (lookingForStart) {
				// Precondition: startTime holds a candidate start time
				//  (var1 value is at startTime == t)
				if (var1Valid[t] && var2Valid[t] && condValid[t]) {
					// This point is OK at the variables
					// Set a candidate endTime
					endTime = t;
					lookingForStart = false;
					if (t == var1Valid.length - 1) {
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
				if (var1Valid[t] && var2Valid[t] && condValid[t]) {
					// We can extend
					endTime = t;
				} else {
					terminateSequence = true;
				}
				if (t == var2Valid.length - 1) {
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

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate#getAddedMoreThanOneObservationSet()
	 */
	@Override
	public boolean getAddedMoreThanOneObservationSet() {
		return addedMoreThanOneObservationSet;
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
