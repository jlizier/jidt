/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2024, Joseph T. Lizier
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

import java.util.Random;
import java.util.Vector;

import infodynamics.utils.MatrixUtils;

/**
 * Implements {@link EntropyCalculatorMultiVariate} to provide a base
 * class with common functionality for child class implementations of
 * {@link EntropyCalculatorMultiVariate}
 * via various estimators. 
 * 
 * <p>These various estimators include: e.g. box-kernel estimation, Kozachenko, etc
 * (see the child classes linked above).
 * </p>
 * 
 * <p>Usage is as outlined in {@link EntropyCalculatorMultiVariate}.</p>
 * 
  * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class EntropyCalculatorMultiVariateCommon implements EntropyCalculatorMultiVariate {
	/**
	 * Whether we're in debug mode
	 */
	protected boolean debug = false;
	/**
	 * Number of observations supplied
	 */
	protected int totalObservations = 0;
	/**
	 * Number of joint variables/dimensions
	 */
	protected int dimensions = 1;
	/**
	 * Track whether we've computed the average for the supplied
	 *  observations yet
	 */
	protected boolean isComputed = false;
	/**
	 * Last computed average entropy
	 */
	protected double lastAverage = 0.0;
	/**
	 * Whether we normalise the incoming observations to mean 0,
	 * standard deviation 1.
	 */
	protected boolean normalise = true;
	/**
	 * Property for whether we normalise the incoming observations to mean 0,
	 * standard deviation 1.
	 */
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	/**
	 * Property name for an amount of random Gaussian noise to be
	 *  added to the data (default is 1e-8, matching the MILCA toolkit).
	 */
	public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";
	/**
	 * Whether to add an amount of random noise to the incoming data
	 */
	protected boolean addNoise = false;
	/**
	 * Amount of random Gaussian noise to add to the incoming data
	 */
	protected double noiseLevel = (double) 0.0;
	/**
	 * Storage for observations supplied via {@link #addObservations(double[][])}
	 * type calls
	 */
	protected Vector<double[][]> vectorOfObservations;
	/**
	 * The set of observations, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 */
	protected double[][] observations;

	/**
	 * Default constructor
	 */
	public EntropyCalculatorMultiVariateCommon() {
		// Nothing to do
	}
	
	@Override
	public void initialise() throws Exception {
		initialise(dimensions);
	}

	@Override
	public void initialise(int dimensions) {
		this.dimensions = dimensions;
		observations = null;
		totalObservations = 0;
		isComputed = false;
		lastAverage = 0;
	}

	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception {

		boolean propertySet = true;

		if (propertyName.equalsIgnoreCase(NUM_DIMENSIONS_PROP_NAME)) {
			dimensions = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
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
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}
	
	@Override
	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(NUM_DIMENSIONS_PROP_NAME)) {
			return Integer.toString(dimensions);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			return Boolean.toString(normalise);
	    } else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
	        return Double.toString(noiseLevel);
		} else {
			// No property was set, and no superclass to call:
			return null;
		}
	}

	@Override
	public void startAddObservations() {
		isComputed = false;
		totalObservations = 0;
		observations = null;
		vectorOfObservations = new Vector<double[][]>();
	}

	@Override
	public void finaliseAddObservations() throws Exception {

		observations = new double[totalObservations][dimensions];
		
		// Construct the joint vectors from the given observations
		int startObservation = 0;
		for (double[][] obs : vectorOfObservations) {
			if ((obs == null) || (obs.length < 1)) {
				continue;
			}
			// Dimension was checked earlier by addObservations.
			// Copy the data from these given observations into our master array
			MatrixUtils.arrayCopy(obs, 0, 0,
					observations, startObservation, 0,
					obs.length, dimensions);
			startObservation += obs.length;
		}

		// We don't need to keep the vector of observation sets anymore:
		vectorOfObservations = null;

		if (addNoise) {
			Random random = new Random();
			// Add Gaussian noise of std dev noiseLevel to the data
			for (int r = 0; r < totalObservations; r++) {
				for (int c = 0; c < dimensions; c++) {
					observations[r][c] += random.nextGaussian()*noiseLevel;
				}
			}
		}
	}

	@Override
	public void setObservations(double[][] observations) throws Exception {
		startAddObservations();
		addObservations(observations);
		finaliseAddObservations();
	}

	/*
	 * (non-Javadoc)
	 * @see infodynamics.measures.continuous.EntropyCalculator#setObservations(double[])
	 * 
	 * This method here to ensure we make compatibility with the 
	 *  EntropyCalculator interface; only valid if dimension > 1
	 */
	@Override
	public void setObservations(double[] observations) throws Exception {
		startAddObservations();
		addObservations(observations);
		finaliseAddObservations();
	}

	/**
	 * Shortcut method to join two sets of samples into a joint multivariate 
	 * 
	 * Each row of the data is an observation; each column of
	 *  the row is a new variable in the multivariate observation.
	 * This method signature allows the user to call setObservations for
	 *  joint time series without combining them into a single joint time
	 *  series (we do the combining for them).
	 * 
	 * @param data1 observations1 few variables in the joint data
	 * @param observations2 the other variables in the joint data
	 * @throws Exception When the length of the two arrays of observations do not match.
	 * @see #setObservations(double[][])
	 */
	@Deprecated
	public void setObservations(double[][] observations1, double[][] observations2)
			throws Exception {
		startAddObservations();
		addObservations(observations1, observations2);
		finaliseAddObservations();
	}

	@Override
	public void addObservations(double[][] observations) throws Exception {
		if (vectorOfObservations == null) {
		// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if ((observations != null) && (observations.length > 0)) {
			// Check the dimension:
			if (observations[0].length != dimensions) {
				throw new Exception(String.format("Cannot supply observations with dimension %d when expected dimension = %d",
						observations[0].length, dimensions));
			}
			totalObservations += observations.length;
		}
		// Add the observations whether empty or else valid:
		vectorOfObservations.add(observations);
	}

	/*
	 * (non-Javadoc)
	 * @see infodynamics.measures.continuous.EntropyCalculator#addObservations(double[])
	 * 
	 * This method here to ensure we make compatibility with the 
	 *  EntropyCalculator interface; only valid if dimension > 1
	 */
	@Override
	public void addObservations(double[] observations) throws Exception {
		if (dimensions != 1) {
			throw new Exception(String.format("Cannot set univariate observations when expected dimension = %d", dimensions));
		}
		double[][] observations2D = MatrixUtils.reshape(observations, observations.length, 1);
		addObservations(observations2D);
	}
	
	/**
	 * Shortcut method to join two sets of samples into a joint multivariate. 
	 * 
	 * Each row of the data is an observation; each column of
	 *  the row is a new variable in the multivariate observation.
	 * This method signature allows the user to call addObservations for
	 *  joint time series without combining them into a single joint time
	 *  series (we do the combining for them).
	 * 
	 * @param data1 first few variables in the joint data
	 * @param data2 the other variables in the joint data
	 * @throws Exception When the length of the two arrays of observations do not match.
	 * @see #addObservations(double[][])
	 */
	@Deprecated
	public void addObservations(double[][] data1,
			double[][] data2) throws Exception {
		int timeSteps = data1.length;
		if ((data1 == null) || (data2 == null)) {
			throw new Exception("Cannot have null data arguments");
		}
		if (data1.length != data2.length) {
			throw new Exception("Length of data1 (" + data1.length + ") is not equal to the length of data2 (" +
					data2.length + ")");
		}
		int data1Variables = data1[0].length;
		int data2Variables = data2[0].length;
		double[][] data = new double[timeSteps][data1Variables + data2Variables];
		for (int t = 0; t < timeSteps; t++) {
			System.arraycopy(data1[t], 0, data[t], 0, data1Variables);
			System.arraycopy(data2[t], 0, data[t], data1Variables, data2Variables);
		}
		// Now defer to the normal setObservations method
		addObservations(data);
	}

	@Override
	public int getNumObservations() throws Exception {
		return totalObservations;
	}

	@Override
	public double getLastAverage() {
		return lastAverage;
	}
	
	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

}
