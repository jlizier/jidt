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

package infodynamics.measures.continuous.kraskov;

import java.util.Random;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.MutualInfoMultiVariateCommon;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential mutual information of two given multivariate sets of
 *  observations (implementing {@link MutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below).
 *  This is an abstract class to gather common functionality between the two
 *  algorithms defined by Kraskov et al.
 *  Two child classes {@link MutualInfoCalculatorMultiVariateKraskov1} and
 *  {@link MutualInfoCalculatorMultiVariateKraskov2} then
 *  actually implement the two algorithms in the Kraskov et al. paper</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link MutualInfoCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>For constructors see the child classes.</li>
 *  <li>Further properties are defined in {@link #setProperty(String, String)}.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 *  
 * <p>
 * TODO Add fast nearest neighbour searches to the child classes
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class MutualInfoCalculatorMultiVariateKraskov
	extends MutualInfoMultiVariateCommon
	implements MutualInfoCalculatorMultiVariate {

	/**
	 * we compute distances to the kth nearest neighbour
	 */
	protected int k = 4;
	
	/**
	 * Calculator for the norm between data points
	 */
	protected EuclideanUtils normCalculator;
	/**
	 * Cache for the norms between x (source) points
	 */
	protected double[][] xNorms;
	/**
	 * Cache for the norms between x (dest) points
	 */
	protected double[][] yNorms;
	/**
	 * Whether we cache the norms each time (making reordering very quick).
	 *  (Should only be set to false for testing)
	 */
	public static boolean tryKeepAllPairsNorms = true;
	/**
	 * An upper limit on the number of samples for which
	 * we will cache the norms between data points.
	 */
	public static int MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM = 2000;
	
	/**
	 * Property name for the number of K nearest neighbours used in
	 * the KSG algorithm in the full joint space (default 4).
	 */
	public final static String PROP_K = "k";
	/**
	 * Property name for what type of norm to use between data points
	 *  for each marginal variable -- Options are defined by 
	 *  {@link EuclideanUtils#setNormToUse(String)} and the
	 *  default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 */
	public final static String PROP_NORM_TYPE = "NORM_TYPE";
	/**
	 * Property name for whether to normalise the incoming data to 
	 * mean 0, standard deviation 1 (default true)
	 */
	public static final String PROP_NORMALISE = "NORMALISE";
	/**
	 * Property name for an amount of random Gaussian noise to be
	 *  added to the data (default is 0).
	 */
	public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";
	
	/**
	 * Whether to normalise the incoming data 
	 */
	protected boolean normalise = true;
	/**
	 * Whether to add an amount of random noise to the incoming data
	 */
	protected boolean addNoise = false;
	/**
	 * Amount of random Gaussian noise to add to the incoming data
	 */
	protected double noiseLevel = 0.0;

	/**
	 * Construct an instance of the KSG MI calculator
	 */
	public MutualInfoCalculatorMultiVariateKraskov() {
		super();
		normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
	}

	public void initialise(int sourceDimensions, int destDimensions) {
		super.initialise(sourceDimensions, destDimensions);
		xNorms = null;
		yNorms = null;
	}

	/**
	 * Sets properties for the KSG MI calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 *  <li>{@link #PROP_K} -- number of k nearest neighbours to use in joint kernel space
	 *      in the KSG algorithm (default is 4).</li>
	 * 	<li>{@link #PROP_NORM_TYPE}</li> -- normalization type to apply to 
	 * 		working out the norms between the points in each marginal space.
	 * 		Options are defined by {@link EuclideanUtils#setNormToUse(String)} -
	 * 		default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 *  <li>{@link #PROP_NORMALISE} -- whether to normalise the incoming individual
	 *      variables to mean 0 and standard deviation 1 (true by default)</li>
	 *  <li>{@link #PROP_ADD_NOISE} -- an amount of random noise to add to
	 *       each variable, to avoid having neighbourhoods with artificially
	 *       large counts. The amount is added in before any normalisation. 
	 *       (Recommended by Kraskov. MILCA uses 1e-8; but adds in
	 *       a random amount of noise in [0,noiseLevel) ). Default 0.</li>
	 *  <li>any valid properties for {@link MutualInfoMultiVariateCommon#setProperty(String, String)}.</li>
	 * </ul>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normCalculator.setNormToUse(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
			addNoise = true;
			noiseLevel = Double.parseDouble(propertyValue);
		} else {
			// No property was set here
			propertySet = false;
			// try the superclass:
			super.setProperty(propertyName, propertyValue);
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.MutualInfoMultiVariateCommon#finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		// Allow the parent to generate the data for us first
		super.finaliseAddObservations();
		
		if (addNoise) {
			Random random = new Random();
			// Add Gaussian noise of std dev noiseLevel to the data
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
		
		// Then normalise the data if required
		if (normalise) {
			// Take a copy since we're going to normalise it
			sourceObservations = MatrixUtils.normaliseIntoNewArray(sourceObservations);
			destObservations = MatrixUtils.normaliseIntoNewArray(destObservations);
		}
	}

	/**
	 * Utility function to compute the norms between each pair of points in each marginal time series
	 *
	 */
	protected void computeNorms() {
		int N = sourceObservations.length; // number of observations
		
		xNorms = new double[N][N];
		yNorms = new double[N][N];
		for (int t = 0; t < N; t++) {
			// Compute the norms from t to all other time points
			double[][] xyNormsForT = normCalculator.computeNorms(sourceObservations, destObservations, t);
			for (int t2 = 0; t2 < N; t2++) {
				xNorms[t][t2] = xyNormsForT[t2][0];
				yNorms[t][t2] = xyNormsForT[t2][1];
			}
		}
	}
	
	/**
	 * Compute the average MI from the previously supplied observations.
	 * 
	 * @return the average MI in nats (not bits!)
	 */
	public abstract double computeAverageLocalOfObservations() throws Exception;

	/**
	 * @return the MI under the new ordering, in nats (not bits!).
	 *  Returns NaN if any of the determinants are zero
	 *  (because this will make the denominator of the log 0).
	 */
	public abstract double computeAverageLocalOfObservations(int[] reordering) throws Exception;

	/**
	 * <p>Computes the local values of the MI,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>If the samples were supplied via a single call such as
	 * {@link #setObservations(double[])},
	 * then the return value is a single time-series of local
	 * channel measure values corresponding to these samples.</p>
	 * 
	 * <p>Otherwise where disjoint time-series observations were supplied using several 
	 *  calls such as {@link addObservations(double[])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single "time-series" return array.</p>
	 * 
	 * @return the "time-series" of local MIs in bits
	 * @throws Exception
	 */
	public abstract double[] computeLocalOfPreviousObservations() throws Exception;

	/**
	 * This method, specified in {@link MutualInfoCalculatorMultiVariate}
	 * is not implemented yet here.
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] states1, double[][] states2) throws Exception {
		// TODO If implemented, will need to incorporate any time difference here.
		// Will also need to handle normalisation of the incoming data
		//  appropriately
		throw new Exception("Local method not implemented yet");
	}
	
	/**
	 * Utility function used for debugging, printing digamma constants
	 * 
	 * @param N
	 * @return
	 * @throws Exception
	 */
	public abstract String printConstants(int N) throws Exception ;	
}
