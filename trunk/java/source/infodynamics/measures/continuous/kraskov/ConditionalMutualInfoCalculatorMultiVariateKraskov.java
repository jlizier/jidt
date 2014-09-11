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

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.utils.EuclideanUtils;

/**
 * <p>Computes the differential conditional mutual information of two multivariate
 *  <code>double[][]</code> sets of observations, conditioned on another
 *  (implementing {@link ConditionalMutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see references below).
 *  This is an abstract class, building on the common code base in
 *  {@link ConditionalMutualInfoMultiVariateCommon},
 *  to gather common functionality between the two
 *  algorithms defined by Kraskov et al.
 *  Two child classes {@link ConditionalMutualInfoCalculatorMultiVariateKraskov1} and
 *  {@link ConditionalMutualInfoCalculatorMultiVariateKraskov2} then
 *  actually implement the two KSG algorithms.</p>
 * 
 * <p>Crucially, the calculation is performed by examining
 * neighbours in the full joint space (as specified by Frenzel and Pompe)
 * rather than two MI calculators.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link ConditionalMutualInfoCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>For constructors see the child classes.</li>
 *  <li>Further properties are defined in {@link #setProperty(String, String)}.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 * 
 * Finally, note that {@link Cloneable} is implemented allowing clone()
 *  to produce only an automatic shallow copy, which is fine
 *  for the statistical significance calculation it is intended for
 *  (none of the array 
 *  data will be changed there).
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Frenzel and Pompe, <a href="http://dx.doi.org/10.1103/physrevlett.99.204101">
 * 	"Partial Mutual Information for Coupling Analysis of Multivariate Time Series"</a>,
 * 	Physical Review Letters, <b>99</b>, p. 204101+ (2007).</li>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class ConditionalMutualInfoCalculatorMultiVariateKraskov 
	extends ConditionalMutualInfoMultiVariateCommon
	implements Cloneable { // See comments on clonability above
	
	/**
	 * we compute distances to the kth neighbour in the joint space
	 */
	protected int k = 4;
		
	/**
	 * Calculator for the norm between data points
	 */
	protected EuclideanUtils normCalculator;
	/**
	 * Cache for the norms between x (variable 1) points
	 */
	protected double[][] xNorms;
	/**
	 * Cache for the norms between y (variable 2) points
	 */
	protected double[][] yNorms;
	/**
	 * Cache for the norms between z (conditionals) points
	 */
	protected double[][] zNorms;
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
	 * Property name for an amount of random Gaussian noise to be
	 *  added to the data (default is 0).
	 */
	public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";

	/**
	 * Whether to add an amount of random noise to the incoming data
	 */
	protected boolean addNoise = false;
	/**
	 * Amount of random Gaussian noise to add to the incoming data
	 */
	protected double noiseLevel = 0.0;

	/**
	 * Construct an instance of the KSG conditional MI calculator
	 */
	public ConditionalMutualInfoCalculatorMultiVariateKraskov() {
		super();
		normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
	}

	@Override
	public void initialise(int dimensions1, int dimensions2, int dimensionsCond) {
		super.initialise(dimensions1, dimensions2, dimensionsCond);
		
		xNorms = null;
		yNorms = null;
		zNorms = null;
	}

	/**
	 * Sets properties for the KSG conditional MI calculator.
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
	 *  <li>{@link #PROP_ADD_NOISE} -- a standard deviation for an amount of
	 *  	random Gaussian noise to add to
	 *      each variable, to avoid having neighbourhoods with artificially
	 *      large counts. The amount is added in after any normalisation,
	 *      so can be considered as a number of standard deviations of the data.
	 *      (Recommended by Kraskov. MILCA uses 1e-8; but adds in
	 *      a random amount of noise in [0,noiseLevel) ). Default 0.</li>
	 *  <li>any valid properties for {@link ConditionalMutualInfoMultiVariateCommon#setProperty(String, String)},
	 *     notably including {@link ConditionalMutualInfoMultiVariateCommon#PROP_NORMALISE}.</li>
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
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normCalculator.setNormToUse(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
			addNoise = true;
			noiseLevel = Double.parseDouble(propertyValue);
		} else {
			// Assume this is a property for the common parent class
			super.setProperty(propertyName, propertyValue);
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon#finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		// Allow the parent to generate the data for us first
		super.finaliseAddObservations();
		
		if (addNoise) {
			Random random = new Random();
			// Add Gaussian noise of std dev noiseLevel to the data
			for (int r = 0; r < var1Observations.length; r++) {
				for (int c = 0; c < dimensionsVar1; c++) {
					var1Observations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
				for (int c = 0; c < dimensionsVar2; c++) {
					var2Observations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
				for (int c = 0; c < dimensionsCond; c++) {
					condObservations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
			}
		}
	}

	/**
	 * Utility function to compute the norms between each pair of points in each marginal time series
	 *
	 */
	protected void computeNorms() {
		int N = var1Observations.length; // number of observations
		
		xNorms = new double[N][N];
		yNorms = new double[N][N];
		zNorms = new double[N][N];
		for (int t = 0; t < N; t++) {
			// Compute the norms from t to all other time points
			double[][] xyzNormsForT = normCalculator.computeNorms(var1Observations,
					var2Observations, condObservations, t);
			for (int t2 = 0; t2 < N; t2++) {
				xNorms[t][t2] = xyzNormsForT[t2][0];
				yNorms[t][t2] = xyzNormsForT[t2][1];
				zNorms[t][t2] = xyzNormsForT[t2][2];
			}
		}
	}
	
	/**
	 * @returns the series of local conditional MI values in nats, not bits.
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] states1,
			double[][] states2, double[][] condStates) throws Exception {
		// If implemented, will need to incorporate any normalisation here
		//  (normalising the incoming data the same way the previously
		//   supplied observations were normalised).
		throw new Exception("Local method not implemented yet");
	}
	
	public abstract String printConstants(int N) throws Exception;

	// Note: no extra implementation of clone provided; we're simply
	//  allowing clone() to produce a shallow copy, which is fine
	//  for the statistical significance calculation (none of the array
	//  data will be changed there.
	//
	// public ConditionalMutualInfoCalculatorMultiVariateKraskov clone() {
	//	return this;
	// }
}
