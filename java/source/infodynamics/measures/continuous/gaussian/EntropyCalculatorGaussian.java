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

package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.EntropyCalculator;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential entropy of a given set of observations
 *  (implementing {@link EntropyCalculator}, assuming that
 *  the probability distribution function for these observations is Gaussian.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link EntropyCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #EntropyCalculatorGaussian()}.</li>
 * 	<li>The user can call {@link #setVariance(double)}
 *     instead of supplying observations via {@link #setObservations(double[])}.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991).</li>
    <li>Differential entropy for Gaussian random variables defined at 
 *      <a href="http://mathworld.wolfram.com/DifferentialEntropy.html">MathWorld</a></li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EntropyCalculatorGaussian implements EntropyCalculator {

	/**
	 * Variance of the most recently supplied observations, or set directly
	 */
	protected double variance;
	
	/**
	 * Whether we are in debug mode
	 */
	protected boolean debug;
	
	/**
	 * The set of observations, retained in case the user wants to retrieve the local
	 *  entropy values of these
	 */
	protected double[] observations;
	
	/**
	 * Store the last computed average Entropy
	 */
	protected double lastAverage;
	
	/**
	 * Construct an instance
	 */
	public EntropyCalculatorGaussian() {
		// Nothing to do
	}
	
	@Override
	public void initialise() {
		observations = null;
		variance = 0;
	}

	public void setObservations(double[] observations) {
		variance = MatrixUtils.stdDev(observations);
		variance *= variance;
		this.observations = observations;
	}

	/**
	 * An alternative to {@link #setObservations(double[])}, allowing user to
	 * set the variance of the distribution for which we will compute the
	 *  entropy.
	 * 
	 * @param variance the variance of the univariate distribution.
	 */
	public void setVariance(double variance) throws Exception {
		this.variance = variance;
		if (variance < 0) {
			throw new Exception("Cannot have negative variance");
		}
		observations = null;
	}
	
	/**
	 * Compute the entropy from the previously supplied observations, or
	 * based on the supplied variance.
	 * 
	 * <p>The entropy for a Gaussian-distribution random variable with
	 *  variance \sigma is 0.5*\log_e{2*pi*e*\sigma}.</p>
	 * 
	 * <p>Here we compute the entropy assuming that the recorded estimation of the
	 *  variance is correct (i.e. we will not make a bias correction for limited
	 *  observations here).</p>
	 * 
	 * @return the entropy of the previously provided observations or from the supplied
	 *   covariance matrix. Entropy returned in <b>nats</b>, not bits!
	 */
	@Override
	public double computeAverageLocalOfObservations() {
		lastAverage = 0.5 * Math.log(2.0*Math.PI*Math.E*variance);
		return lastAverage;
	}

	/**
	 * @throws Exception if {@link #setVariance(double)} was used previously instead
	 * of {@link #setObservations(double[][])}
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		if (observations == null) {
			throw new Exception("Cannot compute local values since no observations were supplied");
		}
		
		// Check that the variance was non-zero:
		if (variance == 0) {
			throw new Exception("variance is not positive - cannot compute local entropies");
		}
		// Now we are clear to take the variance inverse
		double invVariance = 1.0 / variance;
		double mean = MatrixUtils.mean(observations);
		
		double[] localValues = new double[observations.length];
		for (int t = 0; t < observations.length; t++) {
			double deviationFromMean = observations[t] - mean;
			// Computing PDF
			// (see the PDF defined at the wikipedia page referenced in the method header)
			double jointExpArg = deviationFromMean * deviationFromMean * invVariance;
			double pJoint = Math.exp(-0.5 * jointExpArg) /
					Math.sqrt(2.0 * Math.PI * variance);
			localValues[t] = - Math.log(pJoint);
		}
		
		// Don't set average if this was the previously supplied observations,
		//  since it won't be the same as what would have been computed 
		//  analytically.
		
		return localValues;
	}
	
	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	/**
	 * No properties are defined here, so this method will have no effect.
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		// No properties to set here
	}

	/**
	 * No properties are defined here, so this method will always return null.
	 */
	@Override
	public String getProperty(String propertyName)
			throws Exception {
		// No properties to return here
		return null;
	}

	@Override
	public int getNumObservations() throws Exception {
		if (observations == null) {
			throw new Exception("Cannot return number of observations because either " +
					"this calculator has not had observations supplied or " +
					"the user supplied the variance instead of observations");
		}
		return observations.length;
	}

	@Override
	public double getLastAverage() {
		return lastAverage;
	}
}
