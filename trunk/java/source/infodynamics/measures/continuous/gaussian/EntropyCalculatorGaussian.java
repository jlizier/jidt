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
 * <p>Computes the differential entropy of a given set of observations, assuming that
 *  the probability distribution function for these observations is Gaussian.</p>
 *  
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct</li>
 *		<li>{@link #initialise()}</li>
 * 		<li>Either set the observations: {@link #setObservations()}, or 
 * 			directly set the variance {@link #setVariance()}.</li> 
 * 		<li>{@link #computeAverageLocalOfObservations()} to return the average differential
 *          entropy based on either the set variance or the variance of
 *          the supplied observations.</li>
 * 	</ol>
 * </p>
 * 
 * @see Differential entropy for Gaussian random variables defined at 
 *      {@link http://mathworld.wolfram.com/DifferentialEntropy.html}
 * @author Joseph Lizier joseph.lizier_at_gmail.com
 *
 */
public class EntropyCalculatorGaussian implements EntropyCalculator {

	/**
	 * Variance of the most recently supplied observations
	 */
	protected double variance;
	
	protected boolean debug;
	
	/**
	 * Constructor
	 */
	public EntropyCalculatorGaussian() {
		// Nothing to do
	}
	
	/**
	 * Initialise the calculator ready for reuse
	 */
	public void initialise() {
		// Nothing to do
	}

	/**
	 * Provide the observations from which to compute the entropy
	 * 
	 * @param observations the observations to compute the entropy from
	 */
	public void setObservations(double[] observations) {
		variance = MatrixUtils.stdDev(observations);
		variance *= variance;
	}

	/**
	 * Set the variance of the distribution for which we will compute the
	 *  entropy.
	 * 
	 * @param variance
	 */
	public void setVariance(double variance) {
		this.variance = variance;
	}
	
	/**
	 * <p>The entropy for a Gaussian-distribution random variable with
	 *  variance \sigma is 0.5*\log_e{2*pi*e*\sigma}.</p>
	 * 
	 * <p>Here we compute the entropy assuming that the recorded estimation of the
	 *  variance is correct (i.e. we will not make a bias correction for limited
	 *  observations here).</p>
	 * 
	 * @return the entropy of the previously provided observations or from the supplied
	 *   covariance matrix. Entropy returned in nats, not bits!
	 */
	public double computeAverageLocalOfObservations() {
		return 0.5 * Math.log(2.0*Math.PI*Math.E*variance);
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		// No properties to set here
	}

}
