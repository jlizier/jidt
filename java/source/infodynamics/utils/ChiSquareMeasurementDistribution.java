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

package infodynamics.utils;

import infodynamics.utils.commonsmath3.distribution.ChiSquaredDistribution;

/**
 * Class to represent analytic distributions of info theoretic measurements under
 * some null hypothesis of a relationship between the variables, where the
 * distribution of <b>a function of</b> those information-theoretic measurements
 * is a Chi Square distribution.
 * Can represent an analytic distribution of raw estimates or bias corrected
 * estimates.
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class ChiSquareMeasurementDistribution extends
		AnalyticMeasurementDistribution {

	/**
	 * Number of degrees of freedom for the distribution
	 */
	protected int degreesOfFreedom;
	
	/**
	 * The number of observations that the information theoretic estimate
	 * is computed from
	 */
	protected int numObservations;
	
	/**
	 * An object of the chi2dist class from commons.math3,
	 *  which we'll use to access the distribution values
	 */
	protected ChiSquaredDistribution chi2dist;
	
	/**
	 * Tracks whether we will return a bias-corrected distribution or not
	 */
	protected boolean isBiasCorrected;
	
	/**
	 * Stores the mean of the uncorrected distribution, useful for bias correction externally
	 */
	protected double meanOfUncorrectedDistribution;
	
	/**
	 * Construct the distribution.
	 * Note: the Chi squared distribution is technically
	 * of 2*numObservations*(the info theoretic estimate), not
	 * of the info theoretic measurement itself.
	 * This constructor assumes that we are not seeking a bias-corrected distribution.
	 * 
	 * @param actualValue actual observed information-theoretic value
	 * @param numObservations the number of observations that the information theoretic estimate
	 * 	is computed from
	 * @param degreesOfFreedom degrees of freedom for the distribution
	 */
	public ChiSquareMeasurementDistribution(double actualValue, 
			int numObservations, int degreesOfFreedom) {
		this(actualValue, numObservations, degreesOfFreedom, false);
	}
	
	/**
	 * Construct the distribution.
	 * Note: the Chi squared distribution is technically
	 * of 2*numObservations*(the info theoretic estimate), not
	 * of the info theoretic measurement itself.
	 * 
	 * @param actualValue actual observed information-theoretic value
	 * @param numObservations the number of observations that the information theoretic estimate
	 * 	is computed from
	 * @param degreesOfFreedom degrees of freedom for the distribution
	 * @param isBiasCorrected whether to bias correct the distribution or not (and indeed
	 *   whether the actualValue is bias corrected)
	 */
	public ChiSquareMeasurementDistribution(double actualValue, 
			int numObservations, int degreesOfFreedom, boolean isBiasCorrected) {
		// Make a dummy initialisation until we can properly get the uncorrected distribution ready:
		super(0, 0);
		this.numObservations = numObservations;
		this.degreesOfFreedom = degreesOfFreedom;
		chi2dist = new ChiSquaredDistribution(degreesOfFreedom); // Uncorrected distribution
		meanOfUncorrectedDistribution = chi2dist.getNumericalMean() / (2.0*((double)numObservations));
		this.isBiasCorrected = isBiasCorrected;
		this.actualValue = actualValue; // will be bias corrected, if we are bias correcting
		// Now we can properly compute the p-value of this potentially bias corrected actual value
		this.pValue = computePValueForGivenEstimate(actualValue);
	}

	public double computePValueForGivenEstimate(double estimate) {
		if (isBiasCorrected) {
			// estimate is biasCorrected, so we need to add back into it the meanOfUncorrectedDistribution
			return 1 - MathsUtils.chiSquareCdf(2.0*((double)numObservations)*(estimate + meanOfUncorrectedDistribution), degreesOfFreedom);
		} else {
			return 1 - MathsUtils.chiSquareCdf(2.0*((double)numObservations)*estimate, degreesOfFreedom);
		}
	}

	public double computeEstimateForGivenPValue(double pValue) {
		double uncorrectedEstimate = chi2dist.inverseCumulativeProbability(1 - pValue) / (2.0*((double)numObservations));
		if (isBiasCorrected) {
			return uncorrectedEstimate - meanOfUncorrectedDistribution;
		} else {
			return uncorrectedEstimate;
		}
		// Could also call the following, but this doesn't re-use our objects:
		// return MathsUtils.chiSquareInv(1 - pValue, degreesOfFreedom);
	}
	
	/**
	 * Compute the mean of the measurement distribution
	 * 
	 * @return the mean
	 */
	public double getMeanOfDistribution() {
		if (isBiasCorrected) {
			return 0;
		} else {
			return meanOfUncorrectedDistribution;
		}
	}
	
	/**
	 * 
	 * @return mean of uncorrected distribution
	 */
	public double getMeanOfUncorrectedDistribution() {
		return meanOfUncorrectedDistribution;
	}

	/**
	 * Compute the standard deviation of the measurement distribution
	 * 
	 * @return the standard deviation
	 */
	public double getStdOfDistribution() {
		return Math.sqrt(chi2dist.getNumericalVariance()) / (2.0*((double)numObservations));
	}
}
