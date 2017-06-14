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
	 * Construct the distribution.
	 * Note: the Chi squared distribution is technically
	 * of 2*numObservations*(the info theoretic estimate), not
	 * of the info theoretic measurement itself.
	 * 
	 * @param actualValue actual observed information-theoretic value
	 * @param numObservations the number of observations that the information theoretic estimate
	 * 	is computed from
	 * @param degreesOfFreedom degrees of freedom for the distribution
	 */
	public ChiSquareMeasurementDistribution(double actualValue, 
			int numObservations, int degreesOfFreedom) {
		super(actualValue, 1 - MathsUtils.chiSquareCdf(2.0*((double)numObservations)*actualValue, degreesOfFreedom));
		this.numObservations = numObservations;
		this.degreesOfFreedom = degreesOfFreedom;
		chi2dist = new ChiSquaredDistribution(degreesOfFreedom);
	}
	
	public double computePValueForGivenEstimate(double estimate) {
		return 1 - MathsUtils.chiSquareCdf(2.0*((double)numObservations)*estimate, degreesOfFreedom);
	}
	
	public double computeEstimateForGivenPValue(double pValue) {
		return chi2dist.inverseCumulativeProbability(1 - pValue) / (2.0*((double)numObservations));
		// Could also call the following, but this doesn't re-use our objects:
		// return MathsUtils.chiSquareInv(1 - pValue, degreesOfFreedom);
	}
}
