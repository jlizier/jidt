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

import infodynamics.measures.continuous.ActiveInfoStorageCalculatorMultiVariate;
import infodynamics.measures.continuous.ActiveInfoStorageCalculatorMultiVariateViaMutualInfo;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;

/**
 * An Active Information Storage (AIS) calculator for multivariate time-series
 * data (implementing {@link ActiveInfoStorageCalculatorMultiVariate})
 * which is affected using a 
 * Gaussian Mutual Information (MI) calculator
 * ({@link MutualInfoCalculatorMultiVariateGaussian}) to make the calculations.
 * 
 * <p>
 * That is, this class implements an AIS calculator using model of 
 * Gaussian variables with linear interactions.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateGaussian}
 * as the calculator into the parent class {@link ActiveInfoStorageCalculatorViaMutualInfo}.
 * </p> 
 * 
 * <p>Usage is as per the paradigm outlined for {@link ActiveInfoStorageCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #ActiveInfoStorageCalculatorMultiVariateGaussian()}.</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties for
 *      {@link MutualInfoCalculatorMultiVariateGaussian#setProperty(String, String)}
 *      (except {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} as outlined
 *      in {@link ActiveInfoStorageCalculatorMultiVariateViaMutualInfo#setProperty(String, String)})</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>J.T. Lizier, M. Prokopenko and A.Y. Zomaya,
 * 		<a href="http://dx.doi.org/10.1016/j.ins.2012.04.016">
 * 		"Local measures of information storage in complex distributed computation"</a>,
 * 		Information Sciences, vol. 208, pp. 39-54, 2012.</li>
 * </ul>
 * 
 * @author Pedro AM Mediano (<a href="pmediano at imperial.ac.uk">email</a>,
 * <a href="https://www.doc.ic.ac.uk/~pam213/">www</a>)
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 *
 * @see ActiveInfoStorageCalculatorMultiVariate
 * @see ActiveInfoStorageCalculatorMultiVariateViaMutualInfo
 * @see MutualInfoCalculatorMultiVariateGaussian
 */
public class ActiveInfoStorageCalculatorMultiVariateGaussian
	extends ActiveInfoStorageCalculatorMultiVariateViaMutualInfo {
	
	public static final String MI_CALCULATOR_GAUSSIAN = MutualInfoCalculatorMultiVariateGaussian.class.getName();
	
	/**
	 * Creates a new instance of the Gaussian-estimate style active info storage calculator
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorMultiVariateGaussian() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(MI_CALCULATOR_GAUSSIAN);
	}

	@Override
	protected double computeAdditionalBiasToRemove() throws Exception {
		boolean biasCorrected = Boolean.getBoolean(getProperty(MutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION));
		if (!biasCorrected) {
			ChiSquareMeasurementDistribution analyticMeasDist =
					((MutualInfoCalculatorMultiVariateGaussian)miCalc).computeSignificance();
			return analyticMeasDist.getMeanOfDistribution();
		}
		// else it was already bias corrected
		return 0;
	}

	/**
	 * Generate an <b>analytic</b> distribution of what the AIS would look like,
	 * under a null hypothesis that our variables had no relation.
	 * This is performed without bootstrapping (which is done in
	 * {@link #computeSignificance(int, int)} and {@link #computeSignificance(int, int[][])}).
	 * The method is implemented using the corresponding method of the
	 *  underlying {@link MutualInfoCalculatorMultiVariateGaussian}
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below, and the other papers referenced in
	 * {@link AnalyticNullDistributionComputer#computeSignificance()}
	 * (in particular Brillinger and Geweke),
	 * for a description of how this is done for MI.
	 * Basically, the null distribution is a chi-square distribution.
	 * </p>
	 * 
	 * @return ChiSquareMeasurementDistribution object which describes
	 * the proportion of AIS scores from the null distribution
	 *  which have higher or equal AISs to our actual value.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception
	 */
	public ChiSquareMeasurementDistribution computeSignificance() throws Exception {
		return ((MutualInfoCalculatorMultiVariateGaussian) miCalc).computeSignificance();
	}
}
