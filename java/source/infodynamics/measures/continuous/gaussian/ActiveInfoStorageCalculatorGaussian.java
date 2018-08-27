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

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.measures.continuous.ActiveInfoStorageCalculatorViaMutualInfo;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * An Active Information Storage (AIS) calculator (implementing {@link ActiveInfoStorageCalculator})
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
 * <p>Usage is as per the paradigm outlined for {@link ActiveInfoStorageCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #ActiveInfoStorageCalculatorGaussian()}.</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties for
 *      {@link MutualInfoCalculatorMultiVariateGaussian#setProperty(String, String)}
 *      (except {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}
 *      and those of the super-class as outlined
 *      in {@link ActiveInfoStorageCalculatorViaMutualInfo#setProperty(String, String)}),
 *      including auto-embedding properties.</li>
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
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see ActiveInfoStorageCalculator
 * @see ActiveInfoStorageCalculatorViaMutualInfo
 * @see MutualInfoCalculatorMultiVariateGaussian
 */
public class ActiveInfoStorageCalculatorGaussian
	extends ActiveInfoStorageCalculatorViaMutualInfo
	implements AnalyticNullDistributionComputer {
	
	public static final String MI_CALCULATOR_GAUSSIAN = MutualInfoCalculatorMultiVariateGaussian.class.getName();
	
	/**
	 * Property name for the number of surrogates to use in computing the bias correction
	 *  if required for the {@link ActiveInfoStorageCalculatorViaMutualInfo#AUTO_EMBED_METHOD_MAX_CORR_AIS}
	 *  auto embedding method. Defaults to 0 meaning that we use analytic bias correction rather than empirical
	 *  surrogates. Note: This is not used for bias correction of the raw values, only for auto-embedding
	 */
	public static final String PROP_MAX_CORR_AIS_NUM_SURROGATES = "AUTO_EMBED_MAX_CORR_AIS_SURROGATES";
	/**
	 * Internal variable for storing the number of surrogates to use for the
	 * auto-embedding {@link ActiveInfoStorageCalculatorViaMutualInfo#AUTO_EMBED_METHOD_MAX_CORR_AIS}
	 * method. 0 mean we use analytic approaches rather than surrogates.
	 */
	protected int auto_embed_num_surrogates = 0;

	/**
	 * Creates a new instance of the Gaussian-estimate style active info storage calculator
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorGaussian() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(MI_CALCULATOR_GAUSSIAN);
	}
	
	/**
	 * Sets properties for the AIS Gaussian calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method.
	 *
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #PROP_MAX_CORR_AIS_NUM_SURROGATES} -- number of surrogates to use
	 * 		to compute the bias correction
	 * 		in the auto-embedding if the property {@link #PROP_AUTO_EMBED_METHOD}
	 * 		has been set to {@link #AUTO_EMBED_METHOD_MAX_CORR_AIS}. Defaults to 0
	 * 		meaning that we use analytic bias correction.
	 * 		Note: this is not used for other bias-correction, only inside auto-embedding</li>
	 * 		<li>Any properties accepted by {@link super#setProperty(String, String)}</li>
	 * 		<li>Or properties accepted by the underlying
	 * 		{@link MutualInfoCalculatorMultiVariateGaussian#setProperty(String, String)} implementation.</li>
	 * </ul>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value).
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_MAX_CORR_AIS_NUM_SURROGATES)) {
			auto_embed_num_surrogates = Integer.parseInt(propertyValue);
		} else {
			propertySet = false;
			// Assume it was a property for the parent class or underlying MI calculator
			super.setProperty(propertyName, propertyValue);
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public String getProperty(String propertyName)
			throws Exception {
		
		if (propertyName.equalsIgnoreCase(PROP_MAX_CORR_AIS_NUM_SURROGATES)) {
			return Integer.toString(auto_embed_num_surrogates);
		} else {
			// Assume it was a property for the parent class or underlying MI calculator
			return super.getProperty(propertyName);
		}
	}

	@Override
	protected double computeAdditionalBiasToRemove() throws Exception {
		boolean biasCorrected = Boolean.getBoolean(getProperty(MutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION));
		if (auto_embed_num_surrogates == 0) {
			// Analytic bias correction:
			if (!biasCorrected) {
				ChiSquareMeasurementDistribution analyticMeasDist =
						((MutualInfoCalculatorMultiVariateGaussian)miCalc).computeSignificance();
				return analyticMeasDist.getMeanOfDistribution();
			} else {
				return 0;
			}
		} else {
			// Empirical bias correction with auto_embed_num_surrogates surrogates:
			EmpiricalMeasurementDistribution measDist =
					miCalc.computeSignificance(auto_embed_num_surrogates);
			return measDist.getMeanOfDistribution();
		}
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
