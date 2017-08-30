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

/**
 * Generic class to represent analytic distributions of info theoretic measurements under
 * some null hypothesis of a relationship between the variables.
 * This class is not designed to be creatable - it's children (which represent
 * various concrete distributions) are those which should be created.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class AnalyticMeasurementDistribution extends MeasurementDistribution {

	/**
	 * Construct the instance
	 * 
	 * @param actualValue observed value of the information-theoretic measure
	 * @param pValue p-value that surrogate measurements are greater than
	 *  the observed value. (1 - CDF of null at the observed value)
	 */
	protected AnalyticMeasurementDistribution(double actualValue, double pValue) {
		super(actualValue, pValue);
	}
	
	/**
	 * Compute the <b>analytic</b> p-value 
	 * for the given estimate (i.e. the input argument, not the 
	 * estimate produced by any particular observations supplied)
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination value (possibly
	 * in the context of a conditional value).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * conditional MI and TE as per the references below.
	 * </p>
	 * 
	 * <p><b>References:</b><br/>
	 *  <ul>
	 *   <li>J.T. Lizier, "JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems", 2014.</li>
	 * 	 <li>Brillinger, <a href="http://www.stat.berkeley.edu/~brill/Papers/MIBJPS.pdf">
	 * 		"Some data analyses using mutual information"</a>,
	 * 		Brazilian Journal of Probability and Statistics, <b>18</b>, p. 163, (2004)</li>
	 * 	 <li>Cheng et al., <a href="http://www.jds-online.com/file_download/112/JDS-369.pdf">
	 * 		"Data Information in Contingency Tables: A Fallacy of Hierarchical Loglinear Models"</a>,
	 * 		Journal of Data Science, <b>4</b>, p. 387 (2006).</li>
	 *   <li>Geweke, <a href="http://dx.doi.org/10.1080/01621459.1982.10477803">
	 *   	"Measurement of Linear Dependence and Feedback between Multiple Time Series"</a>,
	 *   	Journal of the American Statistical Association, <b>77</b>, p. 304-313 (1982).</li>
	 * 	 <li>Barnett and Bossomaier, <a href="http://arxiv.org/abs/1205.6339">
	 * 		"Transfer Entropy as a Log-likelihood Ratio"</a>,
	 * 		Physical Review Letters, <b>109</b>, p. 138105+ (2012).</li>
     * </ul>
	 * 
	 * @return p-value for the given channel measure score under this null hypothesis.
	 * @throws Exception
	 */
	public abstract double computePValueForGivenEstimate(double estimate);

	/**
	 * Computes p-values corresponding to a set of estimates, each done via
	 * {@link #computePValueForGivenEstimate(double)}
	 * 
	 * @param estimates array of estimates to return corresponding p-values for.
	 * @return
	 */
	public double[] computePValuesForGivenEstimates(double[] estimates) {
		double[] pValues = new double[estimates.length];
		for (int i = 0; i < estimates.length; i++) {
			pValues[i] = computePValueForGivenEstimate(estimates[i]);
		}
		return estimates;
	}

	/**
	 * Compute the estimated observed measured value corresponding to
	 * a given p-value
	 * (derived from how the given estimates are <b>analytically</b> distributed  
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination value, possibly
	 * in the context of a conditional value).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * conditional MI and TE as per the references below.
	 * </p>
	 * 
	 * <p><b>References:</b><br/>
	 *  <ul>
	 *   <li>J.T. Lizier, "JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems", 2014.</li>
	 * 	 <li>Brillinger, <a href="http://www.stat.berkeley.edu/~brill/Papers/MIBJPS.pdf">
	 * 		"Some data analyses using mutual information"</a>,
	 * 		Brazilian Journal of Probability and Statistics, <b>18</b>, p. 163, (2004)</li>
	 * 	 <li>Cheng et al., <a href="http://www.jds-online.com/file_download/112/JDS-369.pdf">
	 * 		"Data Information in Contingency Tables: A Fallacy of Hierarchical Loglinear Models"</a>,
	 * 		Journal of Data Science, <b>4</b>, p. 387 (2006).</li>
	 *   <li>Geweke, <a href="http://dx.doi.org/10.1080/01621459.1982.10477803">
	 *   	"Measurement of Linear Dependence and Feedback between Multiple Time Series"</a>,
	 *   	Journal of the American Statistical Association, <b>77</b>, p. 304-313 (1982).</li>
	 * 	 <li>Barnett and Bossomaier, <a href="http://arxiv.org/abs/1205.6339">
	 * 		"Transfer Entropy as a Log-likelihood Ratio"</a>,
	 * 		Physical Review Letters, <b>109</b>, p. 138105+ (2012).</li>
     * </ul>
	 * 
	 * @param pValue the sample p-value for the given channel measure score under this null hypothesis.
	 * @return the estimate of the channel measure score corresponding to
	 * the given p-value under this null hypothesis.
	 * @throws Exception
	 */
	public abstract double computeEstimateForGivenPValue(double pValue);
	
	/**
	 * Computes estimates corresponding to a set of p-values, each done via
	 * {@link #computeEstimateForGivenPValue(double)}
	 * 
	 * @param pValues array of p-values to return corresponding estimates for.
	 * @return
	 */
	public double[] computeEstimatesForGivenPValues(double[] pValues) {
		double[] estimates = new double[pValues.length];
		for (int i = 0; i < pValues.length; i++) {
			estimates[i] = computeEstimateForGivenPValue(pValues[i]);
		}
		return estimates;
	}
}
