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

import infodynamics.measures.continuous.ConditionalTransferEntropyCalculatorViaCondMutualInfo;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;

/**
 * 
 * <p>
 * Implements a conditional transfer entropy calculator using model of 
 * Gaussian variables with linear interactions.
 * This is equivalent (up to a multiplicative constant) to 
 * (a conditional) Granger causality (see Barnett et al., below).
 * This is achieved by plugging in {@link ConditionalMutualInfoCalculatorMultiVariateGaussian}
 * as the calculator into {@link ConditionalTransferEntropyCalculatorViaCondMutualInfo}.
 * </p>
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct: {@link #ConditionalTransferEntropyCalculatorGaussian()}</li>
 * 		<li>Set properties: {@link #setProperty(String, String)} for each relevant property, including those
 * 			of either {@link ConditionalTransferEntropyCalculatorViaCondMutualInfo#setProperty(String, String)}
 * 			or {@link ConditionalMutualInfoCalculatorMultiVariateGaussian#setProperty(String, String)}.</li>
 *		<li>Initialise: by calling one of {@link #initialise()} etc.</li>
 * 		<li>Add observations to construct the PDFs: {@link #setObservations(double[], double[], double[][])},
 * 			or [{@link #startAddObservations()},
 * 			{@link #addObservations(double[], double[], double[][])}*, {@link #finaliseAddObservations()}]
 *   		Note: If not using setObservations(), the results from computeLocal
 *   		will be concatenated directly, and getSignificance will mix up observations 
 *          from separate trials (added in separate {@link #addObservations(double[])} calls.</li> 
 * 		<li>Compute measures: e.g. {@link #computeAverageLocalOfObservations()} or
 * 			{@link #computeLocalOfPreviousObservations()} etc </li>
 * 	</ol>
 * </p>
 * 
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 * 
 * @see "Schreiber, Physical Review Letters 85 (2) pp.461-464 (2000);
 *  <a href='http://dx.doi.org/10.1103/PhysRevLett.85.461'>download</a>
 *  (for definition of transfer entropy)"
 * @see "Lizier, Prokopenko and Zomaya, Physical Review E 77, 026110 (2008);
 * <a href='http://dx.doi.org/10.1103/PhysRevE.77.026110'>download</a>
 *  (for the extension to <i>conditional</i> transfer entropy 
 *  or <i>complete</i> where all other causal sources are conditioned on,
 *  and <i>local</i> transfer entropy)"
 * @see "Lizier, Prokopenko and Zomaya, Chaos 20, 3, 037109 (2010);
 * <a href='http://dx.doi.org/10.1063/1.3486801'>download</a>
 *  (for further clarification on <i>conditional</i> transfer entropy 
 *  or <i>complete</i> where all other causal sources are conditioned on)"
 * @see "Lionel Barnett, Adam B. Barrett, Anil K. Seth, Physical Review Letters 103 (23) 238701, 2009;
 *  <a href='http://dx.doi.org/10.1103/physrevlett.103.238701'>download</a>
 *  (for direct relation between transfer entropy and Granger causality)"
 *  
 * @see ConditionalTransferEntropyCalculator
 *
 */
public class ConditionalTransferEntropyCalculatorGaussian
	extends ConditionalTransferEntropyCalculatorViaCondMutualInfo
	implements AnalyticNullDistributionComputer {
	
	public static final String COND_MI_CALCULATOR_GAUSSIAN = ConditionalMutualInfoCalculatorMultiVariateGaussian.class.getName();
	
	/**
	 * Creates a new instance of the Gaussian-estimate style conditional transfer entropy calculator
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ConditionalTransferEntropyCalculatorGaussian() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(COND_MI_CALCULATOR_GAUSSIAN);
	}

	/**
	 * <p>Set the joint covariance of the distribution for which we will compute the
	 *  transfer entropy.</p>
	 *  
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][])}.</p>
	 * 
	 * @param covariance joint covariance matrix of source, dest, dest history
	 *  and conditional variables, considered together.
	 * @param numObservations the number of observations that the covariance
	 *  was determined from. This is used for later significance calculations
	 * @throws Exception for covariance matrix not matching the expected dimensions,
	 *  being non-square, asymmetric or non-positive definite
	 */
	public void setCovariance(double[][] covariance, int numObservations) throws Exception {
		((ConditionalMutualInfoCalculatorMultiVariateGaussian) condMiCalc).
				setCovariance(covariance, numObservations);
	}

	/**
	 * <p>Compute the statistical significance of the TE 
	 *  result analytically, without creating a distribution
	 *  under the null hypothesis by bootstrapping.
	 *  Computed using the corresponding method of the
	 *  underlying 
	 *  {@link ConditionalMutualInfoCalculatorMultiVariateGaussian}</p>
	 *  
	 * @see {@link ConditionalMutualInfoCalculatorMultiVariateGaussian#computeSignificance()}
	 * @return ChiSquareMeasurementDistribution object 
	 *  This object contains the proportion of TE scores from the distribution
	 *  which have higher or equal TEs to ours.
	 */
	public ChiSquareMeasurementDistribution computeSignificance()
			throws Exception {
		return ((ConditionalMutualInfoCalculatorMultiVariateGaussian) condMiCalc).computeSignificance();
	}
}
