package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariateViaCondMutualInfo;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;

/**
 * 
 * <p>
 * Implements a multivariate transfer entropy calculator using model of 
 * Gaussian variables with linear interactions.
 * This is equivalent (up to a multiplicative constant) to 
 * Granger causality (see Barnett et al., below).
 * This is achieved by plugging in {@link ConditionalMutualInfoCalculatorMultiVariateGaussian}
 * as the calculator into {@link TransferEntropyCalculatorMultiVariateViaCondMutualInfo}.
 * </p>
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct: {@link #TransferEntropyCalculatorMultiVariateGaussian()}</li>
 * 		<li>Set properties: {@link #setProperty(String, String)} for each relevant property, including those
 * 			of either {@link TransferEntropyCalculatorMultiVariateViaCondMutualInfo#setProperty(String, String)}
 * 			or {@link ConditionalMutualInfoCalculatorMultiVariateGaussian#setProperty(String, String)}.</li>
 *		<li>Initialise: by calling one of {@link #initialise()} etc.</li>
 * 		<li>Add observations to construct the PDFs: {@link #setObservations(double[])}, or [{@link #startAddObservations()},
 * 			{@link #addObservations(double[])}*, {@link #finaliseAddObservations()}]
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
 * @see "Lionel Barnett, Adam B. Barrett, Anil K. Seth, Physical Review Letters 103 (23) 238701, 2009;
 *  <a href='http://dx.doi.org/10.1103/physrevlett.103.238701'>download</a>
 *  (for direct relation between transfer entropy and Granger causality)"
 * @see "J.T. Lizier, J. Heinzle, A. Horstmann, J.-D. Haynes, M. Prokopenko,
 * Journal of Computational Neuroscience, vol. 30, pp. 85-107, 2011
 * <a href='http://dx.doi.org/10.1007/s10827-010-0271-2'>download</a>
 * (for definition of <i>multivariate</i> transfer entropy"
 *  
 * @see TransferEntropyCalculatorMultiVariate
 *
 */
public class TransferEntropyCalculatorMultiVariateGaussian
	extends TransferEntropyCalculatorMultiVariateViaCondMutualInfo
	implements AnalyticNullDistributionComputer {
	
	public static final String COND_MI_CALCULATOR_GAUSSIAN = ConditionalMutualInfoCalculatorMultiVariateGaussian.class.getName();
	
	/**
	 * Creates a new instance of the Gaussian-estimate style transfer entropy calculator
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public TransferEntropyCalculatorMultiVariateGaussian() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
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
	 * @param covariance joint covariance matrix of the multivariate
	 *  source, dest, dest history
	 *  variables, considered together.
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
