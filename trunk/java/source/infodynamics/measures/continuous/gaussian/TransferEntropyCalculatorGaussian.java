package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.TransferEntropyCalculatorViaCondMutualInfo;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;

/**
 * 
 * <p>
 * Implements a transfer entropy calculator using model of 
 * Gaussian variables with linear interactions.
 * This is equivalent (up to a multiplicative constant) to 
 * Granger causality (see Barnett et al., below).
 * This is achieved by plugging in {@link ConditionalMutualInfoCalculatorMultiVariateGaussian}
 * as the calculator into {@link TransferEntropyCalculatorViaCondMutualInfo}.
 * </p>
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct: {@link #TransferEntropyCalculatorGaussian()}</li>
 * 		<li>Set properties: {@link #setProperty(String, String)} for each relevant property, including those
 * 			of either {@link TransferEntropyCalculatorViaCondMutualInfo#setProperty(String, String)}
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
 * @author Joseph Lizier
 * @see "Lionel Barnett, Adam B. Barrett, Anil K. Seth, Physical Review Letters 103 (23) 238701, 2009;
 *  <a href='http://dx.doi.org/10.1103/physrevlett.103.238701'>download</a>
 *  (for direct relation between transfer entropy and Granger causality)"
 *  
 * @see TransferEntropyCalculator
 *
 */
public class TransferEntropyCalculatorGaussian
	extends TransferEntropyCalculatorViaCondMutualInfo 
	implements AnalyticNullDistributionComputer {
	
	public static final String COND_MI_CALCULATOR_GAUSSIAN = ConditionalMutualInfoCalculatorMultiVariateGaussian.class.getName();
	
	/**
	 * Creates a new instance of the Gaussian-estimate style transfer entropy calculator
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public TransferEntropyCalculatorGaussian() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
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
	
	/**
	 * Debug method to check the computed determinants
	 * 
	 * @return an array of the four relevant determinants
	 *  computed by the underlying {@link ConditionalMutualInfoCalculatorMultiVariateGaussian}:
	 *  of the whole covariance matrix, of the source and conditionals,
	 *  of the destination and conditionals, and of the conditionals themselves.
	 */
	public double[] getDeterminants() {
		ConditionalMutualInfoCalculatorMultiVariateGaussian condMiGaussian = 
				(ConditionalMutualInfoCalculatorMultiVariateGaussian) condMiCalc;
		double[] determinants = new double[4];
		determinants[0] = condMiGaussian.detCovariance;
		determinants[1] = condMiGaussian.det1cCovariance;
		determinants[2] = condMiGaussian.det2cCovariance;
		determinants[3] = condMiGaussian.detccCovariance;
		return determinants;
	}
}
