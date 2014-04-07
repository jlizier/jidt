package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.TransferEntropyCalculatorViaCondMutualInfo;

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
	extends TransferEntropyCalculatorViaCondMutualInfo {
	
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
}
