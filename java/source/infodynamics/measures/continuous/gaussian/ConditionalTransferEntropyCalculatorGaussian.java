package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.ConditionalTransferEntropyCalculatorViaCondMutualInfo;

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
	extends ConditionalTransferEntropyCalculatorViaCondMutualInfo {
	
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
}
