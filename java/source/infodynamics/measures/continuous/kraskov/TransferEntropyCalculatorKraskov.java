package infodynamics.measures.continuous.kraskov;

/**
 * <p>Compute the Transfer Entropy using the Kraskov estimation method for underlying
 * mutual information calculation, using direct calculation of the mutual information
 * rather than breaking it down into component multi-informations.<p>
 * 
 * @see TransferEntropyCalculatorKraskovByMulti
 * @see MutualInfoCalculatorMultiVariateDirect
 * @see MutualInfoCalculatorMultiVariateDirect1
 * @see MutualInfoCalculatorMultiVariateDirect2
 * @author Joseph Lizier joseph.lizier at gmail.com
 *
 */
public class TransferEntropyCalculatorKraskov
	extends TransferEntropyCalculatorKraskovByMulti {

	/**
	 * Create the underlying Kraskov MI calculators using direct multi-variate MI calculators
	 *
	 */
	protected void createKraskovMiCalculators() {
		if (kraskovAlgorithmNumber == 1) {
			mickPastToSource = new MutualInfoCalculatorMultiVariateKraskov1();
			mickNextPastToSource = new MutualInfoCalculatorMultiVariateKraskov1();
		} else {
			// Algorithm 2
			mickPastToSource = new MutualInfoCalculatorMultiVariateKraskov2();
			mickNextPastToSource = new MutualInfoCalculatorMultiVariateKraskov2();
		}
	}
	
	protected void shareDataBetweenUnderlyingCalculators() {
		// Nothing to do here
	}
}
