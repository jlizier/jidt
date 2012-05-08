package infodynamics.measures.continuous.kraskov;

/**
 * <p>Compute the Transfer Entropy for multi-variates using the Kraskov estimation method.<p>
 * <p>This calculator extends the {@link TransferEntropyCalculatorKraskovByMulti TransferEntropyCalculatorKraskov}
 *  by allowing multivariate sources and destinations. See 
 *  {@link TransferEntropyCalculatorKraskovByMulti TransferEntropyCalculatorKraskov}
 *  for further comments on the implementation of Transfer entropy via the Kraskov 
 *  Mutual Information estimation.</p>
 * </p>
 * 
 * 
 * @author Joseph Lizier; joseph.lizier at gmail.com
 *
 */
public class TransferEntropyCalculatorMultiVariateKraskov
	extends TransferEntropyCalculatorMultiVariateKraskovByMulti {

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
