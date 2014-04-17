package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.ConditionalTransferEntropyAbstractTester;

public class ConditionalTransferEntropyKraskovTester extends
	ConditionalTransferEntropyAbstractTester {

	public void testUnivariateMethodSignatureFails() throws Exception {
		ConditionalTransferEntropyCalculatorKraskov teCalc = new ConditionalTransferEntropyCalculatorKraskov();
		String kraskov_K = "4";
		teCalc.setProperty(
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);
		super.testUnivariateCallFailsIfWrongInitialisation(teCalc);
	}
	
	public void testLocalsAverageCorrectly() throws Exception {
		ConditionalTransferEntropyCalculatorKraskov teCalc = new ConditionalTransferEntropyCalculatorKraskov();
		String kraskov_K = "4";
		teCalc.setProperty(
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);
		super.testLocalsAverageCorrectly(teCalc, 100, 2);
	}

	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		ConditionalTransferEntropyCalculatorKraskov teCalc = new ConditionalTransferEntropyCalculatorKraskov();
		String kraskov_K = "4";
		teCalc.setProperty(
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);
		super.testComputeSignificanceDoesntAlterAverage(teCalc, 100, 2);
	}

	public void testAgainstApparentTE() throws Exception {
		TransferEntropyCalculatorKraskov teCalc = new TransferEntropyCalculatorKraskov();
		String kraskov_K = "4";
		teCalc.setProperty(
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);
		ConditionalTransferEntropyCalculatorKraskov condTeCalc = new ConditionalTransferEntropyCalculatorKraskov();
		condTeCalc.setProperty(
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);
		testConditionalAgainstOrdinaryTE(teCalc, condTeCalc, 100, 2);
		testConditionalAgainstOrdinaryTE(teCalc, condTeCalc, 100, 3);
		testConditionalAgainstOrdinaryTE(teCalc, condTeCalc, 100, 4);
	}
}
