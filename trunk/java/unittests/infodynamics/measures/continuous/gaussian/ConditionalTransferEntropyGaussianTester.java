package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.ConditionalTransferEntropyAbstractTester;

public class ConditionalTransferEntropyGaussianTester extends
	ConditionalTransferEntropyAbstractTester {

	public void testUnivariateMethodSignatureFails() throws Exception {
		ConditionalTransferEntropyCalculatorGaussian teCalc = new ConditionalTransferEntropyCalculatorGaussian();
		super.testUnivariateCallFailsIfWrongInitialisation(teCalc);
	}
	
	public void testLocalsAverageCorrectly() throws Exception {
		ConditionalTransferEntropyCalculatorGaussian teCalc = new ConditionalTransferEntropyCalculatorGaussian();
		super.testLocalsAverageCorrectly(teCalc, 100, 2);
	}

	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		ConditionalTransferEntropyCalculatorGaussian teCalc = new ConditionalTransferEntropyCalculatorGaussian();
		super.testComputeSignificanceDoesntAlterAverage(teCalc, 100, 2);
	}

	public void testAgainstApparentTE() throws Exception {
		TransferEntropyCalculatorGaussian teCalc = new TransferEntropyCalculatorGaussian();
		ConditionalTransferEntropyCalculatorGaussian condTeCalc = new ConditionalTransferEntropyCalculatorGaussian();
		testConditionalAgainstOrdinaryTE(teCalc, condTeCalc, 100, 2);
		testConditionalAgainstOrdinaryTE(teCalc, condTeCalc, 100, 3);
		testConditionalAgainstOrdinaryTE(teCalc, condTeCalc, 100, 4);
	}
}
