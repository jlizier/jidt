package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateAbstractTester;

public class ConditionalMutualInfoMultiVariateTester extends
		ConditionalMutualInfoMultiVariateAbstractTester {

	public void testLocalsAverageCorrectly() throws Exception {
		ConditionalMutualInfoCalculatorMultiVariateGaussian condMiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateGaussian();
		super.testLocalsAverageCorrectly(condMiCalc, 2, 100);
	}

	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		ConditionalMutualInfoCalculatorMultiVariateGaussian condMiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateGaussian();
		super.testComputeSignificanceDoesntAlterAverage(condMiCalc, 2, 100);
	}
	
}
