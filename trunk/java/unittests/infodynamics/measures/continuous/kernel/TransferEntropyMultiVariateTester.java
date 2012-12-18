package infodynamics.measures.continuous.kernel;

public class TransferEntropyMultiVariateTester
	extends infodynamics.measures.continuous.TransferEntropyMultiVariateAbstractTester {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testLocalsAverageCorrectly() throws Exception {
		
		TransferEntropyCalculatorMultiVariateKernel teCalc =
				new TransferEntropyCalculatorMultiVariateKernel();
		
		String kernelWidth = "1";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.NORMALISE_PROP_NAME,
				"true");
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.EPSILON_PROP_NAME,
				kernelWidth);

		super.testLocalsAverageCorrectly(teCalc, 2, 100, 1);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		TransferEntropyCalculatorMultiVariateKernel teCalc =
				new TransferEntropyCalculatorMultiVariateKernel();
		
		String kernelWidth = "1";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.NORMALISE_PROP_NAME,
				"true");
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.EPSILON_PROP_NAME,
				kernelWidth);

		super.testComputeSignificanceDoesntAlterAverage(teCalc, 2, 100, 1);
	}

}
