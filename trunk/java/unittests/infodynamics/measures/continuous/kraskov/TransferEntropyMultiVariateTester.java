package infodynamics.measures.continuous.kraskov;

public class TransferEntropyMultiVariateTester
	extends infodynamics.measures.continuous.TransferEntropyMultiVariateAbstractTester {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testLocalsAverageCorrectly() throws Exception {
		
		TransferEntropyCalculatorMultiVariateKraskov teCalc =
				new TransferEntropyCalculatorMultiVariateKraskov();
		
		String kraskov_K = "4";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKraskov.PROP_KRASKOV_ALG_NUM,
				"2");
		teCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testLocalsAverageCorrectly(teCalc, 2, 100, 1);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		TransferEntropyCalculatorMultiVariateKraskov teCalc =
				new TransferEntropyCalculatorMultiVariateKraskov();
		
		String kraskov_K = "4";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKraskov.PROP_KRASKOV_ALG_NUM,
				"2");
		teCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testComputeSignificanceDoesntAlterAverage(teCalc, 2, 100, 1);
	}

}
