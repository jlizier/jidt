package infodynamics.measures.continuous.kraskov;

public class ConditionalMutualInfoMultiVariateTester
	extends infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateAbstractTester {

	/**
	 * Utility function to create a calculator for the given algorithm number
	 * 
	 * @param algNumber
	 * @return
	 */
	public ConditionalMutualInfoCalculatorMultiVariateKraskov getNewCalc(int algNumber) {
		ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc = null;
		if (algNumber == 1) {
			condMiCalc = new ConditionalMutualInfoCalculatorMultiVariateKraskov1();
		} else if (algNumber == 2) {
			condMiCalc = new ConditionalMutualInfoCalculatorMultiVariateKraskov2();
		}
		return condMiCalc;
	}
	
	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void checkLocalsAverageCorrectly(int algNumber) throws Exception {
		
		ConditionalMutualInfoCalculatorMultiVariateKraskov miCalc = getNewCalc(algNumber);
		
		String kraskov_K = "4";
		
		miCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testLocalsAverageCorrectly(miCalc, 2, 100);
	}
	public void testLocalsAverageCorrectly() throws Exception {
		checkLocalsAverageCorrectly(1);
		checkLocalsAverageCorrectly(2);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @throws Exception
	 */
	public void checkComputeSignificanceDoesntAlterAverage(int algNumber) throws Exception {
		
		ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc = getNewCalc(algNumber);
		
		String kraskov_K = "4";
		
		condMiCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testComputeSignificanceDoesntAlterAverage(condMiCalc, 2, 100);
	}
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		checkComputeSignificanceDoesntAlterAverage(1);
		checkComputeSignificanceDoesntAlterAverage(2);
	}
	
}
