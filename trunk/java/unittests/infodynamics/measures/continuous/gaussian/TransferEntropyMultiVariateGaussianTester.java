package infodynamics.measures.continuous.gaussian;

public class TransferEntropyMultiVariateGaussianTester
	extends infodynamics.measures.continuous.TransferEntropyMultiVariateAbstractTester {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testLocalsAverageCorrectly() throws Exception {
		
		TransferEntropyCalculatorMultiVariateGaussian teCalc =
				new TransferEntropyCalculatorMultiVariateGaussian();
		
		super.testLocalsAverageCorrectly(teCalc, 2, 100, 1);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		TransferEntropyCalculatorMultiVariateGaussian teCalc =
				new TransferEntropyCalculatorMultiVariateGaussian();
		
		super.testComputeSignificanceDoesntAlterAverage(teCalc, 2, 100, 1);
	}

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testUnivariateSignatureMatchesMultivariate() throws Exception {
		
		TransferEntropyCalculatorMultiVariateGaussian teCalc =
				new TransferEntropyCalculatorMultiVariateGaussian();
		
		super.testUnivariateMatchesMultivariateRoute(teCalc, 100, 1);
	}
}
