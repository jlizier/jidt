package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;

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
	
	/**
	 * Utility function to run Kraskov conditional MI algorithm 1
	 *  as transfer entropy for data with known results
	 *  from TRENTOOL.
	 * 
	 * @param var1 source multivariate data set
	 * @param var2 dest multivariate data set
	 * @param kNNs array of Kraskov k nearest neighbours parameter to check
	 * @param expectedResults array of expected results for each k
	 */
	protected void checkTEForGivenData(double[][] var1, double[][] var2,
			int[] kNNs, double[] expectedResults) throws Exception {
				
		ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc = getNewCalc(1);
		
		// Normalise the data ourselves rather than letting the calculator do it -
		//  this ensures the extra values in the time series (e.g. last value in source)
		//  are taken into account, in line with TRENTOOL
		var1 = MatrixUtils.normaliseIntoNewArray(var1);
		var2 = MatrixUtils.normaliseIntoNewArray(var2);
		
		for (int kIndex = 0; kIndex < kNNs.length; kIndex++) {
			int k = kNNs[kIndex];
			condMiCalc.setProperty(
					ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
					Integer.toString(k));
			// We already normalised above, and this will do a different
			//  normalisation without taking the extra values in to account if we did it
			condMiCalc.setProperty(
					ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE,
					Boolean.toString(false));
			// No longer need to set this property as it's set by default:
			//condMiCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
			//		EuclideanUtils.NORM_MAX_NORM_STRING);
			condMiCalc.setObservations(MatrixUtils.selectRows(var1, 0, var1.length - 1),
					MatrixUtils.selectRows(var2, 1, var2.length - 1),
					MatrixUtils.selectRows(var2, 0, var2.length - 1));
			double condMi = condMiCalc.computeAverageLocalOfObservations();
			//miCalc.setDebug(false);
			
			System.out.printf("k=%d: Average MI %.8f (expected %.8f)\n",
					k, condMi, expectedResults[kIndex]);
			// 6 decimal places is Matlab accuracy
			assertEquals(expectedResults[kIndex], condMi, 0.000001);			
		}
	}

	/**
	 * Test the computed univariate TE as a conditional MI
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 1, 1, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=1,
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testUnivariateTEforCoupledVariablesFromFile() throws Exception {
		
		// Test set 1:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {0.3058006};
		
		System.out.println("Kraskov Cond MI as TE comparison 1 - univariate coupled data 1");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromTRENTOOL);
		
		// And now in the reverse direction:
		expectedFromTRENTOOL = new double[] {-0.0029744};
		
		System.out.println("  reverse direction:");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {0}),
				kNNs, expectedFromTRENTOOL);

	}

	/**
	 * Test the computed univariate TE as a conditional MI
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 1, 1, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=1,
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testUnivariateTEforCoupledLogisticMapFromFile() throws Exception {
		
		// Test set 1:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/coupledLogisticMapXY.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {0.508417};
		
		System.out.println("Kraskov Cond MI as TE comparison 1 - univariate coupled logistic map data 1");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromTRENTOOL);
		
		// And now in the reverse direction:
		expectedFromTRENTOOL = new double[] {0.016257};
		
		System.out.println("  reverse direction:");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {0}),
				kNNs, expectedFromTRENTOOL);

	}

	/**
	 * Test the computed univariate TE as a conditional MI
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 1, 1, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=1,
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testUnivariateTEforRandomDataFromFile() throws Exception {
		
		// Test set 1:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {-0.0096556};
		
		System.out.println("Kraskov Cond MI as TE comparison 1 - univariate random data 1");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromTRENTOOL);
		
		// And now for other columns
		expectedFromTRENTOOL = new double[] {0.0175389};
		
		System.out.println("  reverse direction:");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {2}),
				kNNs, expectedFromTRENTOOL);

	}
}
