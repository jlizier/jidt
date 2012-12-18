package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;

public class MutualInfoMultiVariateTester
	extends infodynamics.measures.continuous.MutualInfoMultiVariateAbstractTester {

	/**
	 * Utility function to create a calculator for the given algorithm number
	 * 
	 * @param algNumber
	 * @return
	 */
	public MutualInfoCalculatorMultiVariateKraskov getNewCalc(int algNumber) {
		MutualInfoCalculatorMultiVariateKraskov miCalc = null;;
		if (algNumber == 1) {
			miCalc = new MutualInfoCalculatorMultiVariateKraskov1();
		} else if (algNumber == 2) {
			miCalc = new MutualInfoCalculatorMultiVariateKraskov2();
		}
		return miCalc;
	}
	
	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void checkLocalsAverageCorrectly(int algNumber) throws Exception {
		
		MutualInfoCalculatorMultiVariateKraskov miCalc = getNewCalc(algNumber);
		
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
		
		MutualInfoCalculatorMultiVariateKraskov miCalc = getNewCalc(algNumber);
		
		String kraskov_K = "4";
		
		miCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testComputeSignificanceDoesntAlterAverage(miCalc, 2, 100);
	}
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		checkComputeSignificanceDoesntAlterAverage(1);
		checkComputeSignificanceDoesntAlterAverage(2);
	}
	
	/**
	 * Test the computed MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIhigherdim <dataFile> 2 1 1 3000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMIforRandomVariablesFromFile() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/2randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {1, 2, 3, 4, 5, 6, 10, 15};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA = {-0.05294175, -0.03944338, -0.02190217,
				0.00120807, -0.00924771, -0.00316402, -0.00778205, -0.00565778};
		
		// The Kraskov MILCA toolkit MIhigherdim executable 
		//  uses algorithm 2 by default (this is what it means by rectangular):
		MutualInfoCalculatorMultiVariateKraskov miCalc = getNewCalc(2);
		
		System.out.println("Kraskov comparison 1 - random data");
		for (int kIndex = 0; kIndex < kNNs.length; kIndex++) {
			int k = kNNs[kIndex];
			miCalc.setProperty(
					MutualInfoCalculatorMultiVariateKraskov.PROP_K,
					Integer.toString(k));
			miCalc.setObservations(MatrixUtils.selectColumns(data, new int[] {0}),
					MatrixUtils.selectColumns(data, new int[] {1}));
			double mi = miCalc.computeAverageLocalOfObservations();
			//miCalc.setDebug(false);
			
			System.out.printf("k=%d: Average MI for 2 column %d row Kraskov data was %.8f (expected %.8f)\n",
					k, miCalc.getNumObservations(), mi, expectedFromMILCA[kIndex]);
			assertEquals(expectedFromMILCA[kIndex], mi, 0.00000001);			
		}
	}
}
