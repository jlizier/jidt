package infodynamics.measures.continuous;

import junit.framework.TestCase;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

public abstract class TransferEntropyMultiVariateAbstractTester extends TestCase {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * @param teCalc a pre-constructed TransferEntropyCalculatorMultiVariate object
	 * @param dimensions number of dimensions for the source and dest data to use
	 * @param timeSteps number of time steps for the random data
	 * @param k history length for the TE calculator to use
	 */
	public void testLocalsAverageCorrectly(TransferEntropyCalculatorMultiVariate teCalc,
			int dimensions, int timeSteps, int k)
			throws Exception {
		
		teCalc.initialise(k, dimensions, dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		teCalc.setObservations(sourceData, destData);
		
		//teCalc.setDebug(true);
		double te = teCalc.computeAverageLocalOfObservations();
		//teCalc.setDebug(false);
		double[] teLocal = teCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", te);

		assertEquals(te, MatrixUtils.mean(teLocal, k, timeSteps-k), 0.00001);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @param teCalc a pre-constructed TransferEntropyCalculatorMultiVariate object
	 * @param dimensions number of dimensions for the source and dest data to use
	 * @param timeSteps number of time steps for the random data
	 * @param k history length for the TE calculator to use
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage(TransferEntropyCalculatorMultiVariate teCalc,
			int dimensions, int timeSteps, int k) throws Exception {
		
		teCalc.initialise(k, dimensions, dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		teCalc.setObservations(sourceData, destData);
		
		//teCalc.setDebug(true);
		double te = teCalc.computeAverageLocalOfObservations();
		//teCalc.setDebug(false);
		//double[] teLocal = teCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", te);
		
		// Now look at statistical significance tests
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(
				timeSteps - k, 100);

		teCalc.computeSignificance(newOrderings);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double averageCheck1 = teCalc.computeAverageLocalOfObservations();
			assertEquals(te, averageCheck1);
		}
	}

	/**
	 * Confirm that a calculation for univariate data using univariate method signatures
	 *  matches that with multivariate signatures.
	 * 
	 * @param teCalc a pre-constructed TransferEntropyCalculatorMultiVariate object
	 * @param timeSteps number of time steps for the random data
	 * @param k history length for the TE calculator to use
	 */
	public void testUnivariateMatchesMultivariateRoute(TransferEntropyCalculatorMultiVariate teCalc,
			int timeSteps, int k)
			throws Exception {
				
		if (!(teCalc instanceof TransferEntropyCalculator)) {
			throw new Exception("The given calculator does not implement univariate TE");
		}
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, 1,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, 1,
				0, 1);
		
		// Compute via univariate signatures:
		TransferEntropyCalculator teCalcUni = (TransferEntropyCalculator) teCalc;
		teCalc.initialise(k, 1, 1);
		teCalcUni.setObservations(MatrixUtils.selectColumn(sourceData, 0),
				MatrixUtils.selectColumn(destData, 0));
		double teUnivariate = teCalc.computeAverageLocalOfObservations();

		// compute via multivariate signatures:
		teCalc.initialise(k, 1, 1);
		teCalc.setObservations(sourceData, destData);
		//teCalc.setDebug(true);
		double teMultivariate = teCalc.computeAverageLocalOfObservations();
		//teCalc.setDebug(false);
		
		assertEquals(teUnivariate, teMultivariate, 0.00001);
	}
}
