package infodynamics.measures.continuous;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public abstract class MutualInfoMultiVariateAbstractTester extends TestCase {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * @param miCalc a pre-constructed MutualInfoCalculatorMultiVariate object
	 * @param dimensions number of dimensions for the source and dest data to use
	 * @param timeSteps number of time steps for the random data
	 */
	public void testLocalsAverageCorrectly(MutualInfoCalculatorMultiVariate miCalc,
			int dimensions, int timeSteps)
			throws Exception {
		
		miCalc.initialise(dimensions, dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		miCalc.setObservations(sourceData, destData);
		
		//teCalc.setDebug(true);
		double mi = miCalc.computeAverageLocalOfObservations();
		//miCalc.setDebug(false);
		double[] miLocal = miCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", mi);

		assertEquals(mi, MatrixUtils.mean(miLocal), 0.00001);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @param miCalc a pre-constructed MutualInfoCalculatorMultiVariate object
	 * @param dimensions number of dimensions for the source and dest data to use
	 * @param timeSteps number of time steps for the random data
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage(MutualInfoCalculatorMultiVariate miCalc,
			int dimensions, int timeSteps) throws Exception {
		
		miCalc.initialise(dimensions, dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		miCalc.setObservations(sourceData, destData);
		
		//miCalc.setDebug(true);
		double mi = miCalc.computeAverageLocalOfObservations();
		//miCalc.setDebug(false);
		//double[] miLocal = miCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", mi);
		
		// Now look at statistical significance tests
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(
				timeSteps, 100);

		miCalc.computeSignificance(newOrderings);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double averageCheck1 = miCalc.computeAverageLocalOfObservations();
			assertEquals(mi, averageCheck1);
		}
	}

}
