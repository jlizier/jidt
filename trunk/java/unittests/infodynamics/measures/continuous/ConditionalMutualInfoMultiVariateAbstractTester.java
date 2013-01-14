package infodynamics.measures.continuous;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public abstract class ConditionalMutualInfoMultiVariateAbstractTester
	extends TestCase {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * @param condMiCalc a pre-constructed ConditionalMutualInfoCalculatorMultiVariate object
	 * @param dimensions number of dimensions for the source and dest data to use
	 * @param timeSteps number of time steps for the random data
	 */
	public void testLocalsAverageCorrectly(ConditionalMutualInfoCalculatorMultiVariate condMiCalc,
			int dimensions, int timeSteps)
			throws Exception {
		
		condMiCalc.initialise(dimensions, dimensions, dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] condData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		condMiCalc.setObservations(sourceData, destData, condData);
		
		//teCalc.setDebug(true);
		double condmi = condMiCalc.computeAverageLocalOfObservations();
		//miCalc.setDebug(false);
		double[] condMiLocal = condMiCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", condmi);

		assertEquals(condmi, MatrixUtils.mean(condMiLocal), 0.00001);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @param condMiCalc a pre-constructed ConditionalMutualInfoCalculatorMultiVariate object
	 * @param dimensions number of dimensions for the source and dest data to use
	 * @param timeSteps number of time steps for the random data
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage(ConditionalMutualInfoCalculatorMultiVariate condMiCalc,
			int dimensions, int timeSteps) throws Exception {
		
		condMiCalc.initialise(dimensions, dimensions, dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] condData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		condMiCalc.setObservations(sourceData, destData, condData);
		
		//condMiCalc.setDebug(true);
		double condMi = condMiCalc.computeAverageLocalOfObservations();
		//condMiCalc.setDebug(false);
		//double[] condMiLocal = miCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", condMi);
		
		// Now look at statistical significance tests
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(
				timeSteps, 100);
		
		// Compute significance for permuting first variable 
		condMiCalc.computeSignificance(1, newOrderings);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double averageCheck1 = condMiCalc.computeAverageLocalOfObservations();
			assertEquals(condMi, averageCheck1);
		}

		// Compute significance for permuting second variable 
		condMiCalc.computeSignificance(2, newOrderings);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double averageCheck1 = condMiCalc.computeAverageLocalOfObservations();
			assertEquals(condMi, averageCheck1);
		}
	}

}
