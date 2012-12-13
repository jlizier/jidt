package infodynamics.measures.continuous.kernel;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class TransferEntropyMultiVariateTester extends TestCase {

	
	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testLocalsAverageCorrectly() throws Exception {
		
		TransferEntropyCalculatorMultiVariateKernel teCalc =
				new TransferEntropyCalculatorMultiVariateKernel();
		
		int dimensions = 2;
		int timeSteps = 100;
		int k = 1;
		String kernelWidth = "1";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.NORMALISE_PROP_NAME,
				"true");
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.EPSILON_PROP_NAME,
				kernelWidth);
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

		assertEquals(te, MatrixUtils.mean(teLocal, k, timeSteps-k), 0.0001);
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
		
		int dimensions = 2;
		int timeSteps = 100;
		int k = 1;
		String kernelWidth = "1";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.NORMALISE_PROP_NAME,
				"true");
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.EPSILON_PROP_NAME,
				kernelWidth);
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

}
