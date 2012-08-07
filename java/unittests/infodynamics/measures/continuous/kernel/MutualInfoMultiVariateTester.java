package infodynamics.measures.continuous.kernel;

import junit.framework.TestCase;
import infodynamics.utils.RandomGenerator;

public class MutualInfoMultiVariateTester extends TestCase {

	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		MutualInfoCalculatorMultiVariateKernel miCalc =
				new MutualInfoCalculatorMultiVariateKernel();
		
		int dimensions = 2;
		int timeSteps = 100;
		double kernelWidth = 1;
		
		miCalc.setProperty(
				MutualInfoCalculatorMultiVariateKernel.NORMALISE_PROP_NAME,
				"true");
		miCalc.initialise(dimensions, dimensions, kernelWidth);
		
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
