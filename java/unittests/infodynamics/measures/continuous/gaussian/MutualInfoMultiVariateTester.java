package infodynamics.measures.continuous.gaussian;

import junit.framework.TestCase;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

public class MutualInfoMultiVariateTester extends TestCase {

	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		
		int dimensions = 2;
		int timeSteps = 100;
		
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

		EmpiricalMeasurementDistribution measDist =
				miCalc.computeSignificance(newOrderings);
		
		System.out.printf("pValue of sig test was %.3f\n", measDist.pValue);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double averageCheck1 = miCalc.computeAverageLocalOfObservations();
			assertEquals(mi, averageCheck1);
		}
	}
}
