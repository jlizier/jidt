package infodynamics.measures.continuous.gaussian;

import junit.framework.TestCase;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

public class MutualInfoMultiVariateTester extends TestCase {

	public void testCovarianceDoesntMatchDimensions(){
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		miCalc.initialise(1, 1);
		boolean caughtException = false;
		// Check that we catch a covariance matrix which doesn't have
		//  the right number of rows
		try {
			miCalc.setCovariance(new double[][] {{2,1,0.5}, {1,2,0.5}, {0.5,0.5,2}});
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		// Check that we catch a covariance matrix which isn't square
		caughtException = false;
		try {
			miCalc.setCovariance(new double[][] {{2,1}, {1,2,0.5}});
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		// Check that we catch a covariance matrix which isn't symmetric
		caughtException = false;
		try {
			miCalc.setCovariance(new double[][] {{2,1}, {1.000001,2}});
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		// Check that no exception is thrown for an ok covariance matrix
		caughtException = false;
		double[][] goodCovariance = new double[][] {{2,1}, {1,2}};
		try {
			miCalc.setCovariance(goodCovariance);
		} catch (Exception e) {
			caughtException = true;
		}
		assertFalse(caughtException);
		// and that this covariance has been set:
		assertEquals(goodCovariance, miCalc.covariance);
	}
	
	public void testAnalyticMatchesCouplingValue() throws Exception {
		double[][] covarianceMatrix = {{1.0, 0}, {0, 1.0}};
		
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		
		for (double covar = 0.0; covar < 1.0; covar += 0.01) {
			// Insert the new covariance into the matrix:
			covarianceMatrix[0][1] = covar;
			covarianceMatrix[1][0] = covar;

			miCalc.initialise(1, 1);
			miCalc.setCovariance(covarianceMatrix);
			assertEquals(-0.5 * Math.log(1.0 - covar*covar),
					miCalc.computeAverageLocalOfObservations(), 0.00000000001);
		}
	}
	
	
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
