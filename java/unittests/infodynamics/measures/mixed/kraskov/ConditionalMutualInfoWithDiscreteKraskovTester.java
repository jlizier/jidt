package infodynamics.measures.mixed.kraskov;

import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class ConditionalMutualInfoWithDiscreteKraskovTester extends TestCase {

	public void testUnivariateArrayMethods() {
		ConditionalMutualInfoCalculatorMultiVariateWithDiscreteKraskov cmiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateWithDiscreteKraskov();
		boolean exceptionCaught = false;
		int dataLength = 1000;
		RandomGenerator rg = new RandomGenerator();
		double[] univariateContData1 = rg.generateNormalData(dataLength, 0, 1);
		double[][] multivariateActualUniContData1 = rg.generateNormalData(dataLength, 1, 0, 1);
		double[][] multivariateContData1 = rg.generateNormalData(dataLength, 2, 0, 1);
	    int[] discData = rg.generateRandomInts(dataLength, 2);
		double[] univariateCondData = rg.generateNormalData(dataLength, 0, 1);
		double[][] multivariateActualUniCondData = rg.generateNormalData(dataLength, 1, 0, 1);
		double[][] multivariateCondData = rg.generateNormalData(dataLength, 2, 0, 1);
		
		// Testing when univariate methods should work
		cmiCalc.initialise(1, 2, 1);
		try {
			cmiCalc.setObservations(univariateContData1, discData, univariateCondData);
			cmiCalc.initialise(1, 2, 1);
			cmiCalc.setObservations(multivariateActualUniContData1, discData, multivariateActualUniCondData);
		} catch (Exception e) {
			exceptionCaught = true;
		}
		assertFalse(exceptionCaught);
		// and now supply actual multivariates here which is wrong:
		cmiCalc.initialise(1, 2, 1);
		try{
			cmiCalc.setObservations(multivariateContData1, discData, multivariateCondData);
		} catch (Exception e) {
			exceptionCaught = true;
		}
		assertTrue(exceptionCaught);
		exceptionCaught = false;

		// Testing when univariate methods should not work
		cmiCalc.initialise(2, 2, 1);
		try {
			cmiCalc.setObservations(univariateContData1, discData, univariateCondData);
		} catch (Exception e) {
			exceptionCaught = true;
		}
		assertTrue(exceptionCaught);
		exceptionCaught = false;
		cmiCalc.initialise(1, 2, 2);
		try {
			cmiCalc.setObservations(univariateContData1, discData, univariateCondData);
		} catch (Exception e) {
			exceptionCaught = true;
		}
		assertTrue(exceptionCaught);
		exceptionCaught = false;
		cmiCalc.initialise(2, 2, 2);
		try {
			cmiCalc.setObservations(univariateContData1, discData, univariateCondData);
		} catch (Exception e) {
			exceptionCaught = true;
		}
		assertTrue(exceptionCaught);
		exceptionCaught = false;
		// and now supply the proper multivariate arrays:
		cmiCalc.initialise(2, 2, 2);
		try {
			cmiCalc.setObservations(multivariateContData1, discData, multivariateCondData);
		} catch (Exception e) {
			exceptionCaught = true;
		}
		assertFalse(exceptionCaught);		
	}

}
