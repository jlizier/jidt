package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

import junit.framework.TestCase;

public class ConditionalMutualInformationTester extends TestCase {

	protected RandomGenerator rand = new RandomGenerator();
	
	protected int numObservations = 1000;
	
	public void testComputeSignificanceInt() {
		// Just making sure that no exception is thrown
		ConditionalMutualInformationCalculator condMiCalc = new ConditionalMutualInformationCalculator(2, 2, 2);
		int[] x1 = rand.generateRandomInts(numObservations, 2);
		int[] x2 = rand.generateRandomInts(numObservations, 2);
		int[] cond = rand.generateRandomInts(numObservations, 2);
		condMiCalc.initialise();
		condMiCalc.addObservations(x1, x2, cond);
		condMiCalc.computeAverageLocalOfObservations();
		condMiCalc.computeSignificance(1000);
	}

	public void testComputeLocalsGivesCorrectAverage() {
		ConditionalMutualInformationCalculator condMiCalc = new ConditionalMutualInformationCalculator(2, 2, 2);
		int[] x1 = rand.generateRandomInts(numObservations, 2);
		int[] x2 = rand.generateRandomInts(numObservations, 2);
		int[] cond = rand.generateRandomInts(numObservations, 2);
		condMiCalc.initialise();
		condMiCalc.addObservations(x1, x2, cond);
		double average = condMiCalc.computeAverageLocalOfObservations();
		double avLastAverage = condMiCalc.getLastAverage();
		double[] locals = condMiCalc.computeLocal(x1, x2, cond);
		double localsLastAverage = condMiCalc.getLastAverage();
		assertEquals(average, MatrixUtils.mean(locals), 0.000001);
		assertEquals(avLastAverage, localsLastAverage, 0.000001);
		assertEquals(average, avLastAverage, 0.000001);
	}

	public void testSetDebug() {
		ConditionalMutualInformationCalculator condMiCalc = new ConditionalMutualInformationCalculator(2, 2, 2);
		assertFalse(condMiCalc.debug);
		condMiCalc.setDebug(true);
		assertTrue(condMiCalc.debug);
		condMiCalc.setDebug(false);
		assertFalse(condMiCalc.debug);
	}

}
