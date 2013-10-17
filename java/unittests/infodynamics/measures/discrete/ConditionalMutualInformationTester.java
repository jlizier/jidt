package infodynamics.measures.discrete;

import infodynamics.utils.RandomGenerator;

import junit.framework.TestCase;

public class ConditionalMutualInformationTester extends TestCase {

	protected RandomGenerator rand = new RandomGenerator();
	
	protected int numObservations = 100;
	
	public void testComputeSignificanceInt() {
		ConditionalMutualInformationCalculator condMiCalc = new ConditionalMutualInformationCalculator(2, 2, 2);
		int[] x1 = rand.generateRandomInts(numObservations, 2);
		int[] x2 = rand.generateRandomInts(numObservations, 2);
		int[] cond = rand.generateRandomInts(numObservations, 2);
		condMiCalc.initialise();
		condMiCalc.addObservations(x1, x2, cond);
		condMiCalc.computeAverageLocalOfObservations();
		condMiCalc.computeSignificance(1000);
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
