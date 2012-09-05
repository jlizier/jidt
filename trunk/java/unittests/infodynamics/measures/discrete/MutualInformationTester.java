package infodynamics.measures.discrete;

import junit.framework.TestCase;

public class MutualInformationTester extends TestCase {

	public void testFullyDependent() {
		MutualInformationCalculator miCalc = new MutualInformationCalculator(2, 0);
		
		// X2 is a copy of X1 - MI should be 1 bit
		miCalc.initialise();
		miCalc.addObservations(new int[] {0, 0, 1, 1}, new int[] {0, 0, 1, 1});
		double miCopy = miCalc.computeAverageLocalOfObservations();
		assertEquals(1.0, miCopy, 0.000001);
	}
	
	public void testIndependent() {
		MutualInformationCalculator miCalc = new MutualInformationCalculator(2, 0);
		
		// X2 is unrelated to X1 - MI should be 0 bits
		miCalc.initialise();
		miCalc.addObservations(new int[] {0, 0, 1, 1}, new int[] {0, 1, 0, 1});
		double miRand = miCalc.computeAverageLocalOfObservations();
		assertEquals(0.0, miRand, 0.000001);
	}

	public void testXor() {
		MutualInformationCalculator miCalc = new MutualInformationCalculator(2, 0);
		
		int[] X1 = new int[] {0, 0, 1, 1};
		int[] X2 = new int[] {0, 1, 0, 1};
		int[] Y  = new int[] {0, 1, 1, 0};
		
		// Y is independent of X1 - MI should be 0 bits
		miCalc.initialise();
		miCalc.addObservations(X1, Y);
		double miX1Y = miCalc.computeAverageLocalOfObservations();
		assertEquals(0.0, miX1Y, 0.000001);

		// Y is independent of X2 - MI should be 0 bits
		miCalc.initialise();
		miCalc.addObservations(X2, Y);
		double miX2Y = miCalc.computeAverageLocalOfObservations();
		assertEquals(0.0, miX2Y, 0.000001);
		
		// Y is fully determined from X1, X2 - MI should be 0 bits
		MutualInformationCalculator miCalcBase4 = new MutualInformationCalculator(4, 0);
		int[] X12 = new int[] {0, 1, 2, 3};
		miCalcBase4.initialise();
		miCalcBase4.addObservations(X12, Y);
		miCalcBase4.setDebug(true);
		double miX12Y = miCalcBase4.computeAverageLocalOfObservations();
		assertEquals(1.0, miX12Y, 0.000001);
	}
}
