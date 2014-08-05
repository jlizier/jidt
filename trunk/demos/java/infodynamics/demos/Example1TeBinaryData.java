package infodynamics.demos;

import infodynamics.utils.RandomGenerator;
import infodynamics.measures.discrete.TransferEntropyCalculator;

/**
 * 
 * = Example 1 - Transfer entropy on binary data =
 * 
 * Simple transfer entropy (TE) calculation on binary data using the discrete TE calculator.
 * 
 * @author Joseph Lizier
 *
 */
public class Example1TeBinaryData {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		int arrayLengths = 100;
		RandomGenerator rg = new RandomGenerator();
		
		// Generate some random binary data:
		int[] sourceArray = rg.generateRandomInts(arrayLengths, 2);
		int[] destArray = new int[arrayLengths];
		destArray[0] = 0;
		System.arraycopy(sourceArray, 0, destArray, 1, arrayLengths - 1);
		int[] sourceArray2 = rg.generateRandomInts(arrayLengths, 2);
		
		// Create a TE calculator and run it:
		TransferEntropyCalculator teCalc=
				new TransferEntropyCalculator(2, 1);
		teCalc.initialise();
		teCalc.addObservations(destArray, sourceArray);
		double result = teCalc.computeAverageLocalOfObservations();
		System.out.printf("For copied source, result should be close to 1 bit : %.3f bits\n", result);
		teCalc.initialise();
		teCalc.addObservations(destArray, sourceArray2);
		double result2 = teCalc.computeAverageLocalOfObservations();
		System.out.printf("For random source, result should be close to 0 bits: %.3f bits\n", result2);
	}

}
