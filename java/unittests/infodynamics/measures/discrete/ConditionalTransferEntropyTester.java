
package infodynamics.measures.discrete;

import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class ConditionalTransferEntropyTester extends TestCase {

	protected RandomGenerator rand = new RandomGenerator();
  protected double tol = 1e-5;

  /**
   * Source and other are random binary ints, and dest is equal to other.
   * Since other doesn't explain anything apart from dest's past, the result
   * should be the same as ignoring other and using TransferEntropyCalculator.
   */
	public void testNoInfoForPredictableDest() {
    int nb_samples = 100;
		int[] source = rand.generateRandomInts(nb_samples, 2);
		int[] other = rand.generateRandomInts(nb_samples, 2);
		int[] dest = other;
    boolean exceptionCaught = false;
    double cteResult = Double.NaN, teResult = Double.NaN;
		
    try {
      ConditionalTransferEntropyCalculatorDiscrete cteCalc =
          new ConditionalTransferEntropyCalculatorDiscrete(2, 1, 1);
      cteCalc.initialise();
      cteCalc.addObservations(source, dest, other);
      cteResult = cteCalc.computeAverageLocalOfObservations();
    } catch (Throwable e) {
      exceptionCaught = true;
    }

    assertFalse(exceptionCaught);

		TransferEntropyCalculatorDiscrete teCalc =
				new TransferEntropyCalculatorDiscrete(2, 1);
		teCalc.initialise();
		teCalc.addObservations(source, dest);
		teResult = teCalc.computeAverageLocalOfObservations();

		assertEquals(cteResult, teResult, tol);
	}


  /**
   * Generate a time series such that dest = XOR(source, other)
   */
  public void testOthersSameBase() {

    int nb_samples = 10000;
    int[] source = rand.generateRandomInts(nb_samples, 2);
    int[] other = rand.generateRandomInts(nb_samples, 2);
    int[] dest = new int[nb_samples];
    boolean exceptionCaught = false;
    double result = Double.NaN;

    for (int i = 1; i < nb_samples; i++) {
      dest[i] = (source[i-1] + other[i-1]) % 2;
    }

    try {
      ConditionalTransferEntropyCalculatorDiscrete cteCalc =
          new ConditionalTransferEntropyCalculatorDiscrete(2, 1, 1);
      cteCalc.initialise();
      cteCalc.addObservations(source, dest, other);
      result = cteCalc.computeAverageLocalOfObservations();
    } catch (Throwable e) {
      exceptionCaught = true;
    }

    assertFalse(exceptionCaught);
    assertEquals(1, result, 1e-3);

  }

  /**
   * Go through a few examples that should/shouldn't throw exceptions.
   */
  public void testDifferentBases() {

    assertFalse(checkForException(2, 2, 1, 1, 100));
    assertFalse(checkForException(2, 3, 1, 1, 100));
    assertFalse(checkForException(2, 2, 3, 3, 100));
    assertFalse(checkForException(2, 3, 3, 3, 100));

    assertTrue(checkForException(2, 2, 3, 4, 100));
    assertTrue(checkForException(2, 2, 4, 3, 100));

  }

  /**
   * Basic usage of the calculator on random data.
   */
  private boolean checkForException(int base, int base_others, int col_others,
      int num_info_contributors, int rows) {

    boolean exceptionCaught = false;
    
    try {
      ConditionalTransferEntropyCalculatorDiscrete cteCalc =
          new ConditionalTransferEntropyCalculatorDiscrete(base, 1, num_info_contributors, base_others);
      int[][] other = rand.generateRandomInts(rows, col_others, base_others);
      int[] src = rand.generateRandomInts(rows, base);
      int[] tgt = rand.generateRandomInts(rows, base);
      cteCalc.initialise();
      cteCalc.addObservations(src, tgt, other);
      cteCalc.computeAverageLocalOfObservations();

    } catch (Throwable e) {
      exceptionCaught = true;
    }

    return exceptionCaught;

  }

}
