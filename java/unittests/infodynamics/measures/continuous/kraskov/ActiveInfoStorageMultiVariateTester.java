package infodynamics.measures.continuous.kraskov;

import java.util.Arrays;

import infodynamics.measures.continuous.ActiveInfoStorageCalculatorMultiVariate;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

/**
 * @author Pedro AM Mediano (<a href="pmediano at imperial.ac.uk">email</a>,
 * <a href="https://www.doc.ic.ac.uk/~pam213/">www</a>)
 */
public class ActiveInfoStorageMultiVariateTester extends TestCase {

	protected String NUM_THREADS_TO_USE_DEFAULT = MutualInfoCalculatorMultiVariateKraskov.USE_ALL_THREADS;
	protected String NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;
	
	public void testEmbedding() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();

		int[] ks = {1,2,3,4};
		int[] taus = {1,2,3,4};
		
		for (int kIndex = 0; kIndex < ks.length; kIndex++) {
			int k = ks[kIndex];
			for (int tauIndex = 0; tauIndex < taus.length; tauIndex++) {
				int tau = taus[tauIndex];
		
				ActiveInfoStorageCalculatorMultiVariateKraskov ais =
          new ActiveInfoStorageCalculatorMultiVariateKraskov();
				ais.initialise(data[0].length, k, tau);
				ais.setObservations(data);

				assertEquals(data.length -  (k-1)*tau - 1,
						ais.getNumObservations());
			}
		}
	}

  public void testUnivariateMethods() throws Exception {

    ActiveInfoStorageCalculatorMultiVariateKraskov ais =
      new ActiveInfoStorageCalculatorMultiVariateKraskov();

    boolean caughtException = false;

    double[] oneDArray = new double[] {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double[][] twoDArray = new double[][] {{1, 1.5}, {2, 2.5}, {3, 3.5}, {4, 4.5}, {5, 5.5}, {6, 6.5}, {7, 7.5}, {8, 8.5} ,{9, 9.5}};

    // Add 1D observations to a 1D calculator. Should not throw exception.
    try {
      ais.initialise(1);
      ais.startAddObservations();
      ais.addObservations(oneDArray);
      ais.finaliseAddObservations();
      ais.setObservations(oneDArray);
    } catch (Throwable e) {
      caughtException = true;
    }
    assertFalse(caughtException);

    // Add 2D observations to a 2D calculator. Should not throw exception.
    try {
      ais.initialise(2);
      ais.startAddObservations();
      ais.addObservations(twoDArray);
      ais.finaliseAddObservations();
      ais.setObservations(twoDArray);
    } catch (Throwable e) {
      caughtException = true;
    }
    assertFalse(caughtException);

    // Add 1D observations to a 2D calculator. Should throw exception.
    try {
      ais.initialise(2);
      ais.startAddObservations();
      ais.addObservations(oneDArray);
      ais.finaliseAddObservations();
      ais.setObservations(oneDArray);
    } catch (Throwable e) {
      caughtException = true;
    }
    assertTrue(caughtException);
    caughtException = false;

    // Add 2D observations to a 1D calculator. Should throw exception.
    try {
      ais.initialise(1);
      ais.startAddObservations();
      ais.addObservations(twoDArray);
      ais.finaliseAddObservations();
      ais.setObservations(twoDArray);
    } catch (Throwable e) {
      caughtException = true;
    }
    assertTrue(caughtException);
    caughtException = false;

    // Add a vector with too few observations for the given k, tau. Should not
    // throw exception.
    try {
      double[][] shortArray = new double[][] {{1, 2}, {3, 4}};
      int numObs_before = ais.getNumObservations();
      ais.initialise(2, 2, 2);
      ais.startAddObservations();
      ais.addObservations(shortArray);
      assertEquals(numObs_before, ais.getNumObservations());
    } catch (Throwable e) {
      caughtException = true;
    }
    assertFalse(caughtException);

    // Set a vector with too few observations for the given k, tau. Should
    // throw exception.
    try {
      double[][] shortArray = new double[][] {{1, 2}, {3, 4}};
      ais.initialise(2, 2, 2);
      ais.setObservations(shortArray);
    } catch (Throwable e) {
      caughtException = true;
    }
    assertTrue(caughtException);

  }

  public void testSyntheticData() throws Exception {
    // TODO: test something about the value that comes out

    System.out.println("Generate data for AIS multivariate Kraskov test");
    int T = 5500, discard = 500;
    int dim = 2;
    double[][] data = new double[T][dim];
    double[][] A = new double[][] {{0.4, 0.4}, {0.4, 0.4}};
    RandomGenerator rg = new RandomGenerator();
    // TODO: This is an AR process, so we could calculate the true stationary
    // covariance in Matlab and test we get the same result here
    for (int t = 1; t < T; t++) {
      data[t] = MatrixUtils.add(MatrixUtils.matrixProduct(A, data[t-1]),
                                rg.generateNormalData(dim, 0, 1));
    }
    data = MatrixUtils.selectRows(data, 0, T - discard);

    ActiveInfoStorageCalculatorMultiVariateKraskov ais =
      new ActiveInfoStorageCalculatorMultiVariateKraskov();

    boolean caughtException = false;
		int[] ks = {1,2,3,4};
		int[] taus = {1,2,3,4};

		for (int kIndex = 0; kIndex < ks.length; kIndex++) {
			int k = ks[kIndex];
			for (int tauIndex = 0; tauIndex < taus.length; tauIndex++) {
				int tau = taus[tauIndex];

        try {
          ais.initialise(dim, k, tau);
          ais.setObservations(data);
          ais.computeAverageLocalOfObservations();
        } catch (Throwable e) {
          caughtException = true;
        }

        assertFalse(caughtException);

      }
    }

  }

  public void testChangeAlgorithms() throws Exception {

    boolean caughtException = false;

    try {
      // Test everything is ok when changing KSG algorithms
      ActiveInfoStorageCalculatorMultiVariateKraskov ais =
        new ActiveInfoStorageCalculatorMultiVariateKraskov();
      ais.initialise(1);
      ais.setProperty("ALG_NUM", "1");
      ais.initialise(1);
    } catch (Throwable e) {
      caughtException = true;
    }
    assertFalse(caughtException);

    try {
      // Instantiate calculator with algorithm 1 from the start
      ActiveInfoStorageCalculatorMultiVariateKraskov ais =
        new ActiveInfoStorageCalculatorMultiVariateKraskov(1);
      ais.initialise(1);
      ais.setProperty("ALG_NUM", "2");
      ais.initialise(1);
    } catch (Throwable e) {
      caughtException = true;
    }
    assertFalse(caughtException);

    // Test that after changing algorithms the properties remain
    String newKNN = "6";
    ActiveInfoStorageCalculatorMultiVariateKraskov ais =
      new ActiveInfoStorageCalculatorMultiVariateKraskov();
    ais.setProperty("K", newKNN);
    ais.initialise(2);
    ais.setProperty("ALG_NUM", "1");
    ais.initialise(2);
    assertEquals(newKNN, ais.getProperty("K"));
    
  }

  public void testUnivariateComputations() throws Exception {

		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();

    ActiveInfoStorageCalculatorKraskov ais1 =
      new ActiveInfoStorageCalculatorKraskov();

    ActiveInfoStorageCalculatorMultiVariateKraskov ais2 =
      new ActiveInfoStorageCalculatorMultiVariateKraskov();

    ais1.setProperty("NOISE_TO_ADD", "0");
    ais1.initialise();
    ais1.setObservations(MatrixUtils.selectColumn(data, 0));
    double ais1Result = ais1.computeAverageLocalOfObservations();

    ais2.setProperty("NOISE_TO_ADD", "0");
    ais2.initialise(1);
    ais2.setObservations(MatrixUtils.selectColumn(data, 0));
    double ais2Result = ais2.computeAverageLocalOfObservations();

    assertEquals(ais1Result, ais2Result, 0.000001);

  }

}

