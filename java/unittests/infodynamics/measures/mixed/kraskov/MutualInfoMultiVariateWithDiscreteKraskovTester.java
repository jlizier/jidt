/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.measures.mixed.kraskov;

import junit.framework.TestCase;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.MathsUtils;

public class MutualInfoMultiVariateWithDiscreteKraskovTester extends TestCase {

  int DATA_LENGTH = 100;
  int TEST_DIMENSIONS = 2;
  int TEST_BASE = 2;

  public void testSetObservationsDataIntegrity() throws Exception { 
    
    MutualInfoCalculatorMultiVariateWithDiscreteKraskov miCalc =
        new MutualInfoCalculatorMultiVariateWithDiscreteKraskov();
    miCalc.initialise(TEST_DIMENSIONS, TEST_BASE);
    
    // Generate proper and incorrect data sets
    RandomGenerator rg = new RandomGenerator();
    double[][] contDataTooManyVariables = rg.generateNormalData(DATA_LENGTH,
        TEST_DIMENSIONS + 1, 0, 1);
    double[][] contDataCorrect = rg.generateNormalData(DATA_LENGTH, TEST_DIMENSIONS,
        0, 1);
    int[] discDataCorrect = rg.generateRandomInts(DATA_LENGTH, TEST_BASE);
    int[] discDataValuesOutsideRange = rg.generateRandomInts(100, TEST_BASE + 1);
    int[] discDataWrongLength = rg.generateRandomInts(DATA_LENGTH - 1, TEST_BASE);

    // Check that we catch continuous data which doesn't have
    //  the right number of variables
    boolean caughtException = false;
    try {
      miCalc.setObservations(contDataTooManyVariables, discDataCorrect);
    } catch (Exception e) {
      caughtException = true;
    }
    assertTrue(caughtException);
    
    // Check that we catch discrete data outside the range:
    caughtException = false;
    try {
      miCalc.setObservations(contDataCorrect, discDataValuesOutsideRange);
    } catch (Exception e) {
      caughtException = true;
    }
    assertTrue(caughtException);
    
    // Check that we catch mismatched data lengths:
    caughtException = false;
    try {
      miCalc.setObservations(contDataCorrect, discDataWrongLength);
    } catch (Exception e) {
      caughtException = true;
    }
    assertTrue(caughtException);

    // No observations should have been set yet
    assertEquals(0, miCalc.getNumObservations());
    
    // Check that no exception is thrown for ok data:
    caughtException = false;
    try {
      miCalc.setObservations(contDataCorrect, discDataCorrect);
    } catch (Exception e) {
      caughtException = true;
      e.printStackTrace();
    }
    assertFalse(caughtException);
    
    // and that the observations have now been set:
    assertEquals(DATA_LENGTH, miCalc.getNumObservations());
  }
    
  public void testComputeSignificanceDoesntAlterAverage() throws Exception {
    
    MutualInfoCalculatorMultiVariateWithDiscreteKraskov miCalc =
        new MutualInfoCalculatorMultiVariateWithDiscreteKraskov();
    miCalc.initialise(TEST_DIMENSIONS, TEST_BASE);
    
    // Generate proper data sets
    RandomGenerator rg = new RandomGenerator();
    double[][] contDataCorrect = rg.generateNormalData(DATA_LENGTH, TEST_DIMENSIONS,
        0, 1);
    int[] discDataCorrect = rg.generateRandomInts(DATA_LENGTH, TEST_BASE);

    miCalc.setObservations(contDataCorrect, discDataCorrect);

    double mi = miCalc.computeAverageLocalOfObservations();
    
    System.out.printf("Average was %.5f\n", mi);
    
    EmpiricalMeasurementDistribution measDist =
        miCalc.computeSignificance(100);
    
    System.out.printf("pValue of sig test was %.3f\n", measDist.pValue);
    
    // And compute the average value again to check that it's consistent:
    for (int i = 0; i < 10; i++) {
      double averageCheck1 = miCalc.computeAverageLocalOfObservations();
      assertEquals(mi, averageCheck1);
    }
  }
  
  public void testObservationsRequiredBeforeStatTest() throws Exception {
    MutualInfoCalculatorMultiVariateWithDiscreteKraskov miCalc =
        new MutualInfoCalculatorMultiVariateWithDiscreteKraskov();
    miCalc.initialise(TEST_DIMENSIONS, TEST_BASE);

    boolean caughtException = false;
    try {
      miCalc.computeSignificance(100);
    } catch (Exception e) {
      caughtException = true;
    }
    assertTrue(caughtException);
  }

  public void testLocals() throws Exception {
    MutualInfoCalculatorMultiVariateWithDiscreteKraskov miCalc =
        new MutualInfoCalculatorMultiVariateWithDiscreteKraskov();

    // Generate proper data sets
    RandomGenerator rg = new RandomGenerator();
    double[][] contData = rg.generateNormalData(DATA_LENGTH, TEST_DIMENSIONS,
        0, 1);
    int[] discData = rg.generateRandomInts(DATA_LENGTH, TEST_BASE);

    double[][] contDataNew = rg.generateNormalData(10, TEST_DIMENSIONS, 0, 1);
    int[] discDataNew = rg.generateRandomInts(10, TEST_BASE);

    // Check that result of locals does not depend on other points in the same query
    miCalc.initialise(TEST_DIMENSIONS, TEST_BASE);
    miCalc.setObservations(contData, discData);
    double[] locals1 = miCalc.computeLocalUsingPreviousObservations(MatrixUtils.selectRows(contDataNew, 0, 5), MatrixUtils.select(discDataNew, 0, 5));
    double[] locals2 = miCalc.computeLocalUsingPreviousObservations(contDataNew, discDataNew);
    for (int i = 0; i < 5; i++) {
      assertEquals(locals1[i], locals2[i]);
    }

    // Check that computeAverage and average of computeLocal are the same
    miCalc.initialise(TEST_DIMENSIONS, TEST_BASE);
    miCalc.setObservations(contData, discData);
    double average = miCalc.computeAverageLocalOfObservations();

    miCalc.initialise(TEST_DIMENSIONS, TEST_BASE);
    miCalc.setObservations(contData, discData);
    double averageOfLocals = MatrixUtils.mean(miCalc.computeLocalOfPreviousObservations());

    assertEquals(average, averageOfLocals, 0.01);
    
  }

  public void testCompareAnalyticalValue() throws Exception {
    MutualInfoCalculatorMultiVariateWithDiscreteKraskov miCalc =
        new MutualInfoCalculatorMultiVariateWithDiscreteKraskov();

    double[] separation = {2.0, 4.0, 8.0};

    // These values were computed using the test script in the supplementary
    // material of
    // B. Ross, "Mutual Information Between Discrete and Continuous Data Sets",
    // PLOS ONE (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0087357)
    double[] true_val = {0.3368, 0.6327, 0.6931};
    int N = 10000;

    // Generate data with unit Gaussians separated by fixed distance
    RandomGenerator rg = new RandomGenerator();
    for (int j = 0; j < separation.length; j++) {
      double[][] contData = rg.generateNormalData(N, 1, 0, 1);
      int[] discData = rg.generateRandomInts(N, 2);
      for (int i = 0; i < N; i++) {
        if (discData[i] > 0) {
          contData[i][0] += separation[j];
        }
      }

      miCalc.initialise(1, 2);
      miCalc.setObservations(contData, discData);
      double res = miCalc.computeAverageLocalOfObservations();
      assertEquals(true_val[j], res, 0.05);
    }

    // Using the last dataset added, compare a few local
    // values with the analytical calculation
    int nb_locals_test = 10;
    double delta = separation[separation.length-1];
    for (int i = 0; i < nb_locals_test; i++) {
      int y = rg.generateRandomInts(1, 2)[0];
      double x = rg.generateNormalData(1, 0, 1)[0] + delta*y;
      double analytical = Math.log(MathsUtils.normalPdf(x, delta*y, 1)) - Math.log(0.5*(MathsUtils.normalPdf(x, 0, 1) + MathsUtils.normalPdf(x, delta, 1)));
      double estimated = miCalc.computeLocalUsingPreviousObservations(new double[][] {new double[] {x} }, new int[] {y})[0];
      assertEquals(analytical, estimated, 0.05);
    }
  }

  public void testProperties() throws Exception {
    MutualInfoCalculatorMultiVariateWithDiscreteKraskov miCalc =
        new MutualInfoCalculatorMultiVariateWithDiscreteKraskov();

    // Generate data with unit Gaussians separated by fixed distance
    double separation = 4.0;
    double true_val = 0.6327;
    int N = 1000;
    RandomGenerator rg = new RandomGenerator();
    double[][] contData = rg.generateNormalData(N, 1, 0, 1);
    int[] discData = rg.generateRandomInts(N, 2);
    for (int i = 0; i < N; i++) {
      if (discData[i] > 0) {
        contData[i][0] += separation;
      }
    }

    // Test we get correct result with default properties
    miCalc.initialise(1, 2);
    miCalc.setObservations(contData, discData);
    double res1 = miCalc.computeAverageLocalOfObservations();
    assertEquals(true_val, res1, 0.05);

    // Test that timeDiff has an effect
    miCalc.setProperty("TIME_DIFF", "1");
    miCalc.initialise(1, 2);
    miCalc.setObservations(contData, discData);
    double res2 = miCalc.computeAverageLocalOfObservations();
    assertEquals(0.0, res2, 0.05);
    miCalc.setProperty("TIME_DIFF", "0");

    // Test that K has an effect
    miCalc.setProperty("K", "2");
    miCalc.initialise(1, 2);
    miCalc.setObservations(contData, discData);
    double res3 = miCalc.computeAverageLocalOfObservations();
    assertTrue(Math.abs(res3 - res1) > 0.001);
    miCalc.setProperty("k", "4");

    // Test that noiseLevel has an effect
    miCalc.setProperty("NOISE_LEVEL_TO_ADD", "20.0");
    miCalc.initialise(1, 2);
    miCalc.setObservations(contData, discData);
    double res4 = miCalc.computeAverageLocalOfObservations();
    assertEquals(0.0, res4, 0.05);
    miCalc.setProperty("NOISE_LEVEL_TO_ADD", "false");
    
    // TEST REMOVED. dynCorrExcl not properly implemented yet
    // // Test that dynCorrExcl has an effect
    // miCalc.setProperty("DYN_CORR_EXCL", "10");
    // miCalc.initialise(1, 2);
    // miCalc.setObservations(contData, discData);
    // double res5 = miCalc.computeAverageLocalOfObservations();
    // assertTrue(Math.abs(res5 - res1) > 0.001);

  }

  public void testUnivariateOverloadings() throws Exception {
    MutualInfoCalculatorMultiVariateWithDiscreteKraskov miCalc =
        new MutualInfoCalculatorMultiVariateWithDiscreteKraskov();

    // Generate data sets
    RandomGenerator rg = new RandomGenerator();
    double[][] contDataMatrix = rg.generateRandomData(100, 2);
    double[] contDataVector   = rg.generateRandomData(100);
    int[] discData = rg.generateRandomInts(100, 2);

    // Check that no exception is thrown for ok data:
    boolean caughtException = false;
    try {
      miCalc.initialise(1, 2);
      miCalc.setObservations(contDataVector, discData);
    } catch (Exception e) {
      caughtException = true;
    }
    assertFalse(caughtException);

    // Check that exception is thrown when dimensions do not match
    caughtException = false;
    try {
      miCalc.initialise(1, 2);
      miCalc.setObservations(contDataMatrix, discData);
    } catch (Exception e) {
      caughtException = true;
    }
    assertTrue(caughtException);

    caughtException = false;
    try {
      miCalc.initialise(2, 2);
      miCalc.setObservations(contDataVector, discData);
    } catch (Exception e) {
      caughtException = true;
    }
    assertTrue(caughtException);


  }
}
