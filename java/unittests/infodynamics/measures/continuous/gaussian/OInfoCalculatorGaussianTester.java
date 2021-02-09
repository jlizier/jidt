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

package infodynamics.measures.continuous.gaussian;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class OInfoCalculatorGaussianTester extends TestCase {

  /**
   * For two variables, O-info is zero
   */
  public void testTwoVariables() throws Exception {
    double[][] cov = new double[][] {{1, 0.5}, {0.5, 1}};

    OInfoCalculatorGaussian oCalc = new OInfoCalculatorGaussian();
    oCalc.initialise(2);
    oCalc.setCovariance(cov);
    double oinfo = oCalc.computeAverageLocalOfObservations();

    assertEquals(oinfo, 0, 1e-6);
  }


  /**
   * For factorisable pairwise interactions, O-info is zero
   */
  public void testPairwise() throws Exception {
    double[][] cov = new double[][] {{  1, 0.5,   0,   0},
                                     {0.5,   1,   0,   0},
                                     {  0,   0,   1, 0.5},
                                     {  0,   0, 0.5,   1}};

    OInfoCalculatorGaussian oCalc = new OInfoCalculatorGaussian();
    oCalc.initialise(4);
    oCalc.setCovariance(cov);
    double oinfo = oCalc.computeAverageLocalOfObservations();

    assertEquals(oinfo, 0, 1e-6);
  }

  /**
   * Compare against the direct calculation of O-info as a sum of entropies using the
   * entropy calculator.
   */
  public void testCompareWithEntropy() throws Exception {
    double[][] cov = new double[][] {{1, 0.4, 0.3}, {0.4, 1, 0.2}, {0.3, 0.2, 1}};

    OInfoCalculatorGaussian oCalc = new OInfoCalculatorGaussian();
    oCalc.initialise(3);
    oCalc.setCovariance(cov);
    double oinfo = oCalc.computeAverageLocalOfObservations();

    // Calculate using an entropy calculator and picking submatrices manually
    EntropyCalculatorMultiVariateGaussian hCalc = new EntropyCalculatorMultiVariateGaussian();
    hCalc.initialise(3);
    hCalc.setCovariance(cov);
    double oinfo_hCalc = hCalc.computeAverageLocalOfObservations();

    hCalc.initialise(2);
    hCalc.setCovariance(MatrixUtils.selectRowsAndColumns(cov, new int[] {0,1}, new int[] {0,1}));
    oinfo_hCalc -= hCalc.computeAverageLocalOfObservations();
    hCalc.initialise(2);
    hCalc.setCovariance(MatrixUtils.selectRowsAndColumns(cov, new int[] {0,2}, new int[] {0,2}));
    oinfo_hCalc -= hCalc.computeAverageLocalOfObservations();
    hCalc.initialise(2);
    hCalc.setCovariance(MatrixUtils.selectRowsAndColumns(cov, new int[] {1,2}, new int[] {1,2}));
    oinfo_hCalc -= hCalc.computeAverageLocalOfObservations();

    hCalc.initialise(1);
    hCalc.setCovariance(new double[][] {{cov[0][0]}});
    oinfo_hCalc += hCalc.computeAverageLocalOfObservations();
    hCalc.initialise(1);
    hCalc.setCovariance(new double[][] {{cov[1][1]}});
    oinfo_hCalc += hCalc.computeAverageLocalOfObservations();
    hCalc.initialise(1);
    hCalc.setCovariance(new double[][] {{cov[2][2]}});
    oinfo_hCalc += hCalc.computeAverageLocalOfObservations();

    assertEquals(oinfo, oinfo_hCalc, 1e-6);

  }

  /**
   * Confirm that the local values average correctly back to the average value
   */
  public void testLocalsAverageCorrectly() throws Exception {

    int dimensions = 4;
    int timeSteps = 1000;
    OInfoCalculatorGaussian oCalc = new OInfoCalculatorGaussian();
    oCalc.initialise(dimensions);

    // generate some random data
    RandomGenerator rg = new RandomGenerator();
    double[][] data = rg.generateNormalData(timeSteps, dimensions,
        0, 1);

    oCalc.setObservations(data);

    double oinfo = oCalc.computeAverageLocalOfObservations();
    double[] oLocal = oCalc.computeLocalOfPreviousObservations();

    System.out.printf("Average was %.5f%n", oinfo);

    assertEquals(oinfo, MatrixUtils.mean(oLocal), 0.00001);
  }
  
  /**
   * Confirm that for 2D all local values equal zero
   */
  public void testLocalsEqualZero() throws Exception {

    int dimensions = 2;
    int timeSteps = 100;
    OInfoCalculatorGaussian oCalc = new OInfoCalculatorGaussian();
    oCalc.initialise(dimensions);

    // generate some random data
    RandomGenerator rg = new RandomGenerator();
    double[][] data = rg.generateNormalData(timeSteps, dimensions,
        0, 1);

    oCalc.setObservations(data);
    double[] oLocal = oCalc.computeLocalOfPreviousObservations();

    for (int t = 0; t < timeSteps; t++) {
      assertEquals(oLocal[t], 0, 0.00001);
    }
  }

}

