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

public class DualTotalCorrelationCalculatorGaussianTester extends TestCase {

  /**
   * For two variables, DTC is equal to mutual information.
   */
  public void testTwoVariables() throws Exception {
    double[][] cov = new double[][] {{1, 0.5}, {0.5, 1}};

    DualTotalCorrelationCalculatorGaussian dtcCalc = new DualTotalCorrelationCalculatorGaussian();
    dtcCalc.initialise(2);
    dtcCalc.setCovariance(cov);
    double dtc = dtcCalc.computeAverageLocalOfObservations();

    MutualInfoCalculatorMultiVariateGaussian miCalc = new MutualInfoCalculatorMultiVariateGaussian();
    miCalc.initialise(1,1);
    miCalc.setCovariance(cov, 1);
    double mi  = miCalc.computeAverageLocalOfObservations();

    assertEquals(dtc, mi, 1e-6);
  }

  /**
   * Compare against the direct calculation of DTC as a sum of entropies using the
   * entropy calculator.
   */
  public void testCompareWithEntropy() throws Exception {
    double[][] cov = new double[][] {{1, 0.4, 0.3}, {0.4, 1, 0.2}, {0.3, 0.2, 1}};

    DualTotalCorrelationCalculatorGaussian dtcCalc = new DualTotalCorrelationCalculatorGaussian();
    dtcCalc.initialise(3);
    dtcCalc.setCovariance(cov);
    double dtc = dtcCalc.computeAverageLocalOfObservations();

    // Calculate using an entropy calculator and picking submatrices manually
    EntropyCalculatorMultiVariateGaussian hCalc = new EntropyCalculatorMultiVariateGaussian();
    hCalc.initialise(3);
    hCalc.setCovariance(cov);
    double dtc_hCalc = -2 * hCalc.computeAverageLocalOfObservations();

    hCalc.initialise(2);
    hCalc.setCovariance(MatrixUtils.selectRowsAndColumns(cov, new int[] {0,1}, new int[] {0,1}));
    dtc_hCalc += hCalc.computeAverageLocalOfObservations();
    hCalc.initialise(2);
    hCalc.setCovariance(MatrixUtils.selectRowsAndColumns(cov, new int[] {0,2}, new int[] {0,2}));
    dtc_hCalc += hCalc.computeAverageLocalOfObservations();
    hCalc.initialise(2);
    hCalc.setCovariance(MatrixUtils.selectRowsAndColumns(cov, new int[] {1,2}, new int[] {1,2}));
    dtc_hCalc += hCalc.computeAverageLocalOfObservations();

    assertEquals(dtc, dtc_hCalc, 1e-6);

  }

  /**
   * Confirm that the local values average correctly back to the average value
   */
  public void testLocalsAverageCorrectly() throws Exception {

    int dimensions = 4;
    int timeSteps = 1000;
    DualTotalCorrelationCalculatorGaussian dtcCalc = new DualTotalCorrelationCalculatorGaussian();
    dtcCalc.initialise(dimensions);

    // generate some random data
    RandomGenerator rg = new RandomGenerator();
    double[][] data = rg.generateNormalData(timeSteps, dimensions,
        0, 1);

    dtcCalc.setObservations(data);

    double dtc = dtcCalc.computeAverageLocalOfObservations();
    double[] dtcLocal = dtcCalc.computeLocalOfPreviousObservations();

    System.out.printf("Average was %.5f%n", dtc);

    assertEquals(dtc, MatrixUtils.mean(dtcLocal), 0.00001);
  }

  /**
   * Confirm that for 2D the local values equal the local MI values
   *
   */
  public void testLocalsEqualMI() throws Exception {

    int dimensions = 2;
    int timeSteps = 100;
    DualTotalCorrelationCalculatorGaussian dtcCalc = new DualTotalCorrelationCalculatorGaussian();
    dtcCalc.initialise(dimensions);

    // generate some random data
    RandomGenerator rg = new RandomGenerator();
    double[][] data = rg.generateNormalData(timeSteps, dimensions,
        0, 1);

    dtcCalc.setObservations(data);
    double[] dtcLocal = dtcCalc.computeLocalOfPreviousObservations();

    MutualInfoCalculatorMultiVariateGaussian miCalc = new MutualInfoCalculatorMultiVariateGaussian();
    miCalc.initialise(1, 1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1));
    double[] miLocal = miCalc.computeLocalOfPreviousObservations();

    for (int t = 0; t < timeSteps; t++) {
      assertEquals(dtcLocal[t], miLocal[t], 0.00001);
    }
  }

}

