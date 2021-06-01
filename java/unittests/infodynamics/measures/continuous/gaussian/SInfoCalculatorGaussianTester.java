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

public class SInfoCalculatorGaussianTester extends TestCase {

  /**
   * For two variables, S-info is twice the MI between them
   */
  public void testTwoVariables() throws Exception {
    double[][] cov = new double[][] {{1, 0.5}, {0.5, 1}};

    SInfoCalculatorGaussian sCalc = new SInfoCalculatorGaussian();
    sCalc.initialise(2);
    sCalc.setCovariance(cov);
    double sinfo = sCalc.computeAverageLocalOfObservations();

    MutualInfoCalculatorMultiVariateGaussian miCalc = new MutualInfoCalculatorMultiVariateGaussian();
    miCalc.initialise(1,1);
    miCalc.setCovariance(cov, false);
    double mi = miCalc.computeAverageLocalOfObservations();

    assertEquals(sinfo, 2*mi, 1e-6);
  }


  /**
   * Compare against the direct calculation of S-info as a sum of mutual informations.
   */
  public void testCompareWithEntropy() throws Exception {
    double[][] cov = new double[][] {{1, 0.4, 0.3}, {0.4, 1, 0.2}, {0.3, 0.2, 1}};

    SInfoCalculatorGaussian sCalc = new SInfoCalculatorGaussian();
    sCalc.initialise(3);
    sCalc.setCovariance(cov);
    double sinfo = sCalc.computeAverageLocalOfObservations();

    // Calculate using a mutual info calculator and picking submatrices manually
    MutualInfoCalculatorMultiVariateGaussian miCalc = new MutualInfoCalculatorMultiVariateGaussian();
    miCalc.initialise(2,1);
    miCalc.setCovariance(cov, false);
    double sinfo_miCalc = miCalc.computeAverageLocalOfObservations();

    miCalc.initialise(2,1);
    miCalc.setCovariance(MatrixUtils.selectRowsAndColumns(cov, new int[] {2,0,1}, new int[] {2,0,1}), false);
    sinfo_miCalc += miCalc.computeAverageLocalOfObservations();
    miCalc.initialise(2,1);
    miCalc.setCovariance(MatrixUtils.selectRowsAndColumns(cov, new int[] {1,2,0}, new int[] {1,2,0}), false);
    sinfo_miCalc += miCalc.computeAverageLocalOfObservations();

    assertEquals(sinfo, sinfo_miCalc, 1e-6);

  }

  /**
   * Confirm that the local values average correctly back to the average value
   */
  public void testLocalsAverageCorrectly() throws Exception {

    int dimensions = 4;
    int timeSteps = 1000;
    SInfoCalculatorGaussian sCalc = new SInfoCalculatorGaussian();
    sCalc.initialise(dimensions);

    // generate some random data
    RandomGenerator rg = new RandomGenerator();
    double[][] data = rg.generateNormalData(timeSteps, dimensions,
        0, 1);

    sCalc.setObservations(data);

    double sinfo = sCalc.computeAverageLocalOfObservations();
    double[] sLocal = sCalc.computeLocalOfPreviousObservations();

    System.out.printf("Average was %.5f%n", sinfo);

    assertEquals(sinfo, MatrixUtils.mean(sLocal), 0.00001);
  }

  /**
   * Confirm that for 2D all local values equal zero
   */
  public void testLocalsEqualMI() throws Exception {

    int dimensions = 2;
    int timeSteps = 100;
    SInfoCalculatorGaussian sCalc = new SInfoCalculatorGaussian();
    sCalc.initialise(dimensions);

    // generate some random data
    RandomGenerator rg = new RandomGenerator();
    double[][] data = rg.generateNormalData(timeSteps, dimensions,
        0, 1);

    sCalc.setObservations(data);
    double[] sLocal = sCalc.computeLocalOfPreviousObservations();

    MutualInfoCalculatorMultiVariateGaussian miCalc = new MutualInfoCalculatorMultiVariateGaussian();
    miCalc.initialise(1, 1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1));
    double[] miLocal = miCalc.computeLocalOfPreviousObservations();

    for (int t = 0; t < timeSteps; t++) {
      assertEquals(sLocal[t], 2*miLocal[t], 0.00001);
    }
  }

}

