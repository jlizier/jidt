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

package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.gaussian.DualTotalCorrelationCalculatorGaussian;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class DualTotalCorrelationCalculatorKraskovTester extends TestCase {

  /**
   * For two variables, DTC is equal to mutual information.
   */
	public void testTwoVariables() throws Exception {

    double[][] data;
		ArrayFileReader afr = new ArrayFileReader("demos/data/2randomCols-1.txt");
		data = afr.getDouble2DMatrix();

		DualTotalCorrelationCalculatorKraskov dtcCalc = new DualTotalCorrelationCalculatorKraskov();
    dtcCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    dtcCalc.initialise(2);
    dtcCalc.setObservations(data);
    double dtc = dtcCalc.computeAverageLocalOfObservations();


		MutualInfoCalculatorMultiVariateKraskov1 miCalc = new MutualInfoCalculatorMultiVariateKraskov1();
    miCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    miCalc.initialise(1,1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                           MatrixUtils.selectColumn(data, 1));
    double mi  = miCalc.computeAverageLocalOfObservations();

    assertEquals(dtc, mi, 1e-6);
	}

  /**
   * Compare against the values obtained by the Gaussian DTC calculator when the data
   * is actually Gaussian.
   */
  public void testCompareWithGaussian() throws Exception {

    int N = 100000, D = 3;
		RandomGenerator rg = new RandomGenerator();
    rg.setSeed(1);
    double[][] data = rg.generateNormalData(N, D, 0, 1);
    for (int i = 0; i < N; i++) {
      data[i][2] = data[i][2] + 0.1*data[i][0];
      data[i][1] = data[i][1] + 0.5*data[i][0];
    }

		DualTotalCorrelationCalculatorKraskov dtcCalc_ksg = new DualTotalCorrelationCalculatorKraskov();
    dtcCalc_ksg.setProperty("NOISE_LEVEL_TO_ADD", "0");
    dtcCalc_ksg.initialise(3);
    dtcCalc_ksg.setObservations(data);
    double dtc_ksg = dtcCalc_ksg.computeAverageLocalOfObservations();

		DualTotalCorrelationCalculatorGaussian dtcCalc_gau = new DualTotalCorrelationCalculatorGaussian();
    dtcCalc_gau.initialise(3);
    dtcCalc_gau.setObservations(data);
    double dtc_gau = dtcCalc_gau.computeAverageLocalOfObservations();

    assertEquals(dtc_ksg, dtc_gau, 0.001);

  }

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testLocalsAverageCorrectly() throws Exception {
		
    int dimensions = 4;
    int timeSteps = 1000;
		DualTotalCorrelationCalculatorKraskov dtcCalc = new DualTotalCorrelationCalculatorKraskov();
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
		DualTotalCorrelationCalculatorKraskov dtcCalc = new DualTotalCorrelationCalculatorKraskov();
		dtcCalc.initialise(dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] data = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		dtcCalc.setObservations(data);
		double[] dtcLocal = dtcCalc.computeLocalOfPreviousObservations();

    MutualInfoCalculatorMultiVariateKraskov1 miCalc = new MutualInfoCalculatorMultiVariateKraskov1();
    miCalc.initialise(1, 1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1));
    double[] miLocal = miCalc.computeLocalOfPreviousObservations();
		
		for (int t = 0; t < timeSteps; t++) {
      assertEquals(dtcLocal[t], miLocal[t], 0.00001);
    }
	}
	

}

