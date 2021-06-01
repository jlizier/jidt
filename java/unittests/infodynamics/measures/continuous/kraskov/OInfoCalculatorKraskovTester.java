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

import infodynamics.measures.continuous.gaussian.OInfoCalculatorGaussian;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class OInfoCalculatorKraskovTester extends TestCase {

  /**
   * Compare against the values obtained by the Gaussian O-info calculator when the data
   * is actually Gaussian.
   */
  public void testCompareWithGaussian() throws Exception {

    int N = 100000, D = 3;
		RandomGenerator rg = new RandomGenerator();
    rg.setSeed(3);
    double[][] data = rg.generateNormalData(N, D, 0, 1);
    for (int i = 0; i < N; i++) {
      data[i][2] = data[i][2] + 0.4*data[i][0];
      data[i][1] = data[i][1] + 0.5*data[i][0];
    }

		OInfoCalculatorKraskov oCalc_ksg = new OInfoCalculatorKraskov();
    oCalc_ksg.setProperty("NOISE_LEVEL_TO_ADD", "0");
    oCalc_ksg.initialise(3);
    oCalc_ksg.setObservations(data);
    double o_ksg = oCalc_ksg.computeAverageLocalOfObservations();

		OInfoCalculatorGaussian oCalc_gau = new OInfoCalculatorGaussian();
    oCalc_gau.initialise(3);
    oCalc_gau.setObservations(data);
    double o_gau = oCalc_gau.computeAverageLocalOfObservations();

    assertEquals(o_ksg, o_gau, 0.001);

  }

}


