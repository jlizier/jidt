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

package infodynamics.measures.discrete;

import infodynamics.utils.RandomGenerator;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.MathsUtils;

import java.util.Arrays;
import junit.framework.TestCase;


public class DualTotalCorrelationTester extends TestCase {

  public void testIndependent() throws Exception {
    DualTotalCorrelationCalculatorDiscrete dtcCalc = new DualTotalCorrelationCalculatorDiscrete(2, 3);
    double dtc = dtcCalc.compute(new int[][] {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}});
    assertEquals(0.0, dtc, 0.000001);
  }

  public void testXor() throws Exception {
    // 3 variables
    DualTotalCorrelationCalculatorDiscrete dtcCalc = new DualTotalCorrelationCalculatorDiscrete(2, 3);
    double dtc = dtcCalc.compute(new int[][] {{0,0,1},{0,1,0},{1,0,0},{1,1,1}});
    assertEquals(2.0, dtc, 0.000001);

    // 4 variables
    dtcCalc = new DualTotalCorrelationCalculatorDiscrete(2, 4);
    dtc = dtcCalc.compute(new int[][] {{0,0,0,1},{0,0,1,0},{0,1,0,0},{1,0,0,0},{1,1,1,0},{1,1,0,1},{1,0,1,1},{0,1,1,1}});
    assertEquals(3.0, dtc, 0.000001);
  }

  public void testCopy() throws Exception {
    // 3 variables
    DualTotalCorrelationCalculatorDiscrete dtcCalc = new DualTotalCorrelationCalculatorDiscrete(2, 3);
    double dtc = dtcCalc.compute(new int[][] {{0,0,0},{1,1,1}});
    assertEquals(1.0, dtc, 0.000001);

    // 4 variables
    dtcCalc = new DualTotalCorrelationCalculatorDiscrete(2, 4);
    dtc = dtcCalc.compute(new int[][] {{0,0,0,0},{1,1,1,1}});
    assertEquals(1.0, dtc, 0.000001);
  }

  public void testCompareEntropy() throws Exception {
    // Generate random data and check that it matches the explicit computation
    // using entropy calculators
    RandomGenerator rg = new RandomGenerator();
    int D = 4;
    int[][] data = rg.generateRandomInts(10, D, 2);

    // DTC calculator
    DualTotalCorrelationCalculatorDiscrete dtcCalc = new DualTotalCorrelationCalculatorDiscrete(2, D);
    double dtc_direct = dtcCalc.compute(data);

    // Entropy calculators
    EntropyCalculatorDiscrete hCalc = new EntropyCalculatorDiscrete(MathsUtils.power(2, D));
    hCalc.initialise();
    hCalc.addObservations(MatrixUtils.computeCombinedValues(data, 2));
    double dtc_test = (1 - D) * hCalc.computeAverageLocalOfObservations();

    hCalc = new EntropyCalculatorDiscrete(MathsUtils.power(2, D-1));
    for (int i = 0; i < 4; i++) {
      hCalc.initialise();
      hCalc.addObservations(MatrixUtils.computeCombinedValues(MatrixUtils.selectColumns(data, allExcept(i, D)), 2));
      dtc_test += hCalc.computeAverageLocalOfObservations();
    }

    assertEquals(dtc_direct, dtc_test, 0.000001);

  }

  protected int[] allExcept(int idx, int N) {
    boolean[] v = new boolean[N];
    Arrays.fill(v, true);
    v[idx] = false;

    int[] v2 = new int[N - 1];
    int counter = 0;
    for (int i = 0; i < N; i++) {
      if (v[i]) {
        v2[counter] = i;
        counter++;
      }
    }

    return v2;
  }


}

