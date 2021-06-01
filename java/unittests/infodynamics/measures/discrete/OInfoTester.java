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

import junit.framework.TestCase;

public class OInfoTester extends TestCase {

  public void testIndependent() throws Exception {
    OInfoCalculatorDiscrete oCalc = new OInfoCalculatorDiscrete(2, 3);
    double oinfo = oCalc.compute(new int[][] {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}});
    assertEquals(0.0, oinfo, 0.000001);
  }

  public void testXor() throws Exception {
    // 3 variables
    OInfoCalculatorDiscrete oCalc = new OInfoCalculatorDiscrete(2, 3);
    double oinfo = oCalc.compute(new int[][] {{0,0,1},{0,1,0},{1,0,0},{1,1,1}});
    assertEquals(-1.0, oinfo, 0.000001);

    // 4 variables
    oCalc = new OInfoCalculatorDiscrete(2, 4);
    oinfo = oCalc.compute(new int[][] {{0,0,0,1},{0,0,1,0},{0,1,0,0},{1,0,0,0},{1,1,1,0},{1,1,0,1},{1,0,1,1},{0,1,1,1}});
    assertEquals(-2.0, oinfo, 0.000001);
  }

  public void testCopy() throws Exception {
    // 3 variables
    OInfoCalculatorDiscrete oCalc = new OInfoCalculatorDiscrete(2, 3);
    double oinfo = oCalc.compute(new int[][] {{0,0,0},{1,1,1}});
    assertEquals(1.0, oinfo, 0.000001);

    // 4 variables
    oCalc = new OInfoCalculatorDiscrete(2, 4);
    oinfo = oCalc.compute(new int[][] {{0,0,0,0},{1,1,1,1}});
    assertEquals(2.0, oinfo, 0.000001);
  }

  public void testPairwise() throws Exception {
    // Variables 0 and 1 are correlated and independent from 2 and 3, that are
    // also correlated
    OInfoCalculatorDiscrete oCalc = new OInfoCalculatorDiscrete(2, 4);
    double oinfo = oCalc.compute(new int[][] {{0,0,0,0},{0,0,1,1},{1,1,0,0},{1,1,1,1}});
    assertEquals(0.0, oinfo, 0.000001);
  }

  public void testCompareTCAndDTC() throws Exception {
    // Generate random data and check that it matches the explicit computation
    // using TC and DTC calculators
    RandomGenerator rg = new RandomGenerator();
    int[][] data = rg.generateRandomInts(10, 4, 2);

    // O-info calculator
    OInfoCalculatorDiscrete oCalc = new OInfoCalculatorDiscrete(2, 4);
    double oinfo_direct = oCalc.compute(data);

    // TC and DTC calculators
    MultiInformationCalculatorDiscrete tcCalc = new MultiInformationCalculatorDiscrete(2, 4);
    DualTotalCorrelationCalculatorDiscrete dtcCalc = new DualTotalCorrelationCalculatorDiscrete(2, 4);
    tcCalc.initialise();
    tcCalc.addObservations(data);
    double oinfo_test = tcCalc.computeAverageLocalOfObservations() - dtcCalc.compute(data);

    assertEquals(oinfo_direct, oinfo_test, 0.000001);

  }

}

