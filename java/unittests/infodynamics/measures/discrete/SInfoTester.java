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

public class SInfoTester extends TestCase {

  public void testIndependent() throws Exception {
    SInfoCalculatorDiscrete sCalc = new SInfoCalculatorDiscrete(2, 3);
    double sinfo = sCalc.compute(new int[][] {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}});
    assertEquals(0.0, sinfo, 0.000001);
  }

  public void testXor() throws Exception {
    // 3 variables
    SInfoCalculatorDiscrete sCalc = new SInfoCalculatorDiscrete(2, 3);
    double sinfo = sCalc.compute(new int[][] {{0,0,1},{0,1,0},{1,0,0},{1,1,1}});
    assertEquals(3.0, sinfo, 0.000001);

    // 4 variables
    sCalc = new SInfoCalculatorDiscrete(2, 4);
    sinfo = sCalc.compute(new int[][] {{0,0,0,1},{0,0,1,0},{0,1,0,0},{1,0,0,0},{1,1,1,0},{1,1,0,1},{1,0,1,1},{0,1,1,1}});
    assertEquals(4.0, sinfo, 0.000001);
  }

  public void testCopy() throws Exception {
    // 3 variables
    SInfoCalculatorDiscrete sCalc = new SInfoCalculatorDiscrete(2, 3);
    double sinfo = sCalc.compute(new int[][] {{0,0,0},{1,1,1}});
    assertEquals(3.0, sinfo, 0.000001);

    // 4 variables
    sCalc = new SInfoCalculatorDiscrete(2, 4);
    sinfo = sCalc.compute(new int[][] {{0,0,0,0},{1,1,1,1}});
    assertEquals(4.0, sinfo, 0.000001);
  }

  public void testCompareTCAndDTC() throws Exception {
    // Generate random data and check that it matches the explicit computation
    // using TC and DTC calculators
    RandomGenerator rg = new RandomGenerator();
    int[][] data = rg.generateRandomInts(10, 4, 2);

    // O-info calculator
    SInfoCalculatorDiscrete sCalc = new SInfoCalculatorDiscrete(2, 4);
    double sinfo_direct = sCalc.compute(data);

    // TC and DTC calculators
    MultiInformationCalculatorDiscrete tcCalc = new MultiInformationCalculatorDiscrete(2, 4);
    DualTotalCorrelationCalculatorDiscrete dtcCalc = new DualTotalCorrelationCalculatorDiscrete(2, 4);
    tcCalc.initialise();
    tcCalc.addObservations(data);
    double sinfo_test = tcCalc.computeAverageLocalOfObservations() + dtcCalc.compute(data);

    assertEquals(sinfo_direct, sinfo_test, 0.000001);

  }

}

