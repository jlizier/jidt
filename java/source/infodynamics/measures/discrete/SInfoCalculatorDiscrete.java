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

import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the S-information of a given multivariate
 *  <code>int[][]</code> set of
 *  observations (extending {@link MultiVariateInfoMeasureCalculatorDiscrete}).</p>
 *
 * <p>Usage is as per the paradigm outlined for {@link MultiVariateInfoMeasureCalculatorDiscrete}.
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>Rosas, F., Mediano, P., Gastpar, M, Jensen, H.,
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.100.032305">"Quantifying high-order
 *   interdependencies via multivariate extensions of the mutual information"</a>,
 *   Physical Review E 100, (2019) 032305.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class SInfoCalculatorDiscrete
    extends MultiVariateInfoMeasureCalculatorDiscrete {
  
  /**
   * Construct an instance.
   * 
   * @param base number of symbols for each variable.
   *        E.g. binary variables are in base-2.
   * @param numVars numbers of joint variables that DTC
   *     will be computed over.
   */
  public SInfoCalculatorDiscrete(int base, int numVars) {
    super(base, numVars);
  }

  protected double computeLocalValueForTuple(int[] tuple, int jointValue) {

    if (jointCount[jointValue] == 0) {
      // This joint state does not occur, so it makes no contribution here
      return 0;
    }

    double jointProb = (double) jointCount[jointValue] / (double) observations;

    // Local TC value
    double localTC = Math.log(jointProb);
    for (int i = 0; i < numVars; i++) {
      int marginalState = tuple[i];
      double marginalProb = (double) smallMarginalCounts[i][marginalState] / (double) observations;
      localTC -= Math.log(marginalProb);
    }

    // Local DTC value
    double localDTC = (numVars - 1) * Math.log(jointProb);
    for (int i = 0; i < numVars; i++) {
      int marginalState = computeBigMarginalState(jointValue, i, tuple[i]);
      double marginalProb = (double) bigMarginalCounts[i][marginalState] / (double) observations;
      localDTC -= Math.log(marginalProb);
    }

    // Combine local TC and DTC into S-info
    double logValue = localTC + localDTC;
    double localValue = logValue / log_2;

    if (jointProb > 0.0) {
      checkLocals(localValue);
    }

    return localValue;
  }

}

