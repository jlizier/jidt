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
 * Implements a base class with common functionality for child class
 * implementations of multivariate information measures. 
 *
 * <p>Multivariate information measures are functionals of probability
 * distributions over <code>R^n</code>, and typical examples include multi-information
 * (a.k.a. total correlation), dual total correlation, O-information, and connected
 * information.</p>
 * 
 * <p>Usage of child classes is intended to follow this paradigm:</p>
 * <ol>
 *    <li>Construct the calculator;</li>
 *    <li>Initialise the calculator using {@link #initialise()};</li>
 *    <li>Provide the observations/samples for the calculator
 *        to set up the PDFs, using one or more calls to
 *      sets of {@link #addObservations(int[][], int[])} methods, then</li>
 *    <li>Compute the required quantities, being one or more of:
 *      <ul>
 *        <li>the average measure: {@link #computeAverageLocalOfObservations()};</li>
 *      </ul>
 *    </li>
 *    <li>
 *    Return to step 2 to re-use the calculator on a new data set.
 *    </li>
 *  </ol>
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
public abstract class MultiVariateInfoMeasureCalculatorDiscrete
    extends InfoMeasureCalculatorDiscrete {
  
  /**
   * Count of occurrences of each joint state in the provided observations.
   */
  protected int[] jointCount = null;

  /**
   * Count of occurrences of each state of each variable in the provided
   * observations.
   *
   * For a given variable <code>v</code> and state <code>i</code>,
   * <code>smallMarginalCounts[v][i]</code> counts how many times variable
   * <code>v</code> was observed in state <code>i</code>,
   */
  protected int[][] smallMarginalCounts = null; // marginalCounts[marginalIndex][state]

  /**
   * Count of occurrences of each state of each (D-1)-dimensional marginal in
   * the provided observations.
   *
   * For a given variable <code>v</code> and state <code>i</code>,
   * <code>bigMarginalCounts[v][i]</code> counts how many times the _rest_ of
   * the system, excluding variable <code>v</code>, was observed in state <code>i</code>,
   */
  protected int[][] bigMarginalCounts = null;
  
  /**
   * Number of variables in the system.
   */
  protected int numVars;

  /**
   * Number of possible states of the whole system.
   */
  protected int jointStates;

  /**
   * Whether the first local value has been checked. (Used to initialise some variables
   * related to computation of local values.
   */
  protected boolean checkedFirst = false;

  /**
   * Abstract constructor (to be called by child classes).
   * 
   * @param base number of symbols for each variable.
   *        E.g. binary variables are in base-2.
   * @param numVars numbers of joint variables that the measure
   *     will be computed over.
   */
  protected MultiVariateInfoMeasureCalculatorDiscrete(int base, int numVars) {
    super(base);
    this.numVars = numVars;
    jointStates = MathsUtils.power(base, numVars);
    try {
      jointCount = new int[jointStates];
      smallMarginalCounts = new int[numVars][base];
      bigMarginalCounts = new int[numVars][jointStates];
    } catch (OutOfMemoryError e) {
      // Allow any Exceptions to be thrown, but catch and wrap
      //  Error as a RuntimeException
      throw new RuntimeException("Requested memory for the base " +
          base + " with " + numVars +
          " variables is too large for the JVM at this time", e);
    }
  }

  @Override
  public void initialise(){
    super.initialise();
    MatrixUtils.fill(jointCount, 0);
    MatrixUtils.fill(smallMarginalCounts, 0);
    MatrixUtils.fill(bigMarginalCounts, 0);
  }
  
  /**
   * Given multiple time samples of a homogeneous array of variables (states),
   * add the observations of all sets of numVars of these
   * Do this for every time point
   *  
   * @param states 2D array of values of an array of variables
   *  at many observations (first index is time, second is variable index)
   */
  public void addObservations(int[][] states) throws Exception {
    int[] jointStates = MatrixUtils.computeCombinedValues(states, base);
    for (int t = 0; t < states.length; t++) {
      for (int i = 0; i < numVars; i++) {
        // Extract values of the 1D and the (N-1)D marginals
        int thisValue = states[t][i];
        int bigMarginalState = computeBigMarginalState(jointStates[t], i, thisValue);

        // Update counts
        bigMarginalCounts[i][bigMarginalState]++;
        smallMarginalCounts[i][thisValue]++;
      }
      jointCount[jointStates[t]]++;
      observations++;
    }
  }

  @Override
  public double computeAverageLocalOfObservations() {
    
    int[] jointTuple = new int[numVars];
    checkedFirst = false;
    try {
      average = computeForGivenTupleFromVarIndex(jointTuple, 0);
    } catch (Exception e) {
      System.out.println("Something went wrong during the calculation.");
      average = -1;
    }
        
    return average;
  }
  
  /**
   * Private utility to compute the contribution to the measure for all tuples
   * starting with tuple[0..(fromIndex-1)].
   * 
   * @param tuple
   * @param fromIndex
   * @return
   */
  public double computeForGivenTupleFromVarIndex(int[] tuple, int fromIndex) throws Exception {
    double miCont = 0;
    if (fromIndex == numVars) {
      // The whole tuple is filled in, so compute the contribution to the MI from this tuple
      int jointValue = MatrixUtils.computeCombinedValues(new int[][] {tuple}, base)[0];

      if (jointCount[jointValue] == 0) {
        // This joint state does not occur, so it makes no contribution here
        return 0;
      }

      double jointProb = (double) jointCount[jointValue] / (double) observations;
      double localValue = computeLocalValueForTuple(tuple, jointValue);
      miCont = jointProb * localValue;

    } else {
      // Fill out the next part of the tuple and make the recursive calls
      for (int v = 0; v < base; v++) {
        tuple[fromIndex] = v;
        miCont += computeForGivenTupleFromVarIndex(tuple, fromIndex + 1);
      }
    }
    return miCont;
  }

  /**
   * Shortcut method to initialise the calculator, add observations and compute
   * the average measure in one line.
   *
   * @param state series of multivariate observations
   *  (first index is time or observation index, second is variable number)
   */
  public double compute(int[][] states) throws Exception {
    initialise();
    addObservations(states);
    return computeAverageLocalOfObservations();
  }

  /**
   * Internal method to update maximum and minimum values of local information
   * measures.
   *
   * @param localValue instance of computed local information measure
   */
  protected void checkLocals(double localValue) {
    if (!checkedFirst) {
      max = localValue;
      min = localValue;
      checkedFirst = true;
    } else {
      if (localValue > max) {
        max = localValue;
      }
      if (localValue < min) {
        min = localValue;
      }
    }
  }

  /**
   * Method to be implemented by all child classes to compute the local value
   * of the measure for a given tuple.
   *
   * @param tuple state of the system at a given time (index is variable number)
   * @param jointValue <code>int</code> representing the state of the system
   */
  protected abstract double computeLocalValueForTuple(int[] tuple, int jointValue)
      throws Exception;

  /**
   * Method to be implemented by all child classes to compute the local value
   * of the measure for a given tuple.
   *
   * @param tuple state of the system at a given time (index is variable number)
   */
  protected double computeLocalValueForTuple(int[] tuple) throws Exception {
    int jointValue = MatrixUtils.computeCombinedValues(new int[][] {tuple}, base)[0];
    return computeLocalValueForTuple(tuple, jointValue);
  }

  /**
   * Small utility function to compute the state of the system excluding one variable.
   *
   * @param jointState state of the full system
   * @param varIdx index of the variable to be excluded
   * @param varValue value of the variable in question in the system state
   */
  protected int computeBigMarginalState(int jointState, int varIdx, int varValue) {
    int bigMarginalState = jointState - varValue*MathsUtils.power(base, numVars - varIdx - 1);
    return bigMarginalState;
  }

}

