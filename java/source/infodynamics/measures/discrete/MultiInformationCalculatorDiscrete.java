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
 * <p>Multi-information or integration calculator for multivariate discrete data.
 * That is, it is applied to <code>double[][]</code> data, where the first index
 * is observation number or time, and the second is variable number.
 * See Tononi et al. below for the definition of multi-information/integration.</p>
 * </p>
 * 
 * <p>Usage of the class is intended to follow this paradigm:</p>
 * <ol>
 * 		<li>Construct the calculator:
 * 			{@link #MultiInformationCalculatorDiscrete(int, int)};</li>
 *		<li>Initialise the calculator using {@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			sets of {@link #addObservations(int[][], int[])} methods, then</li>
 * 		<li>Compute the required quantities, being (at this stage):
 * 			<ul>
 * 				<li>the average MI: {@link #computeAverageLocalOfObservations()};</li>
 * 			</ul>
 * 		</li>
 * 		<li>
 * 		Return to step 2 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>"G. Tononi, O. Sporns, G. M. Edelman,
 *  <a href="http://dx.doi.org/10.1073/pnas.91.11.5033">"A measure for
 *  brain complexity:
 * 	relating functional segregation and integration in the nervous system"</a>
 *  Proceedings of the National Academy of Sciences, Vol. 91, No. 11.
 *  (1994), pp. 5033-5037.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MultiInformationCalculatorDiscrete extends InfoMeasureCalculatorDiscrete {
	
	private int[]	jointCount = null; // jointCount[jointState]
	private int[][] marginalCounts = null; // marginalCounts[marginalIndex][state]
	
	private int numVars;
	private int jointStates;
	private boolean checkedFirst = false;

	/**
	 * Construct new instance with default base 2 and numVars 2
	 */
	public MultiInformationCalculatorDiscrete() {
		this(2, 2);
	}
	
	/**
	 * Construct an instance
	 * 
	 * @alphabetSize number of symbols for each variable.
	 *        E.g. binary variables are in base-2, so alphabetSize would be 2.
	 * @numVars numbers of joint variables that multi-info
	 *     will be computed over.
	 */
	public MultiInformationCalculatorDiscrete(int alphabetSize, int numVars) {
		super(alphabetSize);
		changeParams(alphabetSize, numVars);
	}

	/**
	 * Update the parameters using numVars
	 * 
	 * @param alphabetSize
	 * @param numVars
	 * @return
	 */
	private boolean changeParams(int alphabetSize, int numVars) {
		boolean paramsChanged = (this.alphabetSize != alphabetSize) || (this.numVars != numVars);
		
		this.numVars = numVars;
		jointStates = MathsUtils.power(alphabetSize, numVars);
		return paramsChanged;
	}
	
	/**
	 * Initialise (possible) updating the alphabetSize and number of variables
	 * 
	 * @param alphabetSize
	 * @param numVars
	 */
	public void initialise(int alphabetSize, int numVars) {
		boolean paramsChanged = changeParams(alphabetSize, numVars);
		super.initialise(alphabetSize);

		if(paramsChanged || (jointCount == null)) {
			// Create storage for counts of observations
			try {
				jointCount = new int[jointStates];
				marginalCounts = new int[numVars][alphabetSize];
			} catch (OutOfMemoryError e) {
				// Allow any Exceptions to be thrown, but catch and wrap
				//  Error as a RuntimeException
				throw new RuntimeException("Requested memory for the alphabetSize " +
						alphabetSize + " with " + numVars +
						" variables is too large for the JVM at this time", e);
			}
		} else {
			MatrixUtils.fill(jointCount, 0);
			MatrixUtils.fill(marginalCounts, 0);
		}
	}
	
	@Override
	public void initialise(){
		initialise(alphabetSize, numVars);
	}
	
	/**
	 * Given multiple time samples of a homogeneous array of variables (states),
	 * add the observations of all sets of numVars of these
	 * Do this for every time point
	 *  
	 * @param states 2D array of values of an array of variables
	 *  at many observations (first index is time, second is variable index)
	 */
	public void addObservations(int[][] states) {
		for (int t = 0; t < states.length; t++) {
			// Add the marginal observations in, and compute the joint state value
			//  for this sample at time t
			int jointValue = 0;
			for (int i = 0; i < numVars; i++) {
				int thisValue = states[t][i];
				marginalCounts[i][thisValue]++;
				jointValue *= alphabetSize;
				jointValue += thisValue;
			}
			jointCount[jointValue]++;
			observations++;
		}
	}
	
	/**
	 * Given one time sample of a homogeneous array of variables (states),
	 * add the observations of all sets of numVars of these, defined
	 * by the offsets in groupOffsets from every point in the array.
	 *  
	 * @param states current values of an array of variables 
	 * @param groupOffsets offsets for the numVars set from a given
	 *   data point
	 */
	public void addObservations(int[] states, int[] groupOffsets) {
		for (int c = 0; c < states.length; c++) {
			// Add the marginal observations in, and compute the joint state value
			int jointValue = 0;
			for (int i = 0; i < numVars; i++) {
				int thisValue = states[(c + groupOffsets[i] + states.length) % states.length];
				marginalCounts[i][thisValue]++;
				jointValue *= alphabetSize;
				jointValue += thisValue;
			}
			jointCount[jointValue]++;
			observations++;
		}
	}
	
	/**
	 * Given one time sample of a homogeneous array of variables (states),
	 * add the observations of one set of numVars of these, defined
	 * by the offsets in groupOffsets from destinationIndex.
	 *  
	 * @param states current values of an array of variables 
	 * @param destinationIndex which variable to center
	 *  our offsets from and take the sample,
	 * @param groupOffsets offsets for the numVars set from a given
	 *   data point
	 */
	public void addObservations(int[] states, int destinationIndex, int[] groupOffsets) {
		// Add the marginal observations in, and compute the joint state value
		int jointValue = 0;
		for (int i = 0; i < numVars; i++) {
			int thisValue = states[(destinationIndex + groupOffsets[i] + states.length) % states.length];
			marginalCounts[i][thisValue]++;
			jointValue *= alphabetSize;
			jointValue += thisValue;
		}
		jointCount[jointValue]++;
		observations++;
	}

	/**
	 * Given multiple time samples of a homogeneous array of variables (states),
	 * add the observations of all sets of numVars of these, defined
	 * by the offsets in groupOffsets from every point in the array.
	 * Do this for every time point
	 *  
	 * @param states 2D array of values of an array of variables
	 *  at many observations (first index is time, second is variable index) 
	 * @param groupOffsets offsets for the numVars set from a given
	 *   data point
	 */
	public void addObservations(int[][] states, int[] groupOffsets) {
		for (int t = 0; t < states.length; t++) {
			for (int c = 0; c < states[t].length; c++) {
				// Add the marginal observations in, and compute the joint state value
				int jointValue = 0;
				for (int i = 0; i < numVars; i++) {
					int thisValue = states[t][(c + groupOffsets[i] + states.length) % states.length];
					marginalCounts[i][thisValue]++;
					jointValue *= alphabetSize;
					jointValue += thisValue;
				}
				jointCount[jointValue]++;
				observations++;
			}
		}
	}
	
	@Override
	public double computeAverageLocalOfObservations() {
		
		int[] jointTuple = new int[numVars];
		checkedFirst = false;
		average = computeMiOfGivenTupleFromVarIndex(jointTuple, 0);
				
		return average;
	}
	
	/**
	 * Private utility to compute the contribution to the MI for all tuples starting with tuple[0..(fromIndex-1)].
	 * 
	 * @param tuple
	 * @param fromIndex
	 * @return
	 */
	public double computeMiOfGivenTupleFromVarIndex(int[] tuple, int fromIndex) {
		double miCont = 0;
		if (fromIndex == numVars) {
			// The whole tuple is filled in, so compute the contribution to the MI from this tuple
			double prodMarginalProbs = 1.0;
			int jointValue = 0;
			for (int i = 0; i < numVars; i++) {
				prodMarginalProbs *= (double) marginalCounts[i][tuple[i]] / (double) observations;
				jointValue *= alphabetSize;
				jointValue += tuple[i];
			}
			if (jointCount[jointValue] == 0) {
				// This joint state does not occur, so it makes no contribution here
				return 0;
			}
			double jointProb = (double) jointCount[jointValue] / (double) observations;
			double logValue = jointProb / prodMarginalProbs;
			double localValue = Math.log(logValue) / log_2;
			if (jointProb > 0.0) {
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
			miCont = jointProb * localValue;
		} else {
			// Fill out the next part of the tuple and make the recursive calls
			for (int v = 0; v < alphabetSize; v++) {
				tuple[fromIndex] = v;
				miCont += computeMiOfGivenTupleFromVarIndex(tuple, fromIndex + 1);
			}
		}
		return miCont;
	}
}
