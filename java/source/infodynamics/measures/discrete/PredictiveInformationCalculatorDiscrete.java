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

import infodynamics.utils.MatrixUtils;

/**
 * <p> Predictive information calculator for univariate discrete (int[]) data.
 * See definition of Predictive information (PI) by Bialek et al. below,
 * also is a form of the Excess Entropy (see Crutchfield and Feldman below),
 * Basically, PI is the mutual information between the past <i>state</i>
 * of a time-series process <i>X</i> and its future <i>state</i>.
 * The past <i>state</i> at time <code>n</code>
 * is represented by an embedding vector of <code>k</code> values from <code>X_n</code> backwards,
 * each separated by <code>\tau</code> steps, giving
 * <code><b>X^k_n</b> = [ X_{n-(k-1)\tau}, ... , X_{n-\tau}, X_n]</code>.
 * We call <code>k</code> the embedding dimension, and <code>\tau</code>
 * the embedding delay (only delay = 1 is implemented at the moment).
 * The future <i>state</i> at time <code>n</code>
 * is defined similarly into the future:
 * each separated by <code>\tau</code> steps, giving
 * <code><b>X^k+_n</b> = [ X_{n+1}, X_{n+\tau}, X_{n+(k-1)\tau}]</code>.
 * PI is then the mutual information between <b>X^k_n</b> and <b>X^k+_n</b>.</p>
 * 
 * <p>Usage of the class is intended to follow this paradigm:</p>
 * <ol>
 * 		<li>Construct the calculator:
 * 			{@link #PredictiveInformationCalculatorDiscrete(int, int)};</li>
 *		<li>Initialise the calculator using {@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			sets of {@link #addObservations(int[])} methods, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average entropy: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>local entropy values, such as {@link #computeLocal(int[])};</li>
 * 				<li>and variants of these.</li>
 * 			</ul>
 * 		</li>
 * 		<li>As an alternative to steps 3 and 4, the user may undertake
 * 			standalone computation from a single set of observations, via
 *  		e.g.: {@link #computeLocal(int[])},
 *  		{@link #computeAverageLocal(int[])} etc.</li>
 * 		<li>
 * 		Return to step 2 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * <p>TODO Tidy up the Javadocs for
 *  the methods, which are somewhat preliminary</p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Bialek, W., Nemenman, I., and Tishby, N.,
 *  <a href="http://dx.doi.org/10.1016/S0378-4371(01)00444-7">
 * 	"Complexity through nonextensivity"</a>,
 *  Physica A, 302, 89-99. (2001).</li>
 * 	<li>J. P. Crutchfield, D. P. Feldman,
 *  <a href="http://dx.doi.org/10.1063/1.1530990">
 * 	"Regularities Unseen, Randomness Observed: Levels of Entropy Convergence"</a>,
 *  Chaos, Vol. 13, No. 1. (2003), pp. 25-54.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class PredictiveInformationCalculatorDiscrete extends SingleAgentMeasureDiscreteInContextOfPastCalculator {

	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param numDiscreteValues Number of discrete values (e.g. 2 for binary states)
	 * @param blockLength
	 * @deprecated
	 * @return
	 */
	public static PredictiveInformationCalculatorDiscrete newInstance(int numDiscreteValues, int blockLength) {
		return new PredictiveInformationCalculatorDiscrete(numDiscreteValues, blockLength);
	}
	
	/**
	 * Construct a new instance with default base 2 and history 1
	 */
	public PredictiveInformationCalculatorDiscrete() {
		this(2, 1);
	}
	
	/**
	 * Construct a new instance
	 * 
	 * @param numDiscreteValues number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param blockLength embedded history length of the past and future to use -
	 *        this is k in Schreiber's notation.
	 */
	public PredictiveInformationCalculatorDiscrete(int numDiscreteValues, int blockLength) {
		super(numDiscreteValues, blockLength, true);		
	}

	@Override
	public void initialise(int base, int blockLength) {
		boolean baseOrHistoryChanged = (this.base != base) || (k != blockLength);
		super.initialise(base, blockLength);

		if (baseOrHistoryChanged || (nextPastCount == null)) {
			// Create new storage for counts of observations.
			// Need to manage this here since next is now a block variable rather than univariate in the super.
			try {
				nextPastCount = new int[base_power_k][base_power_k];
				pastCount = new int[base_power_k];
				nextCount = new int[base_power_k];
			} catch (OutOfMemoryError e) {
				// Allow any Exceptions to be thrown, but catch and wrap
				//  Error as a RuntimeException
				throw new RuntimeException("Requested memory for the base " +
						base + " with k=" + blockLength +
						") is too large for the JVM at this time", e);
			}
		} else {
			MatrixUtils.fill(nextPastCount, 0);
			MatrixUtils.fill(pastCount, 0);
			MatrixUtils.fill(nextCount, 0);
		}
	}
	
	@Override
	public void initialise(){
		initialise(base, k);
	}
		
	/**
 	 * Add observations in to our estimates of the pdfs.
	 *
	 * @param timeSeries time series of agent states
	 */
	public void addObservations(int timeSeries[]) {
		int timeSteps = timeSeries.length;
		// increment the count of observations:
		// (we miss out on k observations at the start and k-1 at the end)
		if (timeSteps - k - (k-1) <= 0) {
			// Nothing to do
			return;
		}
		observations += (timeSteps - k - (k-1)); 
		
		// Initialise and store the current previous value 
		//  and next values (next val set for t = k-1)
		int prevVal = 0;
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += timeSeries[p];
			nextVal *= base;
			nextVal += timeSeries[k-1+p];
		}
		
		// 1. Count the tuples observed
		for (int t = k; t < timeSteps - (k-1); t++) {
			// Update the next value:
			nextVal -= maxShiftedValue[timeSeries[t-1]];
			nextVal *= base;
			nextVal += timeSeries[k-1+t];
			// Update the counts
			nextPastCount[nextVal][prevVal]++;
			pastCount[prevVal]++;
			nextCount[nextVal]++;
			// Update the previous value:
			prevVal -= maxShiftedValue[timeSeries[t-k]];
			prevVal *= base;
			prevVal += timeSeries[t];
		}		
	}

	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states 1st index is time, 2nd index is agent number
	 */
	public void addObservations(int states[][]) {
		int rows = states.length;
		int columns = states[0].length;
		// (we miss out on k observations at the start and k-1 at the end)
		if (rows - k - (k-1) <= 0) {
			// Nothing to do
			return;
		}
		observations += (rows - k - (k-1))*columns; 
		
		// Initialise and store the current previous and
		//  next val (next val set for t = k-1) value for each column
		int[] prevVal = new int[columns]; 
		int[] nextVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			prevVal[c] = 0;
			nextVal[c] = 0;
			for (int p = 0; p < k; p++) {
				prevVal[c] *= base;
				prevVal[c] += states[p][c];
				nextVal[c] *= base;
				nextVal[c] += states[k-1+p][c];
			}
		}
		
		// 1. Count the tuples observed
		for (int r = k; r < rows - (k-1); r++) {
			for (int c = 0; c < columns; c++) {
				// Update the next value:
				nextVal[c] -= maxShiftedValue[states[r-1][c]];
				nextVal[c] *= base;
				nextVal[c] += states[k-1+r][c];
				// Update the counts
				nextPastCount[nextVal[c]][prevVal[c]]++;
				pastCount[prevVal[c]]++;
				nextCount[nextVal[c]]++;
				// Update the previous value:
				prevVal[c] -= maxShiftedValue[states[r-k][c]];
				prevVal[c] *= base;
				prevVal[c] += states[r][c];
			}
		}		
	}

	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 */
	public void addObservations(int states[][][]) {
		int timeSteps = states.length;
		// (we miss out on k observations at the start and k-1 at the end)
		if (timeSteps - k - (k-1) <= 0) {
			// Nothing to do
			return;
		}
		int agentRows = states[0].length;
		if (agentRows == 0) {
			return;
		}
		int agentColumns = states[0][0].length;
		// increment the count of observations:
		observations += (timeSteps - k - (k-1)) * agentRows * agentColumns; 

		// Initialise and store the current previous and
		// next (next val set for t = k-1) value for each column
		int[][] prevVal = new int[agentRows][agentColumns];
		int[][] nextVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				prevVal[r][c] = 0;
				nextVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					prevVal[r][c] *= base;
					prevVal[r][c] += states[p][r][c];
					nextVal[r][c] *= base;
					nextVal[r][c] += states[k-1+p][r][c];
				}
			}
		}
		
		// 1. Count the tuples observed
		for (int t = k; t < timeSteps - (k-1); t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					// Update the next value:
					nextVal[r][c] -= maxShiftedValue[states[t-1][r][c]];
					nextVal[r][c] *= base;
					nextVal[r][c] += states[k-1+t][r][c];
					// Update the counts
					nextPastCount[nextVal[r][c]][prevVal[r][c]]++;
					pastCount[prevVal[r][c]]++;
					nextCount[nextVal[r][c]]++;
					// Update the previous value:
					prevVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
					prevVal[r][c] *= base;
					prevVal[r][c] += states[t][r][c];
				}
			}
		}		
	}

	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd index is agent number
	 */
	public void addObservations(int states[][], int col) {
		int rows = states.length;
		// (we miss out on k observations at the start and k-1 at the end)
		if (rows - k - (k-1) <= 0) {
			// Nothing to do
			return;
		}
		observations += (rows - k - (k-1)); 
		
		// Initialise and store the current previous and
		//  next value (next val set for t = k-1) for each column
		int prevVal = 0; 
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][col];
			nextVal *= base;
			nextVal += states[k-1+p][col];
		}
		
		// 1. Count the tuples observed
		for (int r = k; r < rows - (k-1); r++) {
			// Update the next value:
			nextVal -= maxShiftedValue[states[r-1][col]];
			nextVal *= base;
			nextVal += states[k-1+r][col];
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			nextPastCount[nextVal][prevVal]++;
			pastCount[prevVal]++;
			nextCount[nextVal]++;					
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][col]];
			prevVal *= base;
			prevVal += states[r][col];
		}
	}

	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 */
	public void addObservations(int states[][][], int agentIndex1, int agentIndex2) {
		int timeSteps = states.length;
		// (we miss out on k observations at the start and k-1 at the end)
		if (timeSteps - k - (k-1) <= 0) {
			// Nothing to do
			return;
		}
		// increment the count of observations:
		observations += (timeSteps - k - (k-1)); 
		
		// Initialise and store the current previous and
		//  next value (next val set for t = k-1) for each column
		int prevVal = 0; 
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][agentIndex1][agentIndex2];
			nextVal *= base;
			nextVal += states[k-1+p][agentIndex1][agentIndex2];
		}
		
		// 1. Count the tuples observed
		for (int t = k; t < timeSteps - (k-1); t++) {
			// Update the next value:
			nextVal -= maxShiftedValue[states[t-1][agentIndex1][agentIndex2]];
			nextVal *= base;
			nextVal += states[k-1+t][agentIndex1][agentIndex2];
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			nextPastCount[nextVal][prevVal]++;
			pastCount[prevVal]++;
			nextCount[nextVal]++;
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k][agentIndex1][agentIndex2]];
			prevVal *= base;
			prevVal += states[t][agentIndex1][agentIndex2];
		}
	}

	@Override
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		double miCont = 0.0;

		max = 0;
		min = 0;
		for (int nextVal = 0; nextVal < base_power_k; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) pastCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) nextPastCount[nextVal][prevVal] / (double) observations;
				// Compute MI contribution:
				if (p_joint > 0.0) {
					double logTerm = p_joint / (p_next * p_prev);
					double localValue = Math.log(logTerm) / log_base;
					miCont = p_joint * localValue;
					if (localValue > max) {
						max = localValue;
					} else if (localValue < min) {
						min = localValue;
					}
				} else {
					miCont = 0.0;
				}
				mi += miCont;
			}
		}
		
		average = mi;
		return mi;
	}
	
	/**
	 * Computes local predictive info for the given values
	 *  
	 * @param next joint state of future (see code for {@link #addObservations(int[])}
	 *   for how the joint integer value is computed)
	 * @param past joint state of past (see code for {@link #addObservations(int[])}
	 *   for how the joint integer value is computed)
	 * @return
	 */
	public double computeLocalFromPreviousObservations(int next, int past){
		double logTerm = ( (double) nextPastCount[next][past] ) /
		  ( (double) nextCount[next] *
			(double) pastCount[past] );
		logTerm *= (double) observations;
		return Math.log(logTerm) / log_base;
	}

	/**
	 * Computes local predictive information for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 *  
	 * @param timeSeries time series of values
	 * @return array of local predictive information values
	 */
	public double[] computeLocalFromPreviousObservations(int timeSeries[]){
		int timeSteps = timeSeries.length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localPredictive = new double[timeSteps];
		if (timeSteps < k + (k-1)) {
			// Nothing to do
			return localPredictive;
		}

		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous and
		//  next value (next val set for t = k-1)
		int prevVal = 0;
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += timeSeries[p];
			nextVal *= base;
			nextVal += timeSeries[k-1+p];
		}

		double logTerm = 0.0;
		for (int t = k; t < timeSteps - (k-1); t++) {
			// Update the next value:
			nextVal -= maxShiftedValue[timeSeries[t-1]];
			nextVal *= base;
			nextVal += timeSeries[k-1+t];
			logTerm = ( (double) nextPastCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) pastCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localPredictive[t] = Math.log(logTerm) / log_base;
			average += localPredictive[t];
			if (localPredictive[t] > max) {
				max = localPredictive[t];
			} else if (localPredictive[t] < min) {
				min = localPredictive[t];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[timeSeries[t-k]];
			prevVal *= base;
			prevVal += timeSeries[t];
		}
		average = average/(double) (timeSteps - k - (k-1));
		
		return localPredictive;
		
	}

	/**
	 * Computes local predictive information for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only, since
	 *  all observations are pooled together.
	 *  
	 * @param timeSeries 1st index is time, 2nd index is agent number
	 * @return
	 */
	public double[][] computeLocalFromPreviousObservations(int timeSeries[][]){
		int rows = timeSeries.length;
		int columns = timeSeries[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localPredictive = new double[rows][columns];
		if (rows < k + (k-1)) {
			// Nothing to do
			return localPredictive;
		}

		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous and
		//  next value (next val set for t = k-1) for each column
		int[] prevVal = new int[columns];
		int[] nextVal = new int[columns];
		for (int c = 0; c < columns; c++) {
			prevVal[c] = 0;
			nextVal[c] = 0;
			for (int p = 0; p < k; p++) {
				prevVal[c] *= base;
				prevVal[c] += timeSeries[p][c];
				nextVal[c] *= base;
				nextVal[c] += timeSeries[k-1+p][c];
			}
		}
		double logTerm = 0.0;
		for (int r = k; r < rows - (k-1); r++) {
			for (int c = 0; c < columns; c++) {
				// Update the next value:
				nextVal[c] -= maxShiftedValue[timeSeries[r-1][c]];
				nextVal[c] *= base;
				nextVal[c] += timeSeries[k-1+r][c];
				logTerm = ( (double) nextPastCount[nextVal[c]][prevVal[c]] ) /
				  		  ( (double) nextCount[nextVal[c]] *
				  			(double) pastCount[prevVal[c]] );
				// Now account for the fact that we've
				//  just used counts rather than probabilities,
				//  and we've got two counts on the bottom
				//  but one count on the top:
				logTerm *= (double) observations;
				localPredictive[r][c] = Math.log(logTerm) / log_base;
				average += localPredictive[r][c];
				if (localPredictive[r][c] > max) {
					max = localPredictive[r][c];
				} else if (localPredictive[r][c] < min) {
					min = localPredictive[r][c];
				}
				// Update the previous value:
				prevVal[c] -= maxShiftedValue[timeSeries[r-k][c]];
				prevVal[c] *= base;
				prevVal[c] += timeSeries[r][c];
			}
		}
		average = average/(double) (columns * (rows - k - (k-1)));
		
		return localPredictive;
		
	}
	
	/**
	 * Computes local active information storage for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only, since
	 *  all observations are pooled together.
	 *  
	 * @param timeSeries 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @return
	 */
	public double[][][] computeLocalFromPreviousObservations(int timeSeries[][][]){
		int timeSteps = timeSeries.length;
		int agentRows = timeSeries[0].length;
		int agentColumns = timeSeries[0][0].length;

		// Allocate for all rows even though we'll leave the first and
		// last ones as zeros
		double[][][] localPredictive = new double[timeSteps][agentRows][agentColumns];
		if (timeSteps < k + (k-1)) {
			// Nothing to do
			return localPredictive;
		}
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous and
		//  next value (next val set for t = k-1) for each agent
		int[][] prevVal = new int[agentRows][agentColumns]; 
		int[][] nextVal = new int[agentRows][agentColumns]; 
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				prevVal[r][c] = 0;
				nextVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					prevVal[r][c] *= base;
					prevVal[r][c] += timeSeries[p][r][c];
					nextVal[r][c] *= base;
					nextVal[r][c] += timeSeries[k-1+p][r][c];
				}
			}
		}
		double logTerm = 0.0;
		for (int t = k; t < timeSteps - (k-1); t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					// Update the next value:
					nextVal[r][c] -= maxShiftedValue[timeSeries[t-1][r][c]];
					nextVal[r][c] *= base;
					nextVal[r][c] += timeSeries[k-1+t][r][c];
					logTerm = ( (double) nextPastCount[nextVal[r][c]][prevVal[r][c]] ) /
					  		  ( (double) nextCount[nextVal[r][c]] *
					  			(double) pastCount[prevVal[r][c]] );
					// Now account for the fact that we've
					//  just used counts rather than probabilities,
					//  and we've got two counts on the bottom
					//  but one count on the top:
					logTerm *= (double) observations;
					localPredictive[t][r][c] = Math.log(logTerm) / log_base;
					average += localPredictive[t][r][c];
					if (localPredictive[t][r][c] > max) {
						max = localPredictive[t][r][c];
					} else if (localPredictive[t][r][c] < min) {
						min = localPredictive[t][r][c];
					}
					// Update the previous value:
					prevVal[r][c] -= maxShiftedValue[timeSeries[t-k][r][c]];
					prevVal[r][c] *= base;
					prevVal[r][c] += timeSeries[t][r][c];
				}
			}
		}
		average = average/(double) (agentRows * agentColumns * (timeSteps - k - (k-1)));
		
		return localPredictive;
	}

	/**
	 * Computes local predictive information for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method is suitable for heterogeneous agents since the 
	 *  user specifies which agent to take the observations from.
	 *  
	 * @param states 1st index is time, 2nd index is agent number
	 * @param col gives the column index identifying the agent
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][], int col){
		int rows = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localPredictive = new double[rows];
		if (rows < k + (k-1)) {
			// Nothing to do
			return localPredictive;
		}
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous and 
		//  next value (next val set for t = k-1) for each column
		int prevVal = 0; 
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][col];
			nextVal *= base;
			nextVal += states[k-1+p][col];
		}
		double logTerm = 0.0;
		for (int r = k; r < rows - (k-1); r++) {
			// Update the next value:
			nextVal -= maxShiftedValue[states[r-1][col]];
			nextVal *= base;
			nextVal += states[k-1+r][col];
			logTerm = ( (double) nextPastCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) pastCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localPredictive[r] = Math.log(logTerm) / log_base;
			average += localPredictive[r];
			if (localPredictive[r] > max) {
				max = localPredictive[r];
			} else if (localPredictive[r] < min) {
				min = localPredictive[r];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][col]];
			prevVal *= base;
			prevVal += states[r][col];
		}
		average = average/(double) (rows - k - (k-1));
		
		return localPredictive;
		
	}
	
	/**
	 * Computes local predictive information for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method is suitable for heterogeneous agents, since the
	 *  relevant agent is identified
	 *  
	 * @param timeSeries 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param agentIndex1 gives the first index identifying the agent
	 * @param agentIndex2 gives the second index identifying the agent
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int timeSeries[][][], int agentIndex1, int agentIndex2){
		int timeSteps = timeSeries.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localPredictive = new double[timeSteps];
		if (timeSteps < k + (k-1)) {
			// Nothing to do
			return localPredictive;
		}
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous and 
		//  next value (next val set for t = k-1) for each column
		int prevVal = 0;
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += timeSeries[p][agentIndex1][agentIndex2];
			nextVal *= base;
			nextVal += timeSeries[k-1+p][agentIndex1][agentIndex2];
		}
		double logTerm = 0.0;
		for (int t = k; t < timeSteps - (k-1); t++) {
			// Update the next value:
			nextVal -= maxShiftedValue[timeSeries[t-1][agentIndex1][agentIndex2]];
			nextVal *= base;
			nextVal += timeSeries[k-1+t][agentIndex1][agentIndex2];
			logTerm = ( (double) nextPastCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) pastCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localPredictive[t] = Math.log(logTerm) / log_base;
			average += localPredictive[t];
			if (localPredictive[t] > max) {
				max = localPredictive[t];
			} else if (localPredictive[t] < min) {
				min = localPredictive[t];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[timeSeries[t-k][agentIndex1][agentIndex2]];
			prevVal *= base;
			prevVal += timeSeries[t][agentIndex1][agentIndex2];
		}
		average = average/(double) (timeSteps - k - (k-1));
		
		return localPredictive;
		
	}

	/**
	 * Writes the current probability distribution functions 
	 *  
	 * @return
	 */
	public void writePdfs() {
		double mi = 0.0;
		double miCont = 0.0;

		System.out.println("nextVal p(next) prevVal p(prev) p(joint) logTerm localVal");
		for (int nextVal = 0; nextVal < base_power_k; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) pastCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) nextPastCount[nextVal][prevVal] / (double) observations;
				// Compute MI contribution:
				if (p_joint > 0.0) {
					double logTerm = p_joint / (p_next * p_prev);
					double localValue = Math.log(logTerm) / log_base;
					miCont = p_joint * localValue;
					System.out.println(String.format("%7d    %.2f %7d    %.2f     %.2f    %.2f     %.2f",
							nextVal, p_next, prevVal, p_prev, p_joint, logTerm, localValue));
				} else {
					miCont = 0.0;
					System.out.println(String.format("%7d    %.2f %7d    %.2f     %.2f    %.2f     %.2f",
							nextVal, p_next, prevVal, p_prev, p_joint, 0.0, 0.0));
				}
				mi += miCont;
			}
		}
		System.out.println("Average is " + mi);
		
		return;
	}

	/**
	 * Utility function to compute the combined past values of x up to and including time step t
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 * 
	 * @param x time series
	 * @param t time step to compute joint past value up to
	 * @return joint discrete value computed from the (x_{t-k+1}, ... ,x_{t-1},x_{t})
	 */
	public int computePastValue(int[] x, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += x[t - k + 1 + p];
		}
		return pastVal;
	}

	/**
	 * Utility function to compute the combined past values of x up to
	 * and including time step t for a given variable number within x
	 *  (i.e. (x_{t-k+1,i}, ... ,x_{t-1,i},x_{t,i}))
	 * 
	 * @param x 2D multivariate time series
	 * @param i column or agent number within x.
	 * @param t time step to compute joint past value up to
	 * @return joint discrete value computed from the (x_{t-k+1,i}, ... ,x_{t-1,i},x_{t,i}))
	 */
	public int computePastValue(int[][] x, int i, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += x[t - k + 1 + p][i];
		}
		return pastVal;
	}

	/**
	 * Utility function to compute the combined past values of x up to
	 * and including time step t for a given variable number within x
	 *  (i.e. (x_{t-k+1,i,j}, ... ,x_{t-1,i,j},x_{t,i,j}))
	 * 
	 * @param x 3D multivariate time series
	 * @param i 1st index to identify agent within x.
	 * @param i 2nd index to identify agent within x.
	 * @param t time step to compute joint past value up to
	 * @return joint discrete value computed from the (x_{t-k+1,i}, ... ,x_{t-1,i},x_{t,i}))
	 */
	public int computePastValue(int[][][] x, int i, int j, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += x[t - k + 1 + p][i][j];
		}
		return pastVal;
	}
}
