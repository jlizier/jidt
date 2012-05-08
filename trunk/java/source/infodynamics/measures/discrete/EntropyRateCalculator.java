package infodynamics.measures.discrete;

/**
 * Compute average and local entropy rates
 * 
 * Usage:
 * 1. Continuous accumulation of observations:
 *   Call: a. initialise()
 *         b. addObservations() several times over
 *         c. computeLocalFromPreviousObservations()
 * 2. Standalone:
 *   Call: localActiveInformation()
 * 
 * @author Joseph Lizier
 *
 */
public class EntropyRateCalculator extends SingleAgentMeasureInContextOfPastCalculator {

	/**
	 * User to create new instances through this factory method.
	 * This allows us to return an efficient calculator for
	 * base 2, for example, without the user needing to have 
	 * knowledge of this. (This is obselete anyway)
	 * @param base
	 * @param history
	 * 
	 * @return
	 */
	public static EntropyRateCalculator newInstance(int base, int history) {
		return new EntropyRateCalculator(base, history);
	}
	
	public EntropyRateCalculator(int base, int history) {
		
		super(base, history);

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
		// increment the count of observations:
		observations += (rows - k)*columns; 
		
		// Initialise and store the current previous value for each column
		int[] prevVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			prevVal[c] = 0;
			for (int p = 0; p < k; p++) {
				prevVal[c] *= base;
				prevVal[c] += states[p][c];
			}
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int r = k; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				nextVal = states[r][c];
				nextPastCount[nextVal][prevVal[c]]++;
				pastCount[prevVal[c]]++;
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
		if (timeSteps == 0) {
			return;
		}
		int agentRows = states[0].length;
		if (agentRows == 0) {
			return;
		}
		int agentColumns = states[0][0].length;
		// increment the count of observations:
		observations += (timeSteps - k) * agentRows * agentColumns; 
		
		// Initialise and store the current previous value for each column
		int[][] prevVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				prevVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					prevVal[r][c] *= base;
					prevVal[r][c] += states[p][r][c];
				}
			}
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int t = k; t < timeSteps; t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					// Add to the count for this particular transition:
					// (cell's assigned as above)
					nextVal = states[t][r][c];
					nextPastCount[nextVal][prevVal[r][c]]++;
					pastCount[prevVal[r][c]]++;
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
		// increment the count of observations:
		observations += (rows - k); 
		
		// Initialise and store the current previous value for each column
		int prevVal = 0; 
		prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][col];
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int r = k; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			nextVal = states[r][col];
			nextPastCount[nextVal][prevVal]++;
			pastCount[prevVal]++;
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][col]];
			prevVal *= base;
			prevVal += states[r][col];
		}
	}

	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 */
	public void addObservations(int states[][][], int agentIndex1, int agentIndex2) {
		int timeSteps = states.length;
		// increment the count of observations:
		observations += (timeSteps - k); 
		
		// Initialise and store the current previous value for this column
		int prevVal = 0; 
		prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][agentIndex1][agentIndex2];
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int r = k; r < timeSteps; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			nextVal = states[r][agentIndex1][agentIndex2];
			nextPastCount[nextVal][prevVal]++;
			pastCount[prevVal]++;
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][agentIndex1][agentIndex2]];
			prevVal *= base;
			prevVal += states[r][agentIndex1][agentIndex2];
		}
	}

	/**
	 * Returns the average local active information storage from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return
	 */
	public double computeAverageLocalOfObservations() {
		double entRate = 0.0;
		double entRateCont = 0.0;

		max = 0;
		min = 0;
		double logTerm = 0;
		for (int nextVal = 0; nextVal < base; nextVal++) {
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) pastCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) nextPastCount[nextVal][prevVal] / (double) observations;
				// Compute entropy rate contribution:
				if (p_joint > 0.0) {
					logTerm = p_joint / p_prev;
					// Entropy rate takes the negative log:
					double localValue = - Math.log(logTerm) / log_2;
					entRateCont = p_joint * localValue;
					if (localValue > max) {
						max = localValue;
					} else if (localValue < min) {
						min = localValue;
					}
				} else {
					entRateCont = 0.0;
				}
				entRate += entRateCont;
			}
		}
		
		average = entRate;
		return entRate;
	}
	
	/**
	 * Computes local entropy rate for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method to be used for homogeneous agents only
	 *  
	 * @param states 1st index is time, 2nd index is agent number
	 * @return
	 */
	public double[][] computeLocalFromPreviousObservations(int states[][]){
		int rows = states.length;
		int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localEntRate = new double[rows][columns];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int[] prevVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			prevVal[c] = 0;
			for (int p = 0; p < k; p++) {
				prevVal[c] *= base;
				prevVal[c] += states[p][c];
			}
		}
		int nextVal;
		double logTerm = 0.0;
		for (int r = k; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				nextVal = states[r][c];
				logTerm = ( (double) nextPastCount[nextVal][prevVal[c]] ) /
				  		  ( (double) pastCount[prevVal[c]] );
				// Entropy rate takes the negative log:
				localEntRate[r][c] = - Math.log(logTerm) / log_2;
				average += localEntRate[r][c];
				if (localEntRate[r][c] > max) {
					max = localEntRate[r][c];
				} else if (localEntRate[r][c] < min) {
					min = localEntRate[r][c];
				}
				// Update the previous value:
				prevVal[c] -= maxShiftedValue[states[r-k][c]];
				prevVal[c] *= base;
				prevVal[c] += states[r][c];
			}
		}
		average = average/(double) (columns * (rows - k));
		
		return localEntRate;
		
	}
	
	/**
	 * Computes local entropy rate for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method to be used for homogeneous agents only
	 *  
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @return
	 */
	public double[][][] computeLocalFromPreviousObservations(int states[][][]){
		int timeSteps = states.length;
		int agentRows = states[0].length;
		int agentColumns = states[0][0].length;

		// Allocate for all time steps even though we'll leave the first ones as zeros
		double[][][] localEntRate = new double[timeSteps][agentRows][agentColumns];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int[][] prevVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				prevVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					prevVal[r][c] *= base;
					prevVal[r][c] += states[p][r][c];
				}
			}
		}
		
		int nextVal;
		double logTerm = 0.0;
		for (int t = k; t < timeSteps; t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					nextVal = states[t][r][c];
					logTerm = ( (double) nextPastCount[nextVal][prevVal[r][c]] ) /
					  		  ( (double) pastCount[prevVal[r][c]] );
					// Entropy rate takes the negative log:
					localEntRate[t][r][c] = - Math.log(logTerm) / log_2;
					average += localEntRate[t][r][c];
					if (localEntRate[t][r][c] > max) {
						max = localEntRate[t][r][c];
					} else if (localEntRate[t][r][c] < min) {
						min = localEntRate[t][r][c];
					}
					// Update the previous value:
					prevVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
					prevVal[r][c] *= base;
					prevVal[r][c] += states[t][r][c];
				}
			}
		}
		average = average/(double) (agentRows * agentColumns * (timeSteps - k));
		
		return localEntRate;
		
	}

	/**
	 * Computes local entropy rate for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states 1st index is time, 2nd index is agent number
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][], int col){
		int rows = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localEntRate = new double[rows];
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous value for each column
		int prevVal = 0; 
		prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][col];
		}
		int nextVal;
		double logTerm = 0.0;
		for (int r = k; r < rows; r++) {
			nextVal = states[r][col];
			logTerm = ( (double) nextPastCount[nextVal][prevVal] ) /
			  		  ( (double) pastCount[prevVal] );
			// Entropy rate takes the negative log:
			localEntRate[r] = - Math.log(logTerm) / log_2;
			average += localEntRate[r];
			if (localEntRate[r] > max) {
				max = localEntRate[r];
			} else if (localEntRate[r] < min) {
				min = localEntRate[r];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][col]];
			prevVal *= base;
			prevVal += states[r][col];
		}
		average = average/(double) (rows - k);
		
		return localEntRate;
		
	}

	/**
	 * Computes local entropy rate for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][][], int agentIndex1, int agentIndex2){
		int timeSteps = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localEntRate = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous value for each column
		int prevVal = 0; 
		prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][agentIndex1][agentIndex2];
		}
		int nextVal;
		double logTerm = 0.0;
		for (int t = k; t < timeSteps; t++) {
			nextVal = states[t][agentIndex1][agentIndex2];
			logTerm = ( (double) nextPastCount[nextVal][prevVal] ) /
			  		  ( (double) pastCount[prevVal] );
			// Entropy rate takes the negative log:
			localEntRate[t] = - Math.log(logTerm) / log_2;
			average += localEntRate[t];
			if (localEntRate[t] > max) {
				max = localEntRate[t];
			} else if (localEntRate[t] < min) {
				min = localEntRate[t];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k][agentIndex1][agentIndex2]];
			prevVal *= base;
			prevVal += states[t][agentIndex1][agentIndex2];
		}
		average = average/(double) (timeSteps - k);
		
		return localEntRate;
		
	}	
}
