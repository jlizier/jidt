package infodynamics.measures.discrete;

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

/**
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
public class ActiveInformationCalculator {

	private double average = 0.0;
	private double max = 0.0;
	private double min = 0.0;
	private int observations = 0;
	private int k = 0; // history length k.
	private int base = 0; // number of individual states.
	private int[][]	jointCount = null; // Count for (i[t+1], i[t]) tuples
	private int[] prevCount = null; // Count for i[t]		
	private int[] nextCount = null; // Count for i[t+1]
	private int[] maxShiftedValue = null; // states * (base^(history-1))
	
	private int base_power_k = 0;
	private double log_base = 0;
	
	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base
	 * @param history
	 * 
	 * @return
	 */
	public static ActiveInformationCalculator newInstance(int base, int history) {
		return new ActiveInformationCalculator(base, history);
	}
	
	public ActiveInformationCalculator(int base, int history) {
		super();
		
		this.base = base;
		k = history;
		base_power_k = MathsUtils.power(base, k);
		log_base = Math.log(base);
		
		if (history < 1) {
			throw new RuntimeException("History k " + history + " is not >= 1 for Active info storage Calculator");
		}
		if (k > Math.log(Integer.MAX_VALUE) / log_base) {
			throw new RuntimeException("Base and history combination too large");
		}

		// Create storage for counts of observations
		jointCount = new int[base][base_power_k];
		prevCount = new int[base_power_k];
		nextCount = new int[base];

		// Create constants for tracking prevValues
		maxShiftedValue = new int[base];
		for (int v = 0; v < base; v++) {
			maxShiftedValue[v] = v * MathsUtils.power(base, k-1);
		}
	}

	/**
	 * Initialise calculator, preparing to take observation sets in
	 * Should be called prior to any of the addObservations() methods.
	 * You can reinitialise without needing to create a new object.
	 *
	 */
	public void initialise(){
		average = 0.0;
		max = 0.0;
		min = 0.0;
		observations = 0;
		
		MatrixUtils.fill(jointCount, 0);
		MatrixUtils.fill(prevCount, 0);
		MatrixUtils.fill(nextCount, 0);
	}
		
	/**
 	 * Add observations in to our estimates of the pdfs.
	 *
	 * @param states time series of agent states
	 */
	public void addObservations(int states[]) {
		int timeSteps = states.length;
		// increment the count of observations:
		observations += (timeSteps - k); 
		
		// Initialise and store the current previous value for each column
		int prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p];
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int t = k; t < timeSteps; t++) {
			// Add to the count for this particular transition:
			nextVal = states[t];
			jointCount[nextVal][prevVal]++;
			prevCount[prevVal]++;
			nextCount[nextVal]++;					
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k]];
			prevVal *= base;
			prevVal += states[t];
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
				jointCount[nextVal][prevVal[c]]++;
				prevCount[prevVal[c]]++;
				nextCount[nextVal]++;					
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
					jointCount[nextVal][prevVal[r][c]]++;
					prevCount[prevVal[r][c]]++;
					nextCount[nextVal]++;					
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
			jointCount[nextVal][prevVal]++;
			prevCount[prevVal]++;
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
		// increment the count of observations:
		observations += (timeSteps - k); 
		
		// Initialise and store the current previous value for each column
		int prevVal = 0; 
		prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][agentIndex1][agentIndex2];
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int t = k; t < timeSteps; t++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			nextVal = states[t][agentIndex1][agentIndex2];
			jointCount[nextVal][prevVal]++;
			prevCount[prevVal]++;
			nextCount[nextVal]++;					
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k][agentIndex1][agentIndex2]];
			prevVal *= base;
			prevVal += states[t][agentIndex1][agentIndex2];
		}
	}

	/**
	 * Returns the average local active information storage from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return
	 */
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		double miCont = 0.0;

		max = 0;
		min = 0;
		for (int nextVal = 0; nextVal < base; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) prevCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) jointCount[nextVal][prevVal] / (double) observations;
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
	 * Returns the average local entropy rate from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return
	 */
	public double computeAverageLocalEntropyRateOfObservations() {
		double entRate = 0.0;
		double entRateCont = 0.0;

		for (int nextVal = 0; nextVal < base; nextVal++) {
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) prevCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) jointCount[nextVal][prevVal] / (double) observations;
				// Compute entropy rate contribution:
				if (p_joint > 0.0) {
					double logTerm = p_joint / p_prev;
					// Entropy rate takes the negative log:
					double localValue = - Math.log(logTerm) / log_base;
					entRateCont = p_joint * localValue;
				} else {
					entRateCont = 0.0;
				}
				entRate += entRateCont;
			}
		}
		
		return entRate;
	}
	
	/**
	 * Computes local active info storage for the given (single)
	 *  specific values
	 *  
	 * @param destNext
	 * @param destPast
	 * @param sourceCurrent
	 * @return
	 */
	public double computeLocalFromPreviousObservations(int next, int past){
		double logTerm = ( (double) jointCount[next][past] ) /
		  ( (double) nextCount[next] *
			(double) prevCount[past] );
		logTerm *= (double) observations;
		return Math.log(logTerm) / log_base;
	}

	/**
	 * Computes local active information storage for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 *  
	 * @param states time series of states
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[]){
		int timeSteps = states.length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localActive = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p];
		}

		int nextVal;
		double logTerm = 0.0;
		for (int t = k; t < timeSteps; t++) {
			nextVal = states[t];
			logTerm = ( (double) jointCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) prevCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localActive[t] = Math.log(logTerm) / log_base;
			average += localActive[t];
			if (localActive[t] > max) {
				max = localActive[t];
			} else if (localActive[t] < min) {
				min = localActive[t];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k]];
			prevVal *= base;
			prevVal += states[t];
		}
		average = average/(double) (timeSteps - k);
		
		return localActive;
		
	}

	/**
	 * Computes local active information storage for the given
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
		double[][] localActive = new double[rows][columns];
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
				logTerm = ( (double) jointCount[nextVal][prevVal[c]] ) /
				  		  ( (double) nextCount[nextVal] *
				  			(double) prevCount[prevVal[c]] );
				// Now account for the fact that we've
				//  just used counts rather than probabilities,
				//  and we've got two counts on the bottom
				//  but one count on the top:
				logTerm *= (double) observations;
				localActive[r][c] = Math.log(logTerm) / log_base;
				average += localActive[r][c];
				if (localActive[r][c] > max) {
					max = localActive[r][c];
				} else if (localActive[r][c] < min) {
					min = localActive[r][c];
				}
				// Update the previous value:
				prevVal[c] -= maxShiftedValue[states[r-k][c]];
				prevVal[c] *= base;
				prevVal[c] += states[r][c];
			}
		}
		average = average/(double) (columns * (rows - k));
		
		return localActive;
		
	}
	
	/**
	 * Computes local active information storage for the given
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

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][][] localActive = new double[timeSteps][agentRows][agentColumns];
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
					logTerm = ( (double) jointCount[nextVal][prevVal[r][c]] ) /
					  		  ( (double) nextCount[nextVal] *
					  			(double) prevCount[prevVal[r][c]] );
					// Now account for the fact that we've
					//  just used counts rather than probabilities,
					//  and we've got two counts on the bottom
					//  but one count on the top:
					logTerm *= (double) observations;
					localActive[t][r][c] = Math.log(logTerm) / log_base;
					average += localActive[t][r][c];
					if (localActive[t][r][c] > max) {
						max = localActive[t][r][c];
					} else if (localActive[t][r][c] < min) {
						min = localActive[t][r][c];
					}
					// Update the previous value:
					prevVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
					prevVal[r][c] *= base;
					prevVal[r][c] += states[t][r][c];
				}
			}
		}
		average = average/(double) (agentRows * agentColumns * (timeSteps - k));
		
		return localActive;
	}

	/**
	 * Computes local active information storage for the given
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
		double[] localActive = new double[rows];
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
			logTerm = ( (double) jointCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) prevCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localActive[r] = Math.log(logTerm) / log_base;
			average += localActive[r];
			if (localActive[r] > max) {
				max = localActive[r];
			} else if (localActive[r] < min) {
				min = localActive[r];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][col]];
			prevVal *= base;
			prevVal += states[r][col];
		}
		average = average/(double) (rows - k);
		
		return localActive;
		
	}
	
	/**
	 * Computes local active information storage for the given
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
		double[] localActive = new double[timeSteps];
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
			logTerm = ( (double) jointCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) prevCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localActive[t] = Math.log(logTerm) / log_base;
			average += localActive[t];
			if (localActive[t] > max) {
				max = localActive[t];
			} else if (localActive[t] < min) {
				min = localActive[t];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k][agentIndex1][agentIndex2]];
			prevVal *= base;
			prevVal += states[t][agentIndex1][agentIndex2];
		}
		average = average/(double) (timeSteps - k);
		
		return localActive;
		
	}

	/**
	 * Standalone routine to 
	 * compute local active information storage across an
	 *  array of the states of homogeneous agents
	 * Return an array of local values.
	 * First history rows are zeros
	 * 
	 * @param history - parameter k
	 * @param maxEmbeddingLength - base of the states
	 * @param states - time series array of states
	 * @return
	 */
	public double[] computeLocal(int states[]) {
		
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	/**
	 * Standalone routine to 
	 * compute local active information storage across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * 
	 * @param history - parameter k
	 * @param maxEmbeddingLength - base of the states
	 * @param states - 2D array of states
	 * @return
	 */
	public double[][] computeLocal(int states[][]) {
		
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}
	
	/**
	 * Standalone routine to 
	 * compute local active information storage across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * 
	 * @param history - parameter k
	 * @param maxEmbeddingLength - base of the states
	 * @param states - array of states - 1st dimension is time, 2nd and 3rd are agent indices
	 * @return
	 */
	public double[][][] computeLocal(int states[][][]) {
		
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	/**
	 * Standalone routine to 
	 * compute average local active information storage across an
	 *  array of the states of homogeneous agents
	 * Return the averagen
	 * 
	 * @param history - parameter k
	 * @param maxEmbeddingLength - base of the states
	 * @param states - narray of states
	 * @return
	 */
	public double computeAverageLocal(int states[]) {
		
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average local active information storage across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return the average
	 * This method to be called for homogeneous agents only
	 * 
	 * @param history - parameter k
	 * @param maxEmbeddingLength - base of the states
	 * @param states - 2D array of states
	 * @return
	 */
	public double computeAverageLocal(int states[][]) {
		
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average local active information storage across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return the average
	 * This method to be called for homogeneous agents only
	 * 
	 * @param history - parameter k
	 * @param maxEmbeddingLength - base of the states
	 * @param states - array of states - 1st dimension is time, 2nd and 3rd are agent indices
	 * @return
	 */
	public double computeAverageLocal(int states[][][]) {
		
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute local active information storage for one agent in a 2D spatiotemporal
	 *  array of the states of agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method should be used for heterogeneous agents
	 * 
	 * @param history - parameter k
	 * @param maxEmbeddingLength - base of the states
	 * @param states - 2D array of states
	 * @param col - column number of the agent in the states array
	 * @return
	 */
	public double[] computeLocal(int states[][], int col) {
		
		initialise();
		addObservations(states, col);
		return computeLocalFromPreviousObservations(states, col);
	}
	
	/**
	 * Standalone routine to 
	 * compute local active information storage for one agent in a 2D spatiotemporal
	 *  array of the states of agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method should be used for heterogeneous agents
	 * 
	 * @param history - parameter k
	 * @param maxEmbeddingLength - base of the states
	 * @param states - array of states - 1st dimension is time, 2nd and 3rd are agent indices
	 * @param agentIndex1 row index of agent
	 * @param agentIndex2 column index of agent
	 * @return
	 */
	public double[] computeLocal(int states[][][], int agentIndex1, int agentIndex2) {
		
		initialise();
		addObservations(states, agentIndex1, agentIndex2);
		return computeLocalFromPreviousObservations(states, agentIndex1, agentIndex2);
	}

	/**
	 * Standalone routine to 
	 * compute average local active information storage 
	 * for a single agent
	 * Returns the average
	 * This method suitable for heterogeneous agents
	 * 
	 * @param history - parameter k
	 * @param maxEmbeddingLength - base of the states
	 * @param states - 2D array of states
	 * @param col - column number of the agent in the states array
	 * @return
	 */
	public double computeAverageLocal(int states[][], int col) {
		
		initialise();
		addObservations(states, col);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average local active information storage 
	 * for a single agent
	 * Returns the average
	 * This method suitable for heterogeneous agents
	 * 
	 * @param history - parameter k
	 * @param maxEmbeddingLength - base of the states
	 * @param states - array of states - 1st dimension is time, 2nd and 3rd are agent indices
	 * @param agentIndex1 row index of agent
	 * @param agentIndex2 column index of agent
	 * @return
	 */
	public double computeAverageLocal(int states[][][], int agentIndex1, int agentIndex2) {
		
		initialise();
		addObservations(states, agentIndex1, agentIndex2);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Compute the significance of obtaining the given average from the given observations
	 * 
	 * @param numPermutationsToCheck number of new orderings of the source values to compare against
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) {
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(observations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * Compute the significance of obtaining the given average from the given observations
	 * 
	 * @param newOrderings the reorderings to use
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) {
		double actualMI = computeAverageLocalOfObservations();
		
		int numPermutationsToCheck = newOrderings.length;
		
		// Reconstruct the values of the previous and next variables (not necessarily in order)
		int[] prevValues = new int[observations];
		int[] nextValues = new int[observations];
		int t_prev = 0;
		int t_next = 0;
		for (int prevVal = 0; prevVal < prevCount.length; prevVal++) {
			int numberOfSamplesPrev = prevCount[prevVal];
			MatrixUtils.fill(prevValues, prevVal, t_prev, numberOfSamplesPrev);
			t_prev += numberOfSamplesPrev;
		}
		for (int nextVal = 0; nextVal < base; nextVal++) {
			int numberOfSamplesNext = nextCount[nextVal];
			MatrixUtils.fill(nextValues, nextVal, t_next, numberOfSamplesNext);
			t_next += numberOfSamplesNext;
		}
		
		ActiveInformationCalculator ais2;
		ais2 = new ActiveInformationCalculator(base, k);
		ais2.initialise();
		ais2.observations = observations;
		ais2.prevCount = prevCount;
		ais2.nextCount = nextCount;
		int countWhereMIIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the next variable
			int[] newDataNext = MatrixUtils.extractSelectedTimePoints(nextValues, newOrderings[p]);
			// compute the joint probability distribution
			MatrixUtils.fill(ais2.jointCount, 0);
			for (int t = 0; t < observations; t++) {
				ais2.jointCount[newDataNext[t]][prevValues[t]]++;
			}
			// And get an MI value for this realisation:
			double newMI = ais2.computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newMI;
			if (newMI >= actualMI) {
				countWhereMIIsMoreSignificantThanOriginal++;
			}

		}
		
		// And return the significance
		measDistribution.pValue = (double) countWhereMIIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualMI;
		return measDistribution;
	}

	public double getLastAverage() {
		return average;
	}

	public double getLastMax() {
		return max;
	}

	public double getLastMin() {
		return min;
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
		for (int nextVal = 0; nextVal < base; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) prevCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) jointCount[nextVal][prevVal] / (double) observations;
				// Compute MI contribution:
				if (p_joint * p_next * p_prev > 0.0) {
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
	 * @param x
	 * @param t
	 * @return
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
	 * Utility function to compute the combined past values of x up to and including time step t
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 * 
	 * @param x
	 * @param agentNumber
	 * @param t
	 * @return
	 */
	public int computePastValue(int[][] x, int agentNumber, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += x[t - k + 1 + p][agentNumber];
		}
		return pastVal;
	}

	/**
	 * Utility function to compute the combined past values of x up to and including time step t
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 * 
	 * @param x
	 * @param agentNumber
	 * @param t
	 * @return
	 */
	public int computePastValue(int[][][] x, int agentRow, int agentColumn, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += x[t - k + 1 + p][agentRow][agentColumn];
		}
		return pastVal;
	}
}
