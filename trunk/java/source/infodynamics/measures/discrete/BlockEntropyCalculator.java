package infodynamics.measures.discrete;

import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * Block entropy calculator
 * Implemented separately from the single state entropy 
 *  calculator to allow the single state calculator 
 *  to have optimal performance.
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
public class BlockEntropyCalculator extends EntropyCalculator {

	protected int blocksize = 0; // temporal blocksize to compute entropy over. Need initialised to 0 for changedSizes
	protected int[] maxShiftedValue = null; // states * (base^(blocksize-1))

	protected int base_power_blocksize = 0;

	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param blocksize
	 * @param base
	 * @return
	 */
	public static EntropyCalculator newInstance(int blocksize, int base) {
		return new BlockEntropyCalculator(blocksize, base);
	}
	
	/**
	 * 
	 * @param blocksize
	 * @param base
	 */
	public BlockEntropyCalculator(int blocksize, int base) {

		super(base);
		
		this.blocksize = blocksize;
		base_power_blocksize = MathsUtils.power(base, blocksize);

		if (blocksize <= 1) {
			throw new RuntimeException("Blocksize " + blocksize + " is not > 1 for Block Entropy Calculator");
		}
		if (blocksize > Math.log(Integer.MAX_VALUE) / log_base) {
			throw new RuntimeException("Base and blocksize combination too large");
		}

		// Create storage for counts of observations.
		// We're recreating stateCount, which was created in super()
		//  however the super() didn't create it for blocksize > 1
		stateCount = new int[MathsUtils.power(base, blocksize)];

		// Create constants for tracking stateValues
		maxShiftedValue = new int[base];
		for (int v = 0; v < base; v++) {
			maxShiftedValue[v] = v * MathsUtils.power(base, blocksize-1);
		}
	}
	
	/**
	 * Initialise calculator, preparing to take observation sets in
	 * Should be called prior to any of the addObservations() methods.
	 * You can reinitialise without needing to create a new object.
	 */
	public void initialise(){
		super.initialise();
		MatrixUtils.fill(stateCount, 0);
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
		observations += (rows - blocksize + 1)*columns; 
		
		// Initialise and store the current previous value for each column
		int[] stateVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			stateVal[c] = 0;
			for (int p = 0; p < blocksize - 1; p++) {
				// Add the contribution from this observation
				stateVal[c] += states[p][c];
				// And shift up
				stateVal[c] *= base;
			}
		}

		// 1. Now count the tuples observed from the next row onwards
		for (int r = blocksize - 1; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				// Add the contribution from this observation
				stateVal[c] += states[r][c];

				// Add to the count for this particular state:
				stateCount[stateVal[c]]++;
				
				// Remove the oldest observation from the state value
				stateVal[c] -= maxShiftedValue[states[r-blocksize+1][c]];
				stateVal[c] *= base;
			}
		}		
	}
	
	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states 1st index is time, 2nd and 3rd index are agent number
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
		observations += (timeSteps - blocksize + 1) * agentRows * agentColumns; 
		
		// Initialise and store the current previous value for each column
		int[][] stateVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				stateVal[r][c] = 0;
				for (int p = 0; p < blocksize - 1; p++) {
					// Add the contribution from this observation
					stateVal[r][c] += states[p][r][c];
					// And shift up
					stateVal[r][c] *= base;
				}
			}
		}

		// 1. Now count the tuples observed from the next row onwards
		for (int t = blocksize - 1; t < timeSteps; t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					// Add the contribution from this observation
					stateVal[r][c] += states[t][r][c];
	
					// Add to the count for this particular state:
					stateCount[stateVal[r][c]]++;
					
					// Remove the oldest observation from the state value
					stateVal[r][c] -= maxShiftedValue[states[t-blocksize+1][r][c]];
					stateVal[r][c] *= base;
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
		observations += rows - blocksize + 1; 
		
		// Initialise and store the current previous value for each column
		int stateVal = 0;
		for (int p = 0; p < blocksize - 1; p++) {
			// Add the contribution from this observation
			stateVal += states[p][col];
			// And shift up
			stateVal *= base;
		}

		// 1. Count the tuples observed
		for (int r = blocksize - 1; r < rows; r++) {
			// Add the contribution from this observation
			stateVal += states[r][col];

			// Add to the count for this particular state:
			stateCount[stateVal]++;
			
			// Remove the oldest observation from the state value
			stateVal -= maxShiftedValue[states[r-blocksize+1][col]];
			stateVal *= base;
		}
	}

	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 */
	public void addObservations(int states[][][], int agentIndex1, int agentIndex2) {
		int timeSteps = states.length;
		// increment the count of observations:
		observations += timeSteps - blocksize + 1; 
		
		// Initialise and store the current previous value for each column
		int stateVal = 0;
		for (int p = 0; p < blocksize - 1; p++) {
			// Add the contribution from this observation
			stateVal += states[p][agentIndex1][agentIndex2];
			// And shift up
			stateVal *= base;
		}

		// 1. Count the tuples observed
		for (int t = blocksize - 1; t < timeSteps; t++) {
			// Add the contribution from this observation
			stateVal += states[t][agentIndex1][agentIndex2];

			// Add to the count for this particular state:
			stateCount[stateVal]++;
			
			// Remove the oldest observation from the state value
			stateVal -= maxShiftedValue[states[t-blocksize+1][agentIndex1][agentIndex2]];
			stateVal *= base;
		}
	}

	/**
	 * Returns the average local entropy from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return
	 */
	public double computeAverageLocalOfObservations() {
		double ent = 0.0;
		double entCont = 0.0;

		max = 0;
		min = 0;
		for (int stateVal = 0; stateVal < base_power_blocksize; stateVal++) {
			// compute p_state
			double p_state = (double) stateCount[stateVal] / (double) observations;
			if (p_state > 0.0) {
				// Entropy takes the negative log:
				double localValue = - Math.log(p_state) / log_2;
				entCont = p_state * localValue;
				if (localValue > max) {
					max = localValue;
				} else if (localValue < min) {
					min = localValue;
				}
			} else {
				entCont = 0.0;
			}
			ent += entCont;
		}
		
		average = ent;
		return ent;
	}
	
	/**
	 * Computes local entropy for the given
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
		
		// Allocate for all rows even though we may not be assigning all of them
		double[][] localEntropy = new double[rows][columns];
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous value for each column
		int[] stateVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			stateVal[c] = 0;
			for (int p = 0; p < blocksize - 1; p++) {
				stateVal[c] += states[p][c];
				stateVal[c] *= base;
			}
		}
		// StateVal just needs the next value put in before processing

		for (int r = blocksize - 1; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				// Add the most recent observation to the state count
				stateVal[c] += states[r][c];
				// Process the stateVal:
				double p_state = (double) stateCount[stateVal[c]] / (double) observations;
				// Entropy takes the negative log:
				localEntropy[r][c] = - Math.log(p_state) / log_2;
				average += localEntropy[r][c];
				if (localEntropy[r][c] > max) {
					max = localEntropy[r][c];
				} else if (localEntropy[r][c] < min) {
					min = localEntropy[r][c];
				}
				// Subtract out the oldest part of the state value:
				stateVal[c] -= maxShiftedValue[states[r-blocksize+1][c]];
				// And shift all upwards
				stateVal[c] *= base;
			}
		}
		average = average/(double) (columns * (rows - blocksize + 1));
		
		return localEntropy;
		
	}
	
	/**
	 * Computes local entropy for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method to be used for homogeneous agents only
	 *  
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 * @return
	 */
	public double[][][] computeLocalFromPreviousObservations(int states[][][]){
		int timeSteps = states.length;
		int agentRows, agentColumns;
		if (timeSteps == 0) {
			agentRows = 0;
			agentColumns = 0;
		} else {
			agentRows = states[0].length;
			if (agentRows == 0) {
				agentColumns = 0;
			} else {
				agentColumns = states[0][0].length;
			}
		}
		
		// Allocate for all rows even though we may not be assigning all of them
		double[][][] localEntropy = new double[timeSteps][agentRows][agentColumns];
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous value for each column
		int[][] stateVal = new int[agentRows][agentColumns]; 
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				stateVal[r][c] = 0;
				for (int p = 0; p < blocksize - 1; p++) {
					stateVal[r][c] += states[p][r][c];
					stateVal[r][c] *= base;
				}
			}
		}
		// StateVal just needs the next value put in before processing

		for (int t = blocksize - 1; t < timeSteps; t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					// Add the most recent observation to the state count
					stateVal[r][c] += states[t][r][c];
					// Process the stateVal:
					double p_state = (double) stateCount[stateVal[r][c]] / (double) observations;
					// Entropy takes the negative log:
					localEntropy[t][r][c] = - Math.log(p_state) / log_2;
					average += localEntropy[t][r][c];
					if (localEntropy[t][r][c] > max) {
						max = localEntropy[t][r][c];
					} else if (localEntropy[t][r][c] < min) {
						min = localEntropy[t][r][c];
					}
					// Subtract out the oldest part of the state value:
					stateVal[r][c] -= maxShiftedValue[states[t-blocksize+1][r][c]];
					// And shift all upwards
					stateVal[r][c] *= base;
				}
			}
		}
		average = average/(double) (agentRows * agentColumns * (timeSteps - blocksize + 1));
		
		return localEntropy;
		
	}

	/**
	 * Computes local entropy for the given
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
		// First row to consider the destination cell at:
		int startRow = blocksize - 1;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localEntropy = new double[rows];
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous value for each column
		int stateVal = 0; 
		for (int p = 0; p < blocksize - 1; p++) {
			stateVal += states[p][col];
			stateVal *= base;
		}
		// StateVal just needs the next value put in before processing

		for (int r = startRow; r < rows; r++) {
			// Add the most recent observation to the state count
			stateVal += states[r][col];
			// Process the stateVal:
			double p_state = (double) stateCount[stateVal] / (double) observations;
			// Entropy takes the negative log:
			localEntropy[r] = - Math.log(p_state) / log_2;
			average += localEntropy[r];
			if (localEntropy[r] > max) {
				max = localEntropy[r];
			} else if (localEntropy[r] < min) {
				min = localEntropy[r];
			}
			// Subtract out the oldest part of the state value:
			stateVal -= maxShiftedValue[states[r-blocksize+1][col]];
			// And shift all upwards
			stateVal *= base;
		}
		average = average/(double) (rows - blocksize + 1);
		
		return localEntropy;
		
	}
	
	/**
	 * Computes local entropy for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][][], int agentIndex1, int agentIndex2){
		int timeSteps = states.length;
		//int columns = states[0].length;
		// First time to consider the destination cell at:
		int startTime = blocksize - 1;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localEntropy = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous value for each column
		int stateVal = 0; 
		for (int p = 0; p < blocksize - 1; p++) {
			stateVal += states[p][agentIndex1][agentIndex2];
			stateVal *= base;
		}
		// StateVal just needs the next value put in before processing

		for (int t = startTime; t < timeSteps; t++) {
			// Add the most recent observation to the state count
			stateVal += states[t][agentIndex1][agentIndex2];
			// Process the stateVal:
			double p_state = (double) stateCount[stateVal] / (double) observations;
			// Entropy takes the negative log:
			localEntropy[t] = - Math.log(p_state) / log_2;
			average += localEntropy[t];
			if (localEntropy[t] > max) {
				max = localEntropy[t];
			} else if (localEntropy[t] < min) {
				min = localEntropy[t];
			}
			// Subtract out the oldest part of the state value:
			stateVal -= maxShiftedValue[states[t-blocksize+1][agentIndex1][agentIndex2]];
			// And shift all upwards
			stateVal *= base;
		}
		average = average/(double) (timeSteps - blocksize + 1);
		
		return localEntropy;
		
	}
}
