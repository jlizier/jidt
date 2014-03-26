package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;

/**
 * <p>Univariate entropy calculator</p>
 * 
 * <p>Usage:
 * <ol>
 *  <li>Continuous accumulation of observations. Call:
 *  	<ol>
 *  		<li>{@link #initialise()};</li>
 *          <li>then supply the observations using any of the addObservations()
 *          methods (several times over), e.g. {@link #addObservations(int[][])};</li>
 *          <li>then when all observations have been added, call any of the
 *          compute methods (several times over), e.g.: 
 *          {@link #computeLocalFromPreviousObservations(int[][])}</li>
 *      </ol>
 *   </li>
 *   <li>Standalone mode. Call:
 *   	<ol>
 *   		<li>one of the standalone methods which supply observations
 *          and compute at once, e.g. {@link #computeAverageLocal(int[][])}.</li>
 *      </ol>
 *   </li>
 * </ol></p>
 * 
 * 
 * @author Joseph Lizier
 */
public class EntropyCalculator extends InfoMeasureCalculator
				implements SingleAgentMeasure
{

	protected int[] stateCount = null; // Count for i[t]

	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base
	 * @param blocksize
	 * 
	 * @return
	 */
	public static EntropyCalculator newInstance(int base, int blocksize) {
		if (blocksize > 1) {
			return BlockEntropyCalculator.newInstance(blocksize, base);
		} else {
			return EntropyCalculator.newInstance(base);
		}
	}
	public static EntropyCalculator newInstance(int base) {
		return new EntropyCalculator(base);
	}
	
	/**
	 * 
	 * @param base
	 */
	public EntropyCalculator(int base) {

		super(base);
		
		// Create storage for counts of observations
		stateCount = new int[base];
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
	 *
	 * @param states 1st index is time
	 */
	public void addObservations(int states[]) {
		int rows = states.length;
		// increment the count of observations:
		observations += rows; 
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			// Add to the count for this particular state:
			stateCount[states[r]]++;
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
		observations += rows * columns; 
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				// Add to the count for this particular state:
				stateCount[states[r][c]]++;					
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
		observations += timeSteps * agentRows * agentColumns; 
		
		// 1. Count the tuples observed
		for (int t = 0; t < timeSteps; t++) {
			for (int i = 0; i < agentRows; i++) {
				for (int j = 0; j < agentColumns; j++) {
					// Add to the count for this particular state:
					stateCount[states[t][i][j]]++;
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
	 * @param agentNumber
	 */
	public void addObservations(int states[][], int agentNumber) {
		int rows = states.length;
		// increment the count of observations:
		observations += rows; 
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			// Add to the count for this particular state:
			stateCount[states[r][agentNumber]]++;					
		}
	}

	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param agentIndex1
	 * @param agentIndex2
	 */
	public void addObservations(int states[][][], int agentIndex1, int agentIndex2) {
		int timeSteps = states.length;
		// increment the count of observations:
		observations += timeSteps; 
		
		// 1. Count the tuples observed
		for (int r = 0; r < timeSteps; r++) {
			// Add to the count for this particular state:
			stateCount[states[r][agentIndex1][agentIndex2]]++;					
		}
	}

	/**
	 * 
	 * @param stateVal
	 * @return count of observations of the given state
	 */
	public int getStateCount(int stateVal) {
		return stateCount[stateVal];
	}
	
	/**
	 * 
	 * @param stateVal
	 * @return probability of the given state
	 */
	public double getStateProbability(int stateVal) {
		return (double) stateCount[stateVal] / (double) observations;
	}

	/**
	 * Returns the average entropy from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return the average entropy
	 */
	public double computeAverageLocalOfObservations() {
		double ent = 0.0;
		double entCont = 0.0;

		max = 0;
		min = 0;
		for (int stateVal = 0; stateVal < base; stateVal++) {
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
	 *  
	 * @param states index is time
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[]){
		int rows = states.length;
		
		double[] localEntropy = new double[rows];
		average = 0;
		max = 0;
		min = 0;
		for (int r = 0; r < rows; r++) {
			double p_state = (double) stateCount[states[r]] / (double) observations;
			// Entropy takes the negative log:
			localEntropy[r] = - Math.log(p_state) / log_2;
			average += localEntropy[r];
			if (localEntropy[r] > max) {
				max = localEntropy[r];
			} else if (localEntropy[r] < min) {
				min = localEntropy[r];
			}
		}
		average = average/(double) rows;
		
		return localEntropy;
		
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
		
		double[][] localEntropy = new double[rows][columns];
		average = 0;
		max = 0;
		min = 0;
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				double p_state = (double) stateCount[states[r][c]] / (double) observations;
				// Entropy takes the negative log:
				localEntropy[r][c] = - Math.log(p_state) / log_2;
				average += localEntropy[r][c];
				if (localEntropy[r][c] > max) {
					max = localEntropy[r][c];
				} else if (localEntropy[r][c] < min) {
					min = localEntropy[r][c];
				}
			}
		}
		average = average/(double) (columns * rows);
		
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

		double[][][] localEntropy = new double[timeSteps][agentRows][agentColumns];
		average = 0;
		max = 0;
		min = 0;
		for (int r = 0; r < timeSteps; r++) {
			for (int i = 0; i < agentRows; i++) {
				for (int j = 0; j < agentColumns; j++) {
					double p_state = (double) stateCount[states[r][i][j]] / (double) observations;
					// Entropy takes the negative log:
					localEntropy[r][i][j] = - Math.log(p_state) / log_2;
					average += localEntropy[r][i][j];
					if (localEntropy[r][i][j] > max) {
						max = localEntropy[r][i][j];
					} else if (localEntropy[r][i][j] < min) {
						min = localEntropy[r][i][j];
					}
				}
			}
		}
		average = average/(double) (agentRows * agentColumns * timeSteps);
		
		return localEntropy;
		
	}

	/**
	 * Computes local entropy for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states
	 * @param agentNumber
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][], int agentNumber){
		int rows = states.length;
		//int columns = states[0].length;

		double[] localEntropy = new double[rows];
		average = 0;
		max = 0;
		min = 0;
		for (int r = 0; r < rows; r++) {
			double p_state = (double) stateCount[states[r][agentNumber]] / (double) observations;
			// Entropy takes the negative log:
			localEntropy[r] = - Math.log(p_state) / log_2;
			average += localEntropy[r];
			if (localEntropy[r] > max) {
				max = localEntropy[r];
			} else if (localEntropy[r] < min) {
				min = localEntropy[r];
			}
		}
		average = average/(double) (rows);
		
		return localEntropy;
		
	}
	
	/**
	 * Computes local entropy for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 * @param agentIndex1
	 * @param agentIndex2
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][][], int agentIndex1, int agentIndex2){
		int timeSteps = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localEntropy = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;
		for (int r = 0; r < timeSteps; r++) {
			double p_state = (double) stateCount[states[r][agentIndex1][agentIndex2]] / (double) observations;
			// Entropy takes the negative log:
			localEntropy[r] = - Math.log(p_state) / log_2;
			average += localEntropy[r];
			if (localEntropy[r] > max) {
				max = localEntropy[r];
			} else if (localEntropy[r] < min) {
				min = localEntropy[r];
			}
		}
		average = average/(double) (timeSteps);
		
		return localEntropy;
		
	}

	/**
	 * Standalone routine to 
	 * compute local information theoretic measure across a 1D temporal
	 *  array of the states of homogeneous agents
	 * Return a 1D temporal array of local values.
	 * First history rows are zeros
	 * 
	 * @param states - 1D array of states
	 * @return
	 */
	public final double[] computeLocal(int states[]) {
		
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	/**
	 * Standalone routine to 
	 * compute local information theoretic measure across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * 
	 * @param states - 2D array of states
	 * @return
	 */
	public final double[][] computeLocal(int states[][]) {
		
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	/**
	 * Standalone routine to 
	 * compute local information theoretic measure across a 3D spatiotemporal
	 *  array of the states of 2D homogeneous agents
	 * Return a 3D spatiotemporal array of local values.
	 * First history rows are zeros
	 * 
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 * @return
	 */
	public final double[][][] computeLocal(int states[][][]) {
		
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	/**
	 * Standalone routine to 
	 * compute average local information theoretic measure across a 1D temporal
	 *  array of the states of homogeneous agents
	 * Return the average
	 * 
	 * @param states - 1D array of states
	 * @return
	 */
	public final double computeAverageLocal(int states[]) {
		
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average local information theoretic measure across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return the average
	 * This method to be called for homogeneous agents only
	 * 
	 * @param states - 2D array of states
	 * @return
	 */
	public final double computeAverageLocal(int states[][]) {
		
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average local information theoretic measure across a 3D spatiotemporal
	 *  array of the states of 2D homogeneous agents
	 * Return the average
	 * This method to be called for homogeneous agents only
	 * 
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 * @return
	 */
	public final double computeAverageLocal(int states[][][]) {
		
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute local information theoretic measure for one agent in a 2D spatiotemporal
	 *  array of the states of agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method should be used for heterogeneous agents
	 * 
	 * @param states - 2D array of states
	 * @param col - column number of the agent in the states array
	 * @return
	 */
	public final double[] computeLocalAtAgent(int states[][], int col) {
		
		initialise();
		addObservations(states, col);
		return computeLocalFromPreviousObservations(states, col);
	}

	/**
	 * Standalone routine to 
	 * compute local information theoretic measure for one agent in a 2D spatiotemporal
	 *  array of the states of agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method should be used for heterogeneous agents
	 * 
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 * @param agentIndex1
	 * @param agentIndex2
	 * @return
	 */
	public final double[] computeLocalAtAgent(int states[][][], int agentIndex1, int agentIndex2) {
		
		initialise();
		addObservations(states, agentIndex1, agentIndex2);
		return computeLocalFromPreviousObservations(states, agentIndex1, agentIndex2);
	}

	/**
	 * Standalone routine to 
	 * compute average local information theoretic measure 
	 * for a single agent
	 * Returns the average
	 * This method suitable for heterogeneous agents
	 * @param states - 2D array of states
	 * @param col - column number of the agent in the states array
	 * 
	 * @return
	 */
	public final double computeAverageLocalAtAgent(int states[][], int col) {
		initialise();
		addObservations(states, col);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average local information theoretic measure 
	 * for a single agent
	 * Returns the average
	 * This method suitable for heterogeneous agents
	 * @param states - 2D array of states
	 * @param agentIndex1
	 * @param agentIndex2
	 * 
	 * @return
	 */
	public final double computeAverageLocalAtAgent(int states[][][], int agentIndex1, int agentIndex2) {
		initialise();
		addObservations(states, agentIndex1, agentIndex2);
		return computeAverageLocalOfObservations();
	}
}
