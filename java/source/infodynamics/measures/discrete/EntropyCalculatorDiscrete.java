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

import java.util.Hashtable;

import infodynamics.utils.MatrixUtils;

/**
 * <p>Entropy calculator for univariate discrete (int[]) data.</p>
 * 
 * <p>Usage of the class is intended to follow this paradigm:</p>
 * <ol>
 * 		<li>Construct the calculator: {@link #EntropyCalculatorDiscrete(int)};</li>
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
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991).</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EntropyCalculatorDiscrete extends InfoMeasureCalculatorDiscrete
				implements SingleAgentMeasureDiscrete
{

	protected int[] stateCount = null; // Count for i[t]

	protected Hashtable<Object, Integer> hashedStateCount = null;
	
	/**
	 * Construct a new instance with no specified alphabet size.
	 */
	public EntropyCalculatorDiscrete() {
		this(-1);
	}
	
	/**
	 * Contruct a new instance
	 * 
	 * @param alphabetSize number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 */
	public EntropyCalculatorDiscrete(int alphabetSize) {
		super(alphabetSize);
	}
	
	/**
	 * Initialise with new alphabet size
	 * 
	 * @param alphabetSize
	 */
	public void initialise(int alphabetSize){
		boolean sizeChange = (this.alphabetSize != alphabetSize);
		super.initialise(alphabetSize);
		
		if (sizeChange || (stateCount == null)) {
			// Create storage for counts of observations
			try {
				stateCount = new int[alphabetSize];
			} catch (OutOfMemoryError e) {
				// Allow any Exceptions to be thrown, but catch and wrap
				//  Error as a RuntimeException
				throw new RuntimeException("Requested memory for the alphabet size (" +
						alphabetSize + ") is too large for the JVM at this time", e);
			}
		} else {
			MatrixUtils.fill(stateCount, 0);
		}
	}

	@Override
	public void initialise(){
		initialise(alphabetSize);
	}

	@Override
	public void addObservations(Object[] states) {
		int rows = states.length;
		// increment the count of observations:
		observations += rows;
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			// Add to the count for this particular state:

			if (!knownIntegerRange) {
				Object key = states[r];
				Integer value = hashedStateCount.getOrDefault(key, 0) + 1;
				hashedStateCount.put(key, value);
				continue;
			}

			// TODO better way to do this??
			// obviously we can do .getClass() twice in the if statements, I think this is
			// slightly more readable though.
			//	joe's input:
			@SuppressWarnings("rawtypes")
			Class objectClass = states[r].getClass();
			Integer index;
			if (objectClass == String.class) {
				index = Integer.parseInt((String) states[r]);
			} else if (objectClass == Integer.class) {
				index = (Integer) states[r];
			} else {
				throw new NumberFormatException(String.format(
					"Cannot parse %s as an Integer for indexing",
					states[r]));
			}

			// valid index, check to make sure user didn't break the alphabet size they provided.
			if (index >= stateCount.length) {
				// runtime exception I guess...?
				// TODO joe's input:
				throw new RuntimeException(String.format(
					"Observation %s exceeds provided alphabet size %d",
					states[r], this.alphabetSize));
			}
			stateCount[index]++;
		}		
	}

	// I added this so that computeAverageLocal works, which we were considering deprecating
	// Though i think keeping the old int[] format is worthwhile for backwards compat.
	// TODO joe's input:
	// maybe deprecate (with computeAverageLocal)
	/**
	 * Add observations in to our estimates of the pdfs.
	 * 
	 * @param states
	 */
	public void addObservations(int[] states) {
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
	 * Return the current count for the given value
	 * 
	 * @param stateVal given value
	 * @return count of observations of the given state
	 */
	public int getStateCount(int stateVal) {
		return stateCount[stateVal];
	}
	
	/**
	 * Return the current probability for the given value
	 * 
	 * @param stateVal given value
	 * @return probability of the given state
	 */
	public double getStateProbability(int stateVal) {
		return (double) stateCount[stateVal] / (double) observations;
	}

	@Override
	public double computeAverageLocalOfObservations() {
		double ent = 0.0;
		double entCont = 0.0;

		max = 0;
		min = 0;
		for (int stateVal = 0; stateVal < alphabetSize; stateVal++) {
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
	
	@Override
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
	
	@Override // maybe deprecate
	public final double computeAverageLocal(int states[]) throws Exception {
		initialise();
		startAddObservations();
		addObservations(states);
		finaliseAddObservations();
		return computeAverageLocalOfObservations();
	}

	
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		// TODO - I think this should just be the super call, as this class has no new properties.
		//	joe's input: 
		
		switch(propertyName.toLowerCase()) {
			case ALPHABET_SIZE:
				this.alphabetSize = Integer.parseInt(propertyValue);
			break;
			case MAX_ALPHA_SIZE_TO_STORE:
				this.maxAlphabetSize = Integer.parseInt(propertyValue);
			break;
			case KNOWN_INTEGER_RANGE:
				this.knownIntegerRange = Boolean.parseBoolean(propertyValue);
			break;
			default:
				// Assume it's a property of parent class
				super.setProperty(propertyName, propertyValue);
			break;
		}
		
	}

	@Override
	public void setObservations(Object[] observations) throws Exception {
		this.observations = observations.length;
		this.addObservations(observations);
		this.finaliseAddObservations();
	}

	/*
	 ***********************************************************************
	 ******************* ALL METHODS BELOW ARE DEPRECATED ******************
	 *******************        DON'T USE THEM            ******************
	 ***********************************************************************
	*/

	@Override
	@Deprecated
	public void startAddObservations() {
		
		// reinitialise if already finalised
		if (currentState == State.COMPUTING) {
			initialise();
		}
		this.currentState = State.ADDING_OBSERVATIONS;

	}

	@Override
	@Deprecated
	public void finaliseAddObservations() throws Exception {

		if (observations == 0) {
			// TODO better error message. maybe not runtimeExc?
			// joe's input:
			throw new RuntimeException("0 observations not allowed.");
		}

		// I actually think this is useless? It shouldn't ever get here without
		// throwing the error above.
		if (currentState == State.SETTING_PROPERTIES) {
			throw new RuntimeException("Estimator should be initialised before finalised...");
		}

		this.currentState = State.COMPUTING;
	}

	@Override
	@Deprecated
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
	
	@Override
	@Deprecated
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

	@Override
	@Deprecated
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

	@Override
	@Deprecated
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
	
	@Override
	@Deprecated
	public final double computeAverageLocal(int states[][]) {
		
		initialise();
		startAddObservations();
		addObservations(states);
		this.currentState = State.COMPUTING;
		return computeAverageLocalOfObservations();
	}
	
	@Override
	@Deprecated
	public final double computeAverageLocal(int states[][][]) {
		initialise();
		startAddObservations();
		addObservations(states);
		this.currentState = State.COMPUTING;
		return computeAverageLocalOfObservations();
	}
	
	@Override
	@Deprecated
	public final double[] computeLocal(int states[][], int col) {
		initialise();
		startAddObservations();
		addObservations(states, col);
		this.currentState = State.COMPUTING;
		return computeLocalFromPreviousObservations(states, col);
	}
	
	@Override
	@Deprecated
	public final double[] computeLocal(int states[][][],
			int agentIndex1, int agentIndex2) {
		initialise();
		startAddObservations();
		addObservations(states, agentIndex1, agentIndex2);
		this.currentState = State.COMPUTING;
		return computeLocalFromPreviousObservations(states, agentIndex1, agentIndex2);
	}
	
	@Override
	@Deprecated
	public final double computeAverageLocal(int states[][], int col) {
		initialise();
		startAddObservations();
		addObservations(states, col);
		this.currentState = State.COMPUTING;
		return computeAverageLocalOfObservations();
	}
	
	@Override
	@Deprecated
	public final double computeAverageLocal(int states[][][], int agentIndex1, int agentIndex2) {
		initialise();
		startAddObservations();
		addObservations(states, agentIndex1, agentIndex2);
		this.currentState = State.COMPUTING;
		return computeAverageLocalOfObservations();
	}
	
	@Override
	@Deprecated
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
	
	@Override
	@Deprecated
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
	
	@Override
	@Deprecated
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
	
	@Override
	@Deprecated
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
	
	@Override
	@Deprecated
	public final double[] computeLocal(int states[]) {
		initialise();
		startAddObservations();
		addObservations(states);
		this.currentState = State.COMPUTING;
		return computeLocalFromPreviousObservations(states);
	}
	
	@Override
	@Deprecated
	public final double[][] computeLocal(int states[][]) {
		initialise();
		startAddObservations();
		addObservations(states);
		this.currentState = State.COMPUTING;
		return computeLocalFromPreviousObservations(states);
	}
	
	@Override
	@Deprecated
	public final double[][][] computeLocal(int states[][][]) {
		initialise();
		startAddObservations();
		addObservations(states);
		this.currentState = State.COMPUTING;
		return computeLocalFromPreviousObservations(states);
	}

	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param alphabetSize number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param blocksize number of consecutive joint values to include
	 *  in the calculation.
	 * @deprecated
	 * @return a new EntropyCalculator
	 */
	@Deprecated
	public static EntropyCalculatorDiscrete newInstance(int alphabetSize, int blocksize) {
		if (blocksize > 1) {
			return BlockEntropyCalculatorDiscrete.newInstance(blocksize, alphabetSize);
		} else {
			return EntropyCalculatorDiscrete.newInstance(alphabetSize);
		}
	}
	@Deprecated
	public static EntropyCalculatorDiscrete newInstance(int alphabetSize) {
		return new EntropyCalculatorDiscrete(alphabetSize);
	}

}
