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
import java.util.Map;
import java.util.List;
import java.util.Arrays;

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
				implements UnivariateMeasureDiscrete
{

	// TODO -- find or calculate reasonable max value. Currently 100.
	protected int[] stateCount = null; // Count for i[t]

	/**
	 * State counts in a hashtable for observations other than ints
	 * as well as for sparse observations
	 */
	protected Hashtable<Object, Integer> hashedStateCount = null;

	/**
	 * Number of dimensions for multi-dimensional calculations
	 */
	protected int numDimensions = 1;

	/**
	 * Sizes of each dimension's alphabet
	 */
	protected int[] alphabetSizes = null;
	
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
		currentState = State.SETTING_PROPERTIES;
	}
	
	/**
	 * Return the current count for the given value
	 * 
	 * @param stateVal given value
	 * @return count of observations of the given state
	 */
	public int getStateCount(int stateVal) {
		// stateCount could have overflowed
		if (stateCount == null) {
			return hashedStateCount.get((Integer) stateVal);
		}
		return stateCount[stateVal];
	}

	/**
	 * Return the current count for the given value
	 * 
	 * @param stateVal given value
	 * @return count of observations of the given state
	 */
	public int getStateCount(Object stateVal) {
		if (stateCount == null) {
			return hashedStateCount.get(stateVal);
		} else if (stateVal instanceof String) {
			int index = Integer.parseInt((String) stateVal);
			return stateCount[index];
		} else {
			int index = (Integer) stateVal;
			return stateCount[index];
		}
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
	
	/**
	 * Initialise with new alphabet size
	 * 
	 * @param alphabetSize
	 */
	public void initialise(int alphabetSize){

		boolean sizeChange = (this.alphabetSize != alphabetSize);
		super.initialise(alphabetSize);

		// if alphabet size is unknown, use large but managable array
		// only moving to hash table for arrays exceeding this temp limit.
		int temp = alphabetSize;
		if (temp == -1) {
			// TODO -- This is where we set the reasonable max value. See top of class.
			temp = 100;
		}

		if (sizeChange || stateCount == null) {
			// Create storage for counts of observations
			try {
				stateCount = new int[temp];
			} catch (OutOfMemoryError e) {
				// Allow any Exceptions to be thrown, but catch and wrap
				//  Error as a RuntimeException
				throw new RuntimeException("Requested memory for the alphabet size (" +
						alphabetSize + ") is too large for the JVM at this time", e);
			}
		} else {
			MatrixUtils.fill(stateCount, 0);
		}
		hashedStateCount = new Hashtable<>();
		currentState = State.INITIALISED;
	}

	@Override
	public void initialise(){
		initialise(alphabetSize);
	}
	
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		// TODO -- set states here too, enforce correct workflow
		// I believe this has been done
		if (currentState != State.SETTING_PROPERTIES) {
			stateCount = null;
			hashedStateCount = new Hashtable<>();
			knownIntegerRange = false;
			alphabetSizes = null;
			numDimensions = 1;
			currentState = State.SETTING_PROPERTIES;
			
			// TODO -- This is where we set the reasonable max value. See top of class.
			alphabetSize = 100;
		}

		switch(propertyName.toUpperCase()) {
			case "NUM_DIMENSIONS":
				try {
					this.numDimensions = Integer.parseInt(propertyValue);
					if (numDimensions < 0) {
						throw new NumberFormatException("propertyValue for NUM_DIMENSIONS must be at least 0.");
					}
				} catch (NumberFormatException e) {
					throw new NumberFormatException("propertyValue for NUM_DIMENSIONS must be interpretable as an integer.");
				}
			break;
			case "ALPHABET_SIZES":
				/*
				 * Currently accepts '1,2,3' or '[1,2,3]' formats, as this seemed
				 * easiest for importing data as well as manual input.
				 * '1, 2, 3' and '[1, 2, 3]' should also work with parseInt, though I'm unsure.
				 */
				propertyValue.replaceFirst("[", "");
				propertyValue.replaceFirst("]", "");
				String[] temp = propertyValue.split(",");
				this.alphabetSizes = new int[temp.length];
				for (int i = 0; i < temp.length; i++) {
					try {
						alphabetSizes[i] = Integer.parseInt(temp[i]);
					} catch (NumberFormatException e) {
						throw new NumberFormatException("propertyValue for ALPHABET_SIZES must be interpretable as an integer array.\n"
						+ "Input should look like [1,2,3]");
					}
				} 
			break;
			default:
				super.setProperty(propertyName, propertyValue);
			break;
		}
	}

	public void addObservations(Object[] states) throws NumberFormatException, RuntimeException {
		// if user has set numDimensions, they intend for multi-dimensional, and
		// should not be using this function
		if (numDimensions != 1) {
			throw new RuntimeException(
			"numDimensions was not 1. " +
			"If you intend to use multi-dimensional observations, use addObservations(Object[][] states)."
			);
		}

		if (currentState == State.COMPUTING ||
			currentState == State.SETTING_PROPERTIES) {
				initialise();
		}
		currentState = State.ADDING_OBSERVATIONS;
		
		int rows = states.length;
		// increment the count of observations:
		observations += rows;
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			// Add to the count for this particular state:

			if (!hashedStateCount.isEmpty()) {
				Object key = states[r];
				Integer value = hashedStateCount.getOrDefault(key, 0) + 1;
				hashedStateCount.put(key, value);
				continue;
			}
			
			Integer index = this.alphabetSize;
			if (states[r] instanceof String ) {
				if (((String) states[r]).matches("[0-9]+")) {
					index = Integer.parseInt((String) states[r]);
				} else if (!knownIntegerRange) {
					// First observation ever is non-integer string, right to hashtable
					hashedStateCount.put((Object) states[r], 1);
					stateCount = null;
					continue;
				}
			} else if (states[r] instanceof Integer) {
				index = (Integer) states[r];
			} else {
				throw new NumberFormatException(
					"Cannot parse " + states[r] + " as an Integer for indexing.");
			}

			// valid index, check to make sure user didn't break the alphabet size they provided.
			if (knownIntegerRange && index >= alphabetSize) {
				throw new RuntimeException("Observation " + states[r] + 
				" exceeds provided alphabet size " + alphabetSize);
			}
			
			// if array is still a managable size...
			if (index < stateCount.length) {
				stateCount[index]++;
			} else {
				// otherwise array has grown too large, move to hash table.
				for (int i = 0; i < stateCount.length; i++) {
					hashedStateCount.put(i, stateCount[i]);
				}
				hashedStateCount.put(index, 1);
				stateCount = null;
			}
		}		
	}

	/**
	 * Add observations in to our estimates of the pdfs.
	 * 
	 * @param states
	 */
	public void addObservations(int[] states) {

		// if user has set numDimensions, they intend for multi-dimensional, and
		// should not be using this function
		if (numDimensions != 1) {
			throw new RuntimeException(
			"numDimensions was not 1. " +
			"If you intend to use multi-dimensional observations, use addObservations(int[][] states)."
			);
		}
		
		if (currentState == State.COMPUTING ||
		currentState == State.SETTING_PROPERTIES) {
			initialise();
		}
		currentState = State.ADDING_OBSERVATIONS;

		int rows = states.length;
		// increment the count of observations:
		observations += rows; 
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {

			if (!hashedStateCount.isEmpty()) {
				Object key = states[r];
				Integer value = hashedStateCount.getOrDefault(key, 0) + 1;
				hashedStateCount.put(key, value);
				continue;
			}

			Integer index = states[r];
			// check to make sure user didn't break the alphabet size they provided.
			if (knownIntegerRange && index >= this.alphabetSize) {
				throw new RuntimeException("Observation " + states[r] + 
				" exceeds provided alphabet size " + alphabetSize);
			}

			// if array is still a managable size...
			if (index < alphabetSize) {
				stateCount[index]++;
			} else {
				// Otherwise array has grown too large, move to hash table.
				hashedStateCount = new Hashtable<>();
				for (int i = 0; i < stateCount.length; i++) {
					if (stateCount[i] != 0) {
						hashedStateCount.put(i, stateCount[i]);
					}
				}
				hashedStateCount.put(index, 1);
				stateCount = null;
			}
		}
	}

	/**
	 * This method is intended for multi-dimensional use, and adds the observations
	 * into our estimates of the pdfs.
	 * 
	 * @param states - the observations of all dimensions at each time step
	 * 	states[i][j] represents the ith observed state of the jth dimension.
	 *  states[i] is the array of observations at time step i.
	 */
	@Override
	public void addObservations(int[][] states) throws RuntimeException {

		// ensure correct number of dimensions are provided
		for (int i = 0; i < states.length; i++) {
			if (states[i].length != numDimensions) {
				throw new RuntimeException(String.format(
					"Incorrect number of dimensions were given. Expected %d got %d.", numDimensions, states[i].length));
			}
		}

		if (currentState == State.COMPUTING ||
			currentState == State.SETTING_PROPERTIES) {
				initialise();
		}
		currentState = State.ADDING_OBSERVATIONS;

		
		// Unwrap unecessary single dimensional arrays, allows additional calls
		// to the single dimensional version without issue.
		if (numDimensions == 1) {
			for (int i = 0; i < states.length; i++) {
				this.stateCount[i] = states[i][0];
			}
			return;
		}
		
		int[] result;
		if (alphabetSizes != null) {
			result = computeCombinedValues(states, alphabetSizes, alphabetSize);
		} else {
			result = computeCombinedValues(states, alphabetSize);
		}
		
		for (int i = 0; i < result.length; i++) {
			this.stateCount[i] += result[i];
		}
	}

	/**
	 * This method is intended for multi-dimensional use, and adds the observations
	 * into our estimates of the pdfs.
	 * 
	 * @param states - the observations of all dimensions at each time step
	 */
	public void addObservations(Object[][] states) {

		// ensure correct number of dimensions are provided
		for (int i = 0; i < states.length; i++) {
			if (states[i].length != numDimensions) {
				throw new RuntimeException(String.format(
					"Incorrect number of dimensions were given. Expected %d got %d.", numDimensions, states[i].length));
			}
		}

		if (currentState == State.COMPUTING ||
			currentState == State.SETTING_PROPERTIES) {
				initialise();
		}
		currentState = State.ADDING_OBSERVATIONS;

		// Unwrap unecessary single dimensional arrays, allows additional calls
		// to the single dimensional version without issue.
		if (numDimensions == 1) {
			Object[] temp = new Object[states.length];
			for (int i = 0; i < states.length; i++) {
				temp[i] = states[i][0];
			}
			addObservations(temp);
			return;
		}

		for (int i = 0; i < states.length; i++) {
			List<Object> key = Arrays.asList(states[i]);
			Integer value = hashedStateCount.getOrDefault(key, 0) + 1;
			hashedStateCount.put(key, value);
		}
	}

	protected int[] computeCombinedValues(int[][] separatedValues, int maxSize) throws RuntimeException {
		int[] sizes = new int[numDimensions];
		MatrixUtils.fill(sizes, alphabetSize);
		return computeCombinedValues(separatedValues, sizes, maxSize);
	}

	// TODO -- move to MatrixUtils.
	protected int[] computeCombinedValues(int[][] separateValues, int[] alphabetSizes, int maxSize) throws RuntimeException {
		// Make sure we won't get any overflow here
		if (MatrixUtils.combinedValuesOverflow(numDimensions, alphabetSize)) {
			// multiplier has overflown
			throw new RuntimeException("Too many numDimensions " + numDimensions + " for the given alphabetSize " + alphabetSize +
					" for this call to computeCombinedValues");
		}
		
		int size = 1;
		for (int s : alphabetSizes) {
			size *= s;
		}

		if (size >= maxSize) {
			throw new RuntimeException("Size of combined observations exceeds capacity.");
		}

		int[] combinedValues = new int[size];
		for (int r = 0; r < separateValues.length; r++) {
			if (separateValues[r].length != alphabetSizes.length) {
				throw new RuntimeException("alphabetSize and given input size do not match");
			}
			int combinedRowValue = 0;
			int multiplier = 1;
			for (int c = numDimensions - 1; c >= 0; c--) {
				// check alphabet size restriction is maintained.
				if (separateValues[r][c] >= alphabetSizes[r]) {
					throw new RuntimeException("input value " + separateValues[r][c]
					 + " exceeded alphabet size" + alphabetSizes[c]);
				}
				combinedRowValue += separateValues[r][c] * multiplier;
				multiplier *= alphabetSizes[r];
			}
			combinedValues[r] = combinedRowValue;
		}
		return combinedValues;
	}
	
	@Override
	public void setObservations(Object[] observations) throws Exception {
		this.observations = observations.length;
		this.addObservations(observations);
		this.finaliseAddObservations();
	}
	
	@Override
	public double computeAverageLocalOfObservations() {
		double ent = 0.0;
		double p_state;
		max = 0;
		min = 0;
		currentState = State.COMPUTING;
		
		if (hashedStateCount.isEmpty()) {
			for (int count: stateCount) {
				// compute p_state
				p_state = (double) count / (double) observations;
				ent += computeEntropyCont(p_state);
			}
		} else {
			for (Map.Entry<Object, Integer> entry : hashedStateCount.entrySet()) {
				p_state = (double) entry.getValue() / (double) observations;
				ent += computeEntropyCont(p_state);
			}
		}
		
		average = ent;
		return ent;
	}
	
	@Override
	public final double computeAverageLocal(int states[]) throws Exception {
		if (states.length == 0) {
			throw new RuntimeException("States cannot be empty.");
		}
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	protected double computeEntropyCont(double p_state) {
		double entCont = 0.0;
		if (p_state > 0.0) {
			// Entropy takes the negative log:
			double localValue = - Math.log(p_state) / log_2;
			entCont = p_state * localValue;
			if (localValue > max) {
				max = localValue;
			} else if (localValue < min) {
				min = localValue;
			}
		}
		return entCont;
	}
	
	@Override
	public double[] computeLocalFromPreviousObservations(int states[]) {

		this.currentState = State.COMPUTING;
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

	/*
	 ***********************************************************************
	 ******************* ALL METHODS BELOW ARE DEPRECATED ******************
	 *******************        DON'T USE THEM            ******************
	 ***********************************************************************
	*/

	// TODO protected duplicates, keep public call, that does the same as private.
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

		if (currentState == State.SETTING_PROPERTIES) {
			throw new RuntimeException("Estimator should be initialised before finalised...");
		}

		if (observations == 0) {
			throw new RuntimeException("Must have some observations to finalise.");
		}

		this.currentState = State.COMPUTING;
	}
}