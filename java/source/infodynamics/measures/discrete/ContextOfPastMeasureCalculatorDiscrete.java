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
 * A base class for calculators computing measures which
 * require knowledge of the embedded past state of a univariate
 * discrete (ie int[]) variable.
 * 
 * Usage is as per {@link InfoMeasureCalculatorDiscrete}, but with some
 * extra utility functions provided for computing embedding vectors.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class ContextOfPastMeasureCalculatorDiscrete extends
		InfoMeasureCalculatorDiscrete {

	/**
	 * History length for the embedding
	 */
	protected int k = 0;
	/**
	 * Do not create storage
	 * 		for observations of the embedded past
	 */
	protected boolean noObservationStorage = false;
	/**
	 * Counts of (next,embedded_past) tuples
	 */
	protected int[][] nextPastCount = null; // Count for (i[t+1], i[t]) tuples
	/**
	 * Counts of (embedded_past) tuples
	 */
	protected int[] pastCount = null; // Count for i[t]
	/**
	 * Counts of (next) observations
	 */
	protected int[] nextCount = null; // count for i[t+1]
	/**
	 * Cached value maxShiftedValue[i] is i * (base^(k-1))
	 */
	protected int[] maxShiftedValue = null; // 
	/**
	 * Cached value of base^k
	 */
	protected int base_power_k = 0;

	/**
	 * Construct an instance with default base 2 and history 1
	 *  and not keeping sample counts
	 */
	public ContextOfPastMeasureCalculatorDiscrete() {
		this(2, 1, false);
	}

	/**
	 * Construct an instance
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history embedding length
	 */
	public ContextOfPastMeasureCalculatorDiscrete(int base, int history) {
		this(base, history, false);
	}

	/**
	 * Constructor to be used by child classes only.
	 * In general, only needs to be explicitly called if child classes
	 *  do not wish to create the observation arrays.
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history embedding length
	 * @param dontCreateObsStorage do not create storage
	 * 		for observations of the embedded past (as the child
	 * 		class is signalling that it does not need it)
	 */
	protected ContextOfPastMeasureCalculatorDiscrete(int base, int history, boolean dontCreateObsStorage) {
		super(base);
		resetHistory(history);
		
		noObservationStorage = dontCreateObsStorage;
	}

	/**
	 * Should be called after {@link #resetAlphabetSize(int)} has just been called.
	 * 
	 * @param history
	 */
	private void resetHistory(int history) {
		
		k = history;
		base_power_k = MathsUtils.power(alphabetSize, k);
		
		if (k < 0) {
			throw new RuntimeException("History k " + history + " is not >= 0 for a ContextOfPastMeasureCalculator");
		}
		
		// Check that we can convert the history value into an integer ok: 
		if (k > Math.log(Integer.MAX_VALUE) / log_base) {
			throw new RuntimeException("Base and history combination too large");
		}

		// Create constants for tracking prevValues
		maxShiftedValue = new int[alphabetSize];
		for (int v = 0; v < alphabetSize; v++) {
			maxShiftedValue[v] = v * MathsUtils.power(alphabetSize, k-1);
		}
	}
	
	@Override
	public void initialise() {
		initialise(alphabetSize, k);
	}
	
	/**
	 * Initialise the calculator with (potentially new) 
	 * specified base and history
	 * 
	 * @param base
	 * @param history
	 */
	public void initialise(int base, int history) {
		boolean baseOrHistoryChanged = false;
		if ((this.alphabetSize != base) || (k != history)) {
			baseOrHistoryChanged = true;
		}
		
		super.initialise(base);
		if (baseOrHistoryChanged) {
			resetHistory(history);
		}
		
		if (!noObservationStorage) {
			// We manage storage for relevant sample counts here:
			if (baseOrHistoryChanged || (nextPastCount == null)) {
				// Create new storage for counts of observations
				try {
					nextPastCount = new int[base][base_power_k];
					pastCount = new int[base_power_k];
					nextCount = new int[base];
				} catch (OutOfMemoryError e) {
					// Allow any Exceptions to be thrown, but catch and wrap
					//  Error as a RuntimeException
					throw new RuntimeException("Requested memory for the base " +
							base + " and k=" + k + " is too large for the JVM at this time", e);
				}
			} else {
				// Storage for sample counts exists and needs to be reset
				MatrixUtils.fill(nextPastCount, 0);
				MatrixUtils.fill(pastCount, 0);
				MatrixUtils.fill(nextCount, 0);
			}
		}
	}
	
	/**
	 * Utility function to compute the combined embedded
	 * past values of x up to and including time step t
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 * 
	 * @param x time-series data
	 * @param t compute embedding vector up to and
	 *  including index t
	 * @return int value representing the embedding vector
	 *  translated into a unique integer
	 */
	public int computePastValue(int[] x, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= alphabetSize;
			pastVal += x[t - k + 1 + p];
		}
		return pastVal;
	}

	/**
	 * Utility function to compute the combined embedded
	 * past values of x up to and including time step t
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 * where x is a column in data 
	 * 
	 * @param data multivariate time-series data
	 *  (first index is time, second is variable number)
	 * @param columnNumber which column to embed
	 * @param t compute embedding vector up to and
	 *  including index t
	 * @return int value representing the embedding vector
	 *  translated into a unique integer
	 */
	public int computePastValue(int[][] data, int columnNumber, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= alphabetSize;
			pastVal += data[t - k + 1 + p][columnNumber];
		}
		return pastVal;
	}

	/**
	 * Utility function to compute the combined embedded 
	 * past values of x up to and including time step t
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 * where x is a time-series for a given row and
	 * column in data 
	 * 
	 * @param data multivariate time-series data
	 *  (first index is time, second is row number for the variable
	 *   and third is column number for the variable)
	 * @param rowNumber row number of the variable to embed
	 * @param columnNumber column number of the variable to embed
	 * @param t compute embedding vector up to and
	 *  including index t
	 * @return int value representing the embedding vector
	 *  translated into a unique integer
	 */
	public int computePastValue(int[][][] data, int rowNumber, int columnNumber, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= alphabetSize;
			pastVal += data[t - k + 1 + p][rowNumber][columnNumber];
		}
		return pastVal;
	}
}
