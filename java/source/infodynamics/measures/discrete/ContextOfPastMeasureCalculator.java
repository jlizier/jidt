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
 * @author Joseph Lizier
 * 
 * Info theoretic measure calculator base class for 
 *  measures which require the context of the past
 *  history of the destination variable.
 * 
 * Usage:
 * 1. Continuous accumulation of observations before computing :
 *   Call: a. initialise()
 *         b. addObservations() several times over
 *         c. computeLocalFromPreviousObservations() or computeAverageLocalOfObservations()
 * 2. Standalone computation from a single set of observations:
 *   Call: computeLocal() or computeAverageLocal()
 * 
 * @author Joseph Lizier
 * joseph.lizier at gmail.com
 * http://lizier.me/joseph/
 *
 */
public abstract class ContextOfPastMeasureCalculator extends
		InfoMeasureCalculator {

	protected int k = 0; // history length k.
	protected boolean noObservationStorage = false;
	protected int[][] nextPastCount = null; // Count for (i[t+1], i[t]) tuples
	protected int[] pastCount = null; // Count for i[t]
	protected int[] nextCount = null; // count for i[t+1]
	protected int[] maxShiftedValue = null; // states * (base^(k-1))
	protected int base_power_k = 0;

	/**
	 * @param base
	 */
	public ContextOfPastMeasureCalculator(int base, int history) {
		this(base, history, false);
	}

	/**
	 * Constructor to be used by child classes only.
	 * In general, only needs to be explicitly called if child classes
	 *  do not wish to create the observation arrays.
	 * 
	 * @param base
	 * @param history
	 * @param dontCreateObsStorage
	 */
	protected ContextOfPastMeasureCalculator(int base, int history, boolean dontCreateObsStorage) {
		super(base);

		k = history;
		base_power_k = MathsUtils.power(base, k);
		
		// Relax the requirement that k >= 1, so that we can 
		//  eliminate considering the history at will ...
		//if (k < 1) {
		//	throw new RuntimeException("History k " + history + " is not >= 1 a ContextOfPastMeasureCalculator");
		//}
		
		// Check that we can convert the history value into an integer ok: 
		if (k > Math.log(Integer.MAX_VALUE) / log_base) {
			throw new RuntimeException("Base and history combination too large");
		}

		// Create constants for tracking prevValues
		maxShiftedValue = new int[base];
		for (int v = 0; v < base; v++) {
			maxShiftedValue[v] = v * MathsUtils.power(base, k-1);
		}
		
		noObservationStorage = dontCreateObsStorage;
		if (!dontCreateObsStorage) {
			// Create storage for counts of observations
			nextPastCount = new int[base][base_power_k];
			pastCount = new int[base_power_k];
			nextCount = new int[base];
		}
	}
	
	/**
	 * Initialise calculator, preparing to take observation sets in
	 * Should be called prior to any of the addObservations() methods.
	 * You can reinitialise without needing to create a new object.
	 *
	 */
	public void initialise() {
		super.initialise();
		
		if (!noObservationStorage) {
			MatrixUtils.fill(nextPastCount, 0);
			MatrixUtils.fill(pastCount, 0);
			MatrixUtils.fill(nextCount, 0);
		}
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
