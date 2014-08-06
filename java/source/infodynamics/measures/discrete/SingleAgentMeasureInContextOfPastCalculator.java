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

/**
 * Combines functionality for single agents with functionality
 * required in the context of the past.
 * 
 * @author Joseph Lizier
 *
 */
public abstract class SingleAgentMeasureInContextOfPastCalculator extends
		ContextOfPastMeasureCalculator implements SingleAgentMeasure {

	public SingleAgentMeasureInContextOfPastCalculator(int base, int history) {
		super(base, history);
	}

	public final double[][] computeLocal(int[][] states) {
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	public final double[][][] computeLocal(int[][][] states) {
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	public final double computeAverageLocal(int[][] states) {
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	public final double computeAverageLocal(int[][][] states) {
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	public final double[] computeLocalAtAgent(int[][] states, int col) {
		initialise();
		addObservations(states, col);
		return computeLocalFromPreviousObservations(states, col);
	}

	public final double[] computeLocalAtAgent(int[][][] states, int index1, int index2) {
		initialise();
		addObservations(states, index1, index2);
		return computeLocalFromPreviousObservations(states, index1, index2);
	}

	public final double computeAverageLocalAtAgent(int[][] states, int col) {
		initialise();
		addObservations(states, col);
		return computeAverageLocalOfObservations();
	}

}
