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
 * A base class for calculators computing measures for
 * a single variable which
 * require knowledge of the embedded past state of a univariate
 * discrete (ie int[]) variable.
 * 
 * <p>This combines functionality for single agents from
 * {@link UnivariateMeasureDiscrete} with functionality
 * required in the context of the past provided by
 * {@link ContextOfPastMeasureCalculatorDiscrete}.</p>
 * 
 * <p>Usage is as defined in {@link InfoMeasureCalculatorDiscrete}, with
 * extra methods for supplying observations and making 
 * calculations defined in {@link UnivariateMeasureDiscrete}</p>.
 * 
 * <p>Users should not need to deal with this class directly;
 * it is simply used to gather common functionality for several
 * child classes.
 * </p>
 * 
 * TODO Make the Active info storage and entropy calculators inherit from this
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class UnivariateMeasureDiscreteInContextOfPastCalculator extends
		ContextOfPastMeasureCalculatorDiscrete implements UnivariateMeasureDiscrete {

	/**
	 * Construct the calculator with default base of 2 and history 1
	 */
	public UnivariateMeasureDiscreteInContextOfPastCalculator() {
		super(2, 1);
	}
	
	/**
	 * Construct the calculator
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history embedding length
	 */
	public UnivariateMeasureDiscreteInContextOfPastCalculator(int base, int history) {
		super(base, history);
	}
	
	/**
	 * Construct the calculator
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history embedding length
	 * @param dontCreateObsStorage do not create storage
	 * 		for observations of the embedded past (as the child
	 * 		class is signalling that it does not need it)
	 */
	public UnivariateMeasureDiscreteInContextOfPastCalculator(int base, int history, boolean dontCreateObsStorage) {
		super(base, history, dontCreateObsStorage);
	}

	
	public final double[] computeLocal(int[] states) {
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	@Override
	public final double computeAverageLocal(int[] states) {
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}
	
	public final double[][] computeLocal(int[][] states) {
		initialise();
		try {
			addObservations(states);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		return computeLocalFromPreviousObservations(states);
	}

	public final double[][][] computeLocal(int[][][] states) {
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	public final double computeAverageLocal(int[][] states) {
		initialise();
		try {
			addObservations(states);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		return computeAverageLocalOfObservations();
	}

	public final double computeAverageLocal(int[][][] states) {
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}
	
	public final double[] computeLocal(int[][] states, int col) {
		initialise();
		addObservations(states, col);
		return computeLocalFromPreviousObservations(states, col);
	}

	public final double[] computeLocal(int[][][] states, int index1, int index2) {
		initialise();
		addObservations(states, index1, index2);
		return computeLocalFromPreviousObservations(states, index1, index2);
	}

	public final double computeAverageLocal(int[][] states, int col) {
		initialise();
		addObservations(states, col);
		return computeAverageLocalOfObservations();
	}
	
	public final double computeAverageLocal(int[][][] states, int index1, int index2) {
		initialise();
		addObservations(states, index1, index2);
		return computeAverageLocalOfObservations();
	}

	// these are only here to get this thing to compile, they were in
	// the interface, and maybe I could or should have done this in the classes
	// that don't need them at all, but whatever, I think these should be
	// removed or only added to the classes that need them.

	private double[][] computeLocalFromPreviousObservations(int[][] states) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Unimplemented method 'computeLocalFromPreviousObservations'");
	}

	private void addObservations(int[][][] states) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Unimplemented method 'addObservations'");
	}

	private double[][][] computeLocalFromPreviousObservations(int[][][] states) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Unimplemented method 'computeLocalFromPreviousObservations'");
	}
	
	private double[] computeLocalFromPreviousObservations(int[][] states, int col) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Unimplemented method 'computeLocalFromPreviousObservations'");
	}
	
	private double[] computeLocalFromPreviousObservations(int[][][] states, int index1, int index2) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Unimplemented method 'computeLocalFromPreviousObservations'");
	}

	private void addObservations(int[][] states, int col) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Unimplemented method 'addObservations'");
	}

	private void addObservations(int[][][] states, int index1, int index2) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Unimplemented method 'addObservations'");
	} 
}
