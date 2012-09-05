package infodynamics.measures.discrete;

import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * Usage:
 * 1. Continuous accumulation of observations:
 *   Call: a. initialise()
 *         b. addObservations() several times over
 *         c. computeLocalFromPreviousObservations()
 * 2. Standalone:
 *   Call: localMultiInformation()
 * 
 * @author Joseph Lizier
 *
 */
public class MultiInformationCalculator extends InfoMeasureCalculator {
	
	private int[]	jointCount = null; // jointCount[jointState]
	private int[][] marginalCounts = null; // marginalCounts[marginalIndex][state]
	
	private int numVars;
	private int jointStates;
	private boolean checkedFirst = false;

	/**
	 * Constructor
	 *
	 */
	public MultiInformationCalculator(int base, int numVars) {
		super(base);
		this.numVars = numVars;
		jointStates = MathsUtils.power(base, numVars);
		jointCount = new int[jointStates];
		marginalCounts = new int[numVars][base];
	}

	/**
	 * Initialise calculator, preparing to take observation sets in
	 *
	 */
	public void initialise(){
		super.initialise();
		MatrixUtils.fill(jointCount, 0);
		MatrixUtils.fill(marginalCounts, 0);
	}
	
	/**
	 * Add observations of the variables based at every point, with the set defined by 
	 *  the group offsets array
	 *  
	 * @param states
	 * @param groupOffsets
	 */
	public void addObservations(int[] states, int[] groupOffsets) {
		for (int c = 0; c < states.length; c++) {
			// Add the marginal observations in, and compute the joint state value
			int jointValue = 0;
			for (int i = 0; i < numVars; i++) {
				int thisValue = states[(c + groupOffsets[i] + states.length) % states.length];
				marginalCounts[i][thisValue]++;
				jointValue *= base;
				jointValue += thisValue;
			}
			jointCount[jointValue]++;
			observations++;
		}
	}
	
	/**
	 * Add observations of the variables based at every point, with the set defined by 
	 *  the group offsets array
	 *  
	 * @param states
	 * @param groupOffsets
	 */
	public void addObservations(int[] states, int destinationIndex, int[] groupOffsets) {
		// Add the marginal observations in, and compute the joint state value
		int jointValue = 0;
		for (int i = 0; i < numVars; i++) {
			int thisValue = states[(destinationIndex + groupOffsets[i] + states.length) % states.length];
			marginalCounts[i][thisValue]++;
			jointValue *= base;
			jointValue += thisValue;
		}
		jointCount[jointValue]++;
		observations++;
	}

	/**
	 * Add observations of the variables based at every point, with the set defined by 
	 *  the group offsets array
	 *  
	 * @param states
	 * @param groupOffsets
	 */
	public void addObservations(int[][] states, int[] groupOffsets) {
		for (int t = 0; t < states.length; t++) {
			for (int c = 0; c < states[t].length; c++) {
				// Add the marginal observations in, and compute the joint state value
				int jointValue = 0;
				for (int i = 0; i < numVars; i++) {
					int thisValue = states[t][(c + groupOffsets[i] + states.length) % states.length];
					marginalCounts[i][thisValue]++;
					jointValue *= base;
					jointValue += thisValue;
				}
				jointCount[jointValue]++;
				observations++;
			}
		}
	}
	
	/**
	 * Returns the average local multi information storage from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return
	 */
	public double computeAverageLocalOfObservations() {
		
		int[] jointTuple = new int[numVars];
		checkedFirst = false;
		average = computeMiOfGivenTupleFromVarIndex(jointTuple, 0);
				
		return average;
	}
	
	/**
	 * Compute the contribution to the MI for all tuples starting with tuple[0..(fromIndex-1)].
	 * 
	 * @param tuple
	 * @param fromIndex
	 * @return
	 */
	public double computeMiOfGivenTupleFromVarIndex(int[] tuple, int fromIndex) {
		double miCont = 0;
		if (fromIndex == numVars) {
			// The whole tuple is filled in, so compute the contribution to the MI from this tuple
			double prodMarginalProbs = 1.0;
			int jointValue = 0;
			for (int i = 0; i < numVars; i++) {
				prodMarginalProbs *= (double) marginalCounts[i][tuple[i]] / (double) observations;
				jointValue *= base;
				jointValue += tuple[i];
			}
			if (jointCount[jointValue] == 0) {
				// This joint state does not occur, so it makes no contribution here
				return 0;
			}
			double jointProb = (double) jointCount[jointValue] / (double) observations;
			double logValue = jointProb / prodMarginalProbs;
			double localValue = Math.log(logValue) / log_2;
			if (jointProb > 0.0) {
				if (!checkedFirst) {
					max = localValue;
					min = localValue;
					checkedFirst = true;
				} else {
					if (localValue > max) {
						max = localValue;
					}
					if (localValue < min) {
						min = localValue;
					}
				}
			}
			miCont = jointProb * localValue;
		} else {
			// Fill out the next part of the tuple and make the recursive calls
			for (int v = 0; v < base; v++) {
				tuple[fromIndex] = v;
				miCont += computeMiOfGivenTupleFromVarIndex(tuple, fromIndex + 1);
			}
		}
		return miCont;
	}
}
