package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * Implements transfer entropy (see Schreiber, PRL, 2000)
 *  and local transfer entropy (see Lizier et al, PRE, 2008)
 * 
 * Usage:
 * 0. Construct or call newInstance
 * 1. Continuous accumulation of observations before computing :
 *   Call: a. initialise()
 *         b. addObservations() several times over
 *         c. computeLocalFromPreviousObservations() or computeAverageLocalOfObservations()
 * 2. Standalone computation from a single set of observations:
 *   Call: computeLocal() or computeAverageLocal()
 * 
 * @see For transfer entropy: Schreiber, PRL 85 (2) pp.461-464, 2000; http://dx.doi.org/10.1103/PhysRevLett.85.461
 * @see For local transfer entropy: Lizier et al, PRE 77, 026110, 2008; http://dx.doi.org/10.1103/PhysRevE.77.026110
 * 
 * @author Joseph Lizier
 * joseph.lizier at gmail.com
 * http://lizier.me/joseph/
 *
 */
public class ApparentTransferEntropyCalculator extends ContextOfPastMeasureCalculator 
	implements ChannelCalculator {

	protected int[][][] sourceNextPastCount = null;	// count for (source[n],dest[n+1],dest[n]^k) tuples
	protected int[][] sourcePastCount = null;			// count for (source[n],dest[n]^k) tuples
	protected boolean periodicBoundaryConditions = true;

	/**
	 * First time step at which we can take an observation
	 *  (needs to account for k previous steps)
	 */
	protected int startObservationTime = 1;
	
	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base
	 * @param history
	 * 
	 * @return
	 */
	public static ApparentTransferEntropyCalculator newInstance(int base, int history) {
		
		return new ApparentTransferEntropyCalculator(base, history);

		// Old code for an attempted optimisation:
		/*
		if (isPowerOf2(base)) {
			return new ApparentTransferEntropyCalculatorBase2(base, history);
		} else {
			return new ApparentTransferEntropyCalculator(base, history);
		}
		*/
	}
	
	/**
	 * Create a new TE calculator for the given base and history length.
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 */
	public ApparentTransferEntropyCalculator(int base, int history) {

		super(base, history);
		// Create storage for extra counts of observations
		sourceNextPastCount = new int[base][base][base_power_k];
		sourcePastCount = new int[base][base_power_k];
		
		// Which time step do we start taking observations from?
		// Normally this is k (to allow k previous time steps)
		//  but if k==0 (becoming a lagged MI), it's 1.
		startObservationTime = Math.max(k, 1);
	}

	/**
	 * Initialise calculator, preparing to take observation sets in
	 * Should be called prior to any of the addObservations() methods.
	 * You can reinitialise without needing to create a new object.
	 *
	 */
	public void initialise(){
		super.initialise();
		
		MatrixUtils.fill(sourceNextPastCount, 0);
		MatrixUtils.fill(sourcePastCount, 0);
	}
	
	/**
 	 * Add observations for a single source-destination pair of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call is for time series not part of the same 2D array
	 *
	 * @param states
	 */
	public void addObservations(int[] dest, int[] source) {
		int rows = dest.length;
		// increment the count of observations:
		observations += (rows - startObservationTime); 
		
		// Initialise and store the current previous value for each column
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += dest[p];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal;
		for (int r = startObservationTime; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = dest[r];
			sourceVal = source[r-1];
			sourceNextPastCount[sourceVal][destVal][pastVal]++;
			sourcePastCount[sourceVal][pastVal]++;
			nextPastCount[destVal][pastVal]++;
			pastCount[pastVal]++;
			nextCount[destVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[dest[r-k]];
				pastVal *= base;
				pastVal += dest[r];
			}
		}
	}

	/**
 	 * Add observations for a single source-destination pair of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call is for time series not part of the same 2D array.
	 * Start and end time are the (inclusive) indices within which to add the observations.
	 * The start time is from the earliest of the k historical values of the destination (inclusive),
	 *  the end time is the last destination time point to add in.
	 *
	 * @param dest
	 * @param source
	 * @param startTime
	 * @param endTime
	 * 
	 */
	public void addObservations(int[] dest, int[] source, int startTime, int endTime) {
		// increment the count of observations:
		observations += (endTime - startTime) - startObservationTime + 1; 
		
		// Initialise and store the current previous value for each column
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += dest[startTime + p];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal;
		for (int r = startTime + startObservationTime; r <= endTime; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = dest[r];
			sourceVal = source[r-1];
			sourceNextPastCount[sourceVal][destVal][pastVal]++;
			sourcePastCount[sourceVal][pastVal]++;
			nextPastCount[destVal][pastVal]++;
			pastCount[pastVal]++;
			nextCount[destVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[dest[r-k]];
				pastVal *= base;
				pastVal += dest[r];
			}
		}
	}

	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states 1st index is time, 2nd index is agent number
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. src i-j, dest i: transfer is j cells to the right) 
	 */
	public void addObservations(int states[][], int j) {
		int timeSteps = states.length;
		int agents = states[0].length;
		// increment the count of observations:
		if (periodicBoundaryConditions) {
			observations += (timeSteps - startObservationTime)*agents; 			
		} else {
			observations += (timeSteps - startObservationTime)*(agents - Math.abs(j));
		}
		
		// Initialise and store the current previous value for each column
		int[] pastVal = new int[agents]; 
		for (int c = 0; c < agents; c++) {
			pastVal[c] = 0;
			for (int p = 0; p < k; p++) {
				pastVal[c] *= base;
				pastVal[c] += states[p][c];
			}
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal;
		for (int r = startObservationTime; r < timeSteps; r++) {
			for (int c = 0; c < agents; c++) {
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				int sourceAgent = c-j;
				if ((sourceAgent < 0) || (sourceAgent >= agents)) {
					// Source agent is out of bounds unless we are using periodic boundary conditions
					if (periodicBoundaryConditions) {
						sourceAgent = (sourceAgent+agents) % agents;
					} else {
						// Don't add this to our observations
						continue;
					}
				}
				sourceVal = states[r-1][sourceAgent];
				destVal = states[r][c];
				sourceNextPastCount[sourceVal][destVal][pastVal[c]]++;
				sourcePastCount[sourceVal][pastVal[c]]++;
				nextPastCount[destVal][pastVal[c]]++;
				pastCount[pastVal[c]]++;
				nextCount[destVal]++;
				// Update the previous value:
				if (k > 0) {
					pastVal[c] -= maxShiftedValue[states[r-k][c]];
					pastVal[c] *= base;
					pastVal[c] += states[r][c];
				}
			}
		}		
	}

	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param h - number of rows to compute transfer entropy across
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. src (g-h,i-j), dest (g,i): transfer is h cells down, j cells to the right) 
	 */
	public void addObservations(int states[][][], int h, int j) {
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
		if (periodicBoundaryConditions) {
			observations += (timeSteps - startObservationTime) * agentRows * agentColumns;
		} else {
			observations += (timeSteps - startObservationTime) * (agentRows - Math.abs(h)) * (agentColumns - Math.abs(j));
		}
		
		// Initialise and store the current previous value for each agent
		int[][] pastVal = new int[agentRows][agentColumns]; 
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				pastVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					pastVal[r][c] *= base;
					pastVal[r][c] += states[p][r][c];
				}
			}
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal;
		for (int t = startObservationTime; t < timeSteps; t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					// Add to the count for this particular transition:
					// (cell's assigned as above)
					int sourceAgentRow = r-h;
					if ((sourceAgentRow < 0) || (sourceAgentRow >= agentRows)) {
						// Source agent is out of bounds unless we are using periodic boundary conditions
						if (periodicBoundaryConditions) {
							sourceAgentRow = (sourceAgentRow+agentRows) % agentRows;
						} else {
							// Don't add this to our observations
							continue;
						}
					}
					int sourceAgentColumn = c-j;
					if ((sourceAgentColumn < 0) || (sourceAgentColumn >= agentColumns)) {
						// Source agent is out of bounds unless we are using periodic boundary conditions
						if (periodicBoundaryConditions) {
							sourceAgentColumn = (sourceAgentColumn+agentColumns) % agentColumns;
						} else {
							// Don't add this to our observations
							continue;
						}
					}
					sourceVal = states[t-1][sourceAgentRow][sourceAgentColumn];
					destVal = states[t][r][c];
					sourceNextPastCount[sourceVal][destVal][pastVal[r][c]]++;
					sourcePastCount[sourceVal][pastVal[r][c]]++;
					nextPastCount[destVal][pastVal[r][c]]++;
					pastCount[pastVal[r][c]]++;
					nextCount[destVal]++;
					// Update the previous value:
					if (k > 0) {
						pastVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
						pastVal[r][c] *= base;
						pastVal[r][c] += states[t][r][c];
					}
				}
			}
		}		
	}

	/**
 	 * Add observations for a single source-destination pair of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd index is agent number
	 * @param destIndex destination agent index
	 * @param sourceIndex source agent index
	 */
	public void addObservations(int states[][], int destIndex, int sourceIndex) {
		int rows = states.length;
		// increment the count of observations:
		observations += (rows - startObservationTime); 
		
		// Initialise and store the current previous value for each column
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[p][destIndex];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal;
		for (int r = startObservationTime; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = states[r][destIndex];
			sourceVal = states[r-1][sourceIndex];
			sourceNextPastCount[sourceVal][destVal][pastVal]++;
			sourcePastCount[sourceVal][pastVal]++;
			nextPastCount[destVal][pastVal]++;
			pastCount[pastVal]++;
			nextCount[destVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[states[r-k][destIndex]];
				pastVal *= base;
				pastVal += states[r][destIndex];
			}
		}
	}

	/**
 	 * Add observations for a single source-destination pair of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param destRowIndex destination agent row index
	 * @param destColumnIndex destination agent column index
	 * @param sourceRowIndex source agent row index
	 * @param sourceColumnIndex source agent column index
	 */
	public void addObservations(int states[][][], int destRowIndex, int destColumnIndex,
												  int sourceRowIndex, int sourceColumnIndex) {
		int timeSteps = states.length;
		// increment the count of observations:
		observations += (timeSteps - startObservationTime); 
		
		// Initialise and store the current previous value for each column
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[p][destRowIndex][destColumnIndex];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal;
		for (int r = startObservationTime; r < timeSteps; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = states[r][destRowIndex][destColumnIndex];
			sourceVal = states[r-1][sourceRowIndex][sourceColumnIndex];
			sourceNextPastCount[sourceVal][destVal][pastVal]++;
			sourcePastCount[sourceVal][pastVal]++;
			nextPastCount[destVal][pastVal]++;
			pastCount[pastVal]++;
			nextCount[destVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[states[r-k][destRowIndex][destColumnIndex]];
				pastVal *= base;
				pastVal += states[r][destRowIndex][destColumnIndex];
			}
		}
	}

	/**
	 * 
	 * Returns the count of observations of the past given state dest[n]^k.
	 * The past state is indicated by a discrete integer representing the joint variable
	 *  of the k past states: (dest[n-k+1],dest[n-k+2],...,dest[n-1],dest[n]).
	 * The integer is computed as:<br/>
	 * pastVal = dest[n-k+1] * base^(k-1) + dest[n-k+2] * base^(k-2) + ... + dest[n-1] * base + dest[n]
	 * 
	 * 
	 * @param pastVal joint state of the past of the destination dest[n]^k
	 * @return count of observations of the given past state
	 */
	public int getPastCount(int pastVal) {
		return pastCount[pastVal];
	}
	
	/**
	 * 
	 * @param pastVal joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state
	 */
	public double getPastProbability(int pastVal) {
		return (double) pastCount[pastVal] / (double) observations;
	}

	/**
	 * 
	 * Returns the count of observations of the past given state dest[n]^k and next state dest[n+1].
	 * @see getPastCount(int) for how the joint value representing the past is calculated.
	 * 
	 * @param destVal next state of the destination dest[n+1]
	 * @param pastVal joint state of the past of the destination dest[n]^k
	 * @return count of observations of the given past state and next state
	 */
	public int getNextPastCount(int destVal, int pastVal) {
		return nextPastCount[destVal][pastVal];
	}
	
	/**
	 * 
	 * @param destVal next state of the destination dest[n+1]
	 * @param pastVal joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state and next state
	 */
	public double getNextPastProbability(int destVal, int pastVal) {
		return (double) nextPastCount[destVal][pastVal] / (double) observations;
	}

	/**
	 * 
	 * Returns the count of observations of the past given state dest[n]^k and the source state source[n].
	 * @see getPastCount(int) for how the joint value representing the past is calculated.
	 * 
	 * @param sourceVal state of the source source[n]
	 * @param pastVal joint state of the past of the destination dest[n]^k
	 * @return count of observations of the given past state and the source state
	 */
	public int getSourcePastCount(int sourceVal, int pastVal) {
		return sourcePastCount[sourceVal][pastVal];
	}
	
	/**
	 * 
	 * @param sourceVal state of the source source[n]
	 * @param pastVal joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state and the source state
	 */
	public double getSourcePastProbability(int sourceVal, int pastVal) {
		return (double) sourcePastCount[sourceVal][pastVal] / (double) observations;
	}

	/**
	 * 
	 * Returns the count of observations of the past given state dest[n]^k,
	 *  the next state of the destination dest[n+1] and the source state source[n].
	 * @see getPastCount(int) for how the joint value representing the past is calculated.
	 * 
	 * @param sourceVal state of the source source[n]
	 * @param nextVal next state of the destination dest[n+1]
	 * @param pastVal joint state of the past of the destination dest[n]^k
	 * @return count of observations of the given past state, next state of destination and the source state
	 */
	public int getSourceNextPastCount(int sourceVal, int destVal, int pastVal) {
		return sourceNextPastCount[sourceVal][destVal][pastVal];
	}
	
	/**
	 * 
	 * @param sourceVal state of the source source[n]
	 * @param nextVal next state of the destination dest[n+1]
	 * @param pastVal joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state, next state of destination and the source state
	 */
	public double getSourceNextPastProbability(int sourceVal, int destVal, int pastVal) {
		return (double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) observations;
	}

	/**
	 * 
	 * Returns the count of observations of the next state dest[n+1].
	 * 
	 * @param nextVal next state of the destination dest[n+1]
	 * @return count of observations of the given next state
	 */
	public int getNextCount(int destVal) {
		return nextCount[destVal];
	}
	
	/**
	 * 
	 * @param nextVal state of the next destination dest[n+1]
	 * @return probability of the given next state
	 */
	public double getNextProbability(int destVal) {
		return (double) nextCount[destVal] / (double) observations;
	}

	/**
	 * Returns the average local transfer entropy from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return
	 */
	public double computeAverageLocalOfObservations() {
		double te = 0.0;
		double teCont = 0.0;

		max = 0;
		min = 0;
		double meanSqLocals = 0;
		for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
			// compute p(past)
			// double p_past = (double) pastCount[pastVal] / (double) observations;
			for (int destVal = 0; destVal < base; destVal++) {
				// compute p(dest,past)
				// double p_dest_past = (double) destPastCount[destVal][pastVal] / (double) observations;
				for (int sourceVal = 0; sourceVal < base; sourceVal++) {
					// compute p(source,dest,past)
					double p_source_dest_past = (double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) observations;
					// compute p(source,past)
					// double p_source_past = (double) sourcePastCount[sourceVal][pastVal] / (double) observations;
					// Compute TE contribution:
					if (sourceNextPastCount[sourceVal][destVal][pastVal] != 0) {
						/* Double check: should never happen
						if ((sourcePastCount[sourceVal][pastVal] == 0) ||
							(destPastCount[destVal][pastVal] == 0) ||
							(pastCount[pastVal] == 0)) {
							throw new RuntimeException("one subcount was zero!!");
						}
						*/
						
						double logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) sourcePastCount[sourceVal][pastVal]) /
						 	((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
						double localValue = Math.log(logTerm) / log_2;
						teCont = p_source_dest_past * localValue;
						if (localValue > max) {
							max = localValue;
						} else if (localValue < min) {
							min = localValue;
						}
						// Add this contribution to the mean 
						//  of the squared local values
						meanSqLocals += teCont * localValue;
					} else {
						teCont = 0.0;
					}
					te += teCont;
				}
			}
		}
		
		average = te;
		std = Math.sqrt(meanSqLocals - average * average);
		return te;
	}
	
	/**
	 * Returns the average active information storage from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return
	 */
	public double computeAverageActiveInfoStorageOfObservations() {
		double active = 0.0;
		double activeCont = 0.0;

		for (int nextVal = 0; nextVal < base; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// Compute MI contribution:
				if (nextPastCount[nextVal][prevVal] != 0) {
					double logTerm = (double) nextPastCount[nextVal][prevVal] /
								(double) pastCount[prevVal] /
								p_next;
					double localValue = Math.log(logTerm) / log_base;
					activeCont = (nextPastCount[nextVal][prevVal] /
								(double) observations) * localValue;
				} else {
					activeCont = 0.0;
				}
				active += activeCont;
			}
		}
		
		return active;
	}

	/**
	 * Dump a debug print of the PDFs of our observations
	 */
	public void debugPrintObservations() {
		
		System.out.println("Src\tDst\tPast\tc(s,d,p)\tc(s,p)\tc(d,p)\tc(p)");
		for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
			for (int destVal = 0; destVal < base; destVal++) {
				for (int sourceVal = 0; sourceVal < base; sourceVal++) {
					// Compute TE contribution:
					System.out.println(sourceVal + "\t" + destVal + "\t" + pastVal + "\t" +
							sourceNextPastCount[sourceVal][destVal][pastVal] + "\t\t" +
							sourcePastCount[sourceVal][pastVal] + "\t" + 
							nextPastCount[destVal][pastVal] + "\t" +
							pastCount[pastVal]);
				}
			}
		}
	}

	/**
	 * Compute the significance of obtaining the given average TE from the given observations
	 * 
	 * 	This is as per Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128.
	 *
	 * @param numPermutationsToCheck number of new orderings of the source values to compare against
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) {
		double actualTE = computeAverageLocalOfObservations();
		
		// Reconstruct the source values (not necessarily in order)
		int[] sourceValues = new int[observations];
		int t_s = 0;
		for (int sourceVal = 0; sourceVal < base; sourceVal++) {
			// Count up the number of times this source value was observed:
			int numberOfSamples = 0;
			for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
				numberOfSamples += sourcePastCount[sourceVal][pastVal];
			}
			// Now add all of these as unordered observations:
			MatrixUtils.fill(sourceValues, sourceVal, t_s, numberOfSamples);
			t_s += numberOfSamples;
		}
		
		// And construct unordered (dest,past) tuples.
		// It doesn't matter that we've appeared to destroy the ordering here because
		//  the joint distribution nextPastCount is actually preserved in
		//  our construction of pastVal and destValues together here.
		int[] destValues = new int[observations];
		int[] pastValues = new int[observations];
		int t_d = 0;
		int t_p = 0;
		for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
			MatrixUtils.fill(pastValues, pastVal, t_p, pastCount[pastVal]);
			t_p += pastCount[pastVal];
			for (int destVal = 0; destVal < base; destVal++) {
				MatrixUtils.fill(destValues, destVal, t_d, nextPastCount[destVal][pastVal]);
				t_d += nextPastCount[destVal][pastVal];
			}
		}
		
		// Construct new source orderings based on the source probabilities only
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(observations, numPermutationsToCheck);

		ApparentTransferEntropyCalculator ate2 = newInstance(base, k);
		ate2.initialise();
		ate2.observations = observations;
		ate2.pastCount = pastCount;
		ate2.nextPastCount = nextPastCount;
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the source
			int[] newSourceData = MatrixUtils.extractSelectedTimePoints(sourceValues, newOrderings[p]);
			// compute the joint probability distributions
			MatrixUtils.fill(ate2.sourceNextPastCount, 0);
			MatrixUtils.fill(ate2.sourcePastCount, 0);
			for (int t = 0; t < observations; t++) {
				ate2.sourcePastCount[newSourceData[t]][pastValues[t]]++;
				ate2.sourceNextPastCount[newSourceData[t]][destValues[t]][pastValues[t]]++;
			}
			// And get a TE value for this realisation:
			double newTe = ate2.computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newTe;
			if (newTe >= actualTE) {
				countWhereTeIsMoreSignificantThanOriginal++;
			}

		}
		
		// And return the significance
		measDistribution.pValue = (double) countWhereTeIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualTE;
		return measDistribution;
	}
	
	/**
	 * Computes local transfer entropy for the given values
	 *  
	 * @param destNext
	 * @param destPast
	 * @param sourceCurrent
	 * @return
	 */
	public double computeLocalFromPreviousObservations(int destNext, int destPast, int sourceCurrent){

		double logTerm = ((double) sourceNextPastCount[sourceCurrent][destNext][destPast] / (double) sourcePastCount[sourceCurrent][destPast]) /
	 		((double) nextPastCount[destNext][destPast] / (double) pastCount[destPast]);
		return Math.log(logTerm) / log_2;
	}

	/**
	 * Computes local apparent transfer entropy for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only
	 *  
	 * @param states 1st index is time, 2nd index is agent number
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int destStates[], int sourceStates[]){
		int timeSteps = destStates.length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localTE = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int pastVal = 0; 
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += destStates[p];
		}
		int destVal, sourceVal;
		double logTerm;
		for (int t = startObservationTime; t < timeSteps; t++) {
			sourceVal = sourceStates[t-1];
			destVal = destStates[t];
			// Now compute the local value
			logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) sourcePastCount[sourceVal][pastVal]) /
		 		((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
			localTE[t] = Math.log(logTerm) / log_2;
			average += localTE[t];
			if (localTE[t] > max) {
				max = localTE[t];
			} else if (localTE[t] < min) {
				min = localTE[t];
			}
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[destStates[t-k]];
				pastVal *= base;
				pastVal += destStates[t];
			}
		}		

		if (periodicBoundaryConditions) {
			average = average/(double) (timeSteps - startObservationTime);
		} else {
			average = average/(double) (timeSteps - startObservationTime);
		}
		
		return localTE;
	}

	/**
	 * Computes local transfer for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only
	 *  
	 * @param states 1st index is time, 2nd index is agent number
	 * @return
	 */
	public double[][] computeLocalFromPreviousObservations(int states[][], int j){
		int timeSteps = states.length;
		int agents = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localTE = new double[timeSteps][agents];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int[] pastVal = new int[agents]; 
		for (int c = 0; c < agents; c++) {
			pastVal[c] = 0;
			for (int p = 0; p < k; p++) {
				pastVal[c] *= base;
				pastVal[c] += states[p][c];
			}
		}
		int destVal, sourceVal;
		double logTerm;
		for (int t = startObservationTime; t < timeSteps; t++) {
			for (int c = 0; c < agents; c++) {
				int sourceAgentIndex = c-j;
				if ((sourceAgentIndex < 0) || (sourceAgentIndex >= agents)) {
					// Source agent is out of bounds unless we are using periodic boundary conditions
					if (periodicBoundaryConditions) {
						sourceAgentIndex = (sourceAgentIndex+agents) % agents;
					} else {
						// Don't compute a local value for this one
						continue;
					}
				}
				sourceVal = states[t-1][sourceAgentIndex];
				destVal = states[t][c];
				// Now compute the local value
				logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal[c]] / (double) sourcePastCount[sourceVal][pastVal[c]]) /
			 		((double) nextPastCount[destVal][pastVal[c]] / (double) pastCount[pastVal[c]]);
				localTE[t][c] = Math.log(logTerm) / log_2;
				average += localTE[t][c];
				if (localTE[t][c] > max) {
					max = localTE[t][c];
				} else if (localTE[t][c] < min) {
					min = localTE[t][c];
				}
				// Update the previous value:
				if (k > 0) {
					pastVal[c] -= maxShiftedValue[states[t-k][c]];
					pastVal[c] *= base;
					pastVal[c] += states[t][c];
				}
			}
		}		

		if (periodicBoundaryConditions) {
			average = average/(double) ((timeSteps - startObservationTime) * agents);
		} else {
			average = average/(double) ((timeSteps - startObservationTime) * (agents - Math.abs(j)));
		}
		
		return localTE;
	}
	
	/**
	 * Computes local transfer for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only
	 *  
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param h - number of rows to compute transfer entropy across
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. src (g-h,i-j), dest (g,i): transfer is h cells down, j cells to the right) 
	 * @return
	 */
	public double[][][] computeLocalFromPreviousObservations(int states[][][], int h, int j){
		int timeSteps = states.length;
		int agentRows = states[0].length;
		int agentColumns = states[0][0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][][] localTE = new double[timeSteps][agentRows][agentColumns];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int[][] pastVal = new int[agentRows][agentColumns]; 
		for (int r = 0; r < agentRows; r++){
			for (int c = 0; c < agentColumns; c++) {
				pastVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					pastVal[r][c] *= base;
					pastVal[r][c] += states[p][r][c];
				}
			}
		}
		int destVal, sourceVal;
		double logTerm;
		for (int t = startObservationTime; t < timeSteps; t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					int sourceAgentRow = r-h;
					if ((sourceAgentRow < 0) || (sourceAgentRow >= agentRows)) {
						// Source agent is out of bounds unless we are using periodic boundary conditions
						if (periodicBoundaryConditions) {
							sourceAgentRow = (sourceAgentRow+agentRows) % agentRows;
						} else {
							// Don't compute a local value for this one
							continue;
						}
					}
					int sourceAgentColumn = c-j;
					if ((sourceAgentColumn < 0) || (sourceAgentColumn >= agentColumns)) {
						// Source agent is out of bounds unless we are using periodic boundary conditions
						if (periodicBoundaryConditions) {
							sourceAgentColumn = (sourceAgentColumn+agentColumns) % agentColumns;
						} else {
							// Don't compute a local value for this one
							continue;
						}
					}
					sourceVal = states[t-1][sourceAgentRow][sourceAgentColumn];
					destVal = states[t][r][c];
					// Now compute the local value
					logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal[r][c]] / (double) sourcePastCount[sourceVal][pastVal[r][c]]) /
				 		((double) nextPastCount[destVal][pastVal[r][c]] / (double) pastCount[pastVal[r][c]]);
					localTE[t][r][c] = Math.log(logTerm) / log_2;
					average += localTE[t][r][c];
					if (localTE[t][r][c] > max) {
						max = localTE[t][r][c];
					} else if (localTE[t][r][c] < min) {
						min = localTE[t][r][c];
					}
					// Update the previous value:
					if (k > 0) {
						pastVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
						pastVal[r][c] *= base;
						pastVal[r][c] += states[t][r][c];
					}
				}
			}
		}		

		if (periodicBoundaryConditions) {
			average = average/(double) ((timeSteps - startObservationTime) * agentRows * agentColumns);
		} else {
			average = average/(double) ((timeSteps - startObservationTime) * (agentRows - Math.abs(h)) *
												(agentColumns - Math.abs(j)));			
		}
		
		return localTE;
	}

	/**
	 * Computes local transfer for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states 1st index is time, 2nd index is agent number
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][], int destCol, int sourceCol){
		int rows = states.length;
		// int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localTE = new double[rows];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int pastVal = 0; 
		pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[p][destCol];
		}
		int destVal, sourceVal;
		double logTerm;
		for (int r = startObservationTime; r < rows; r++) {
			destVal = states[r][destCol];
			sourceVal = states[r-1][sourceCol];
			// Now compute the local value
			logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) sourcePastCount[sourceVal][pastVal]) /
		 		((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
			localTE[r] = Math.log(logTerm) / log_2;
			average += localTE[r];
			if (localTE[r] > max) {
				max = localTE[r];
			} else if (localTE[r] < min) {
				min = localTE[r];
			}
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[states[r-k][destCol]];
				pastVal *= base;
				pastVal += states[r][destCol];
			}
		}

		average = average/(double) (rows - startObservationTime);
		
		return localTE;
	}

	/**
	 * Computes local transfer for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param destRowIndex destination agent row index
	 * @param destColumnIndex destination agent column index
	 * @param sourceRowIndex source agent row index
	 * @param sourceColumnIndex source agent column index
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][][],
			int destRowIndex, int destColumnIndex, int sourceRowIndex, int sourceColumnIndex){
		int timeSteps = states.length;
		// int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localTE = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int pastVal = 0; 
		pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[p][destRowIndex][destColumnIndex];
		}
		int destVal, sourceVal;
		double logTerm;
		for (int r = startObservationTime; r < timeSteps; r++) {
			destVal = states[r][destRowIndex][destColumnIndex];
			sourceVal = states[r-1][sourceRowIndex][sourceColumnIndex];
			// Now compute the local value
			logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) sourcePastCount[sourceVal][pastVal]) /
		 		((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
			localTE[r] = Math.log(logTerm) / log_2;
			average += localTE[r];
			if (localTE[r] > max) {
				max = localTE[r];
			} else if (localTE[r] < min) {
				min = localTE[r];
			}
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[states[r-k][destRowIndex][destColumnIndex]];
				pastVal *= base;
				pastVal += states[r][destRowIndex][destColumnIndex];
			}
		}

		average = average/(double) (timeSteps - startObservationTime);
		
		return localTE;
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy between two time series
	 * Return a time series of local values.
	 * First history rows are zeros
	 * 
	 * @param destStates time series of destination states
	 * @param sourceStates time series of source states
	 * @return
	 */
	public double[] computeLocal(int destStates[], int sourceStates[]) {
		
		initialise();
		addObservations(destStates, sourceStates);
		return computeLocalFromPreviousObservations(destStates, sourceStates);
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method to be called for homogeneous agents only
	 * 
	 * @param states - 2D array of states
	 * @param j number of columns across which to compute the TE
	 * @return
	 */
	public double[][] computeLocal(int states[][], int j) {
		
		initialise();
		addObservations(states, j);
		return computeLocalFromPreviousObservations(states, j);
	}
	
	/**
	 * Standalone routine to 
	 * compute local transfer entropy across a 3D spatiotemporal
	 *  array of the states of homogeneous agents.
	 * Return a 3D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method to be called for homogeneous agents only
	 * 
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param h - number of rows to compute transfer entropy across
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. src (g-h,i-j), dest (g,i): transfer is h cells down, j cells to the right) 
	 * @return
	 */
	public double[][][] computeLocal(int states[][][], int h, int j) {
		
		initialise();
		addObservations(states, h, j);
		return computeLocalFromPreviousObservations(states, h, j);
	}

	/**
	 * Standalone routine to 
	 * compute average local transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return the average
	 * This method to be called for homogeneous agents only
	 * 
	 * @param states - 2D array of states
	 * @param j - TE across j cells to the right
	 * @return
	 */
	public double computeAverageLocal(int states[][], int j) {
		
		initialise();
		addObservations(states, j);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average local transfer entropy across a 3D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return the average
	 * This method to be called for homogeneous agents only
	 * 
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param h - number of rows to compute transfer entropy across
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. src (g-h,i-j), dest (g,i): transfer is h cells down, j cells to the right) 
	 * @return
	 */
	public double computeAverageLocal(int states[][][], int h, int j) {
		
		initialise();
		addObservations(states, h, j);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method suitable for heterogeneous agents
	 * 
	 * @param states - 2D array of states
	 * @param destCol - column index for the destination agent
	 * @param sourceCol - column index for the source agent
	 * @return
	 */
	public double[] computeLocal(int states[][], int destCol, int sourceCol) {
		
		initialise();
		addObservations(states, destCol, sourceCol);
		return computeLocalFromPreviousObservations(states, destCol, sourceCol);
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy across a 3D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method suitable for heterogeneous agents
	 * 
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param destRowIndex destination agent row index
	 * @param destColumnIndex destination agent column index
	 * @param sourceRowIndex source agent row index
	 * @param sourceColumnIndex source agent column index
	 * @return
	 */
	public double[] computeLocal(int states[][][], int destRowIndex, int destColumnIndex,
											int sourceRowIndex, int sourceColumnIndex) {
		
		initialise();
		addObservations(states, destRowIndex, destColumnIndex, sourceRowIndex, sourceColumnIndex);
		return computeLocalFromPreviousObservations(states, destRowIndex, destColumnIndex,
				sourceRowIndex, sourceColumnIndex);
	}

	/**
	 * Standalone routine to 
	 * compute average local transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Returns the average
	 * This method suitable for heterogeneous agents.
	 * 
	 * @param states - 2D array of states
	 * @param destCol - column index for the destination agent
	 * @param sourceCol - column index for the source agent
	 * @return
	 */
	public double computeAverageLocal(int states[][], int destCol, int sourceCol) {
		
		initialise();
		addObservations(states, destCol, sourceCol);
		return computeAverageLocalOfObservations();
	}
	
	/**
	 * Standalone routine to 
	 * compute average local transfer entropy across a 3D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Returns the average
	 * This method suitable for heterogeneous agents
	 * 
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param destRowIndex destination agent row index
	 * @param destColumnIndex destination agent column index
	 * @param sourceRowIndex source agent row index
	 * @param sourceColumnIndex source agent column index
	 * @return
	 */
	public double computeAverageLocal(int states[][][], int destRowIndex, int destColumnIndex,
			int sourceRowIndex, int sourceColumnIndex) {
		
		initialise();
		addObservations(states, destRowIndex, destColumnIndex, sourceRowIndex, sourceColumnIndex);
		return computeAverageLocalOfObservations();
	}

	public boolean isPeriodicBoundaryConditions() {
		return periodicBoundaryConditions;
	}
	public void setPeriodicBoundaryConditions(boolean periodicBoundaryConditions) {
		this.periodicBoundaryConditions = periodicBoundaryConditions;
	}
}
