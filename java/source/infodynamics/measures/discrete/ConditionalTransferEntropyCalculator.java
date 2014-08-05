package infodynamics.measures.discrete;

import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;


/**
 * <p>Implements <i>conditional</i> transfer entropy,
 * and <i>local conditional</i> transfer entropy
 * (see Lizier et al., PRE, 2008, and Lizier et al., Chaos, 2010).
 * This class can be used for to compute <i>complete</i> transfer entropy 
 *  (see Lizier et al, PRE, 2008) which conditions on <b>all</b> 
 *  other causal information contributors to the destination.</p>
 *  
 * <p>Specifically, this implements the complete transfer entropy for 
 * <i>discrete</i>-valued variables.</p>
 *
 * <p>
 * The conditional sources (specified using either their
 *  offsets from the destination variable or their absolute column numbers
 *  in the multivariate data set)
 *  should be supplied in the same order in every method call, otherwise the answer supplied will
 *  be incorrect.
 * </p>
 * 
 * <p>Ideally, this class would extend ContextOfPastMeasure, however
 *  by conditioning on other info contributors, we need to alter
 *  the arrays pastCount and nextPastCount to consider all
 *  conditioned variables (i.e. other sources) also.
 * </p>
 * 
 * <p>Usage:
 * <ol>
 * 	<li>Construct: {@link #CompleteTransferEntropyCalculator(int, int)}</li>
 * 	<li>Initialise: {@link #initialise()}</li>
 *  <li>Either:
 *  	<ol>
 *  	<li>Continuous accumulation of observations then measurement; call:
 *  		<ol>
 *  			<li>{@link #addObservations(int[][], int, int[])} or related calls
 *         several times over - <b>note:</b> each method call adding 
 *         observations can be viewed as updating the PDFs; they do not
 *         append the separate time series (this would be incorrect behaviour
 *         for the transfer entropy, since the start of one time series
 *         is not necessarily related to the end of the other).</li>
 *   			<li>The compute relevant quantities, e.g.
 *   	   {@link #computeLocalFromPreviousObservations(int[][], int, int[])} or
 *         {@link #computeAverageLocalOfObservations()}</li>
 *       	</ol>
 * 		<li>or Standalone computation from a single set of observations;
 *   call e.g.: {@link #computeLocal(int[][], int, int[])} or
 *   {@link #computeAverageLocal(int[][], int, int, int[])}.>/li>
 *     </ol>
 * </ol>
 * </p>
 * 
 * @see "Schreiber, Physical Review Letters 85 (2) pp.461-464, 2000;
 *  <a href='http://dx.doi.org/10.1103/PhysRevLett.85.461'>download</a>
 *  (for definition of transfer entropy)"
 * @see "Lizier, Prokopenko and Zomaya, Physical Review E 77, 026110, 2008;
 * <a href='http://dx.doi.org/10.1103/PhysRevE.77.026110'>download</a>
 *  (for definition of <i>local</i> transfer entropy and 
 *  <i>complete</i> transfer entropy)"
 * @see "Lizier, Prokopenko and Zomaya, Chaos vol. 20, no. 3, 037109, 2010;
 * <a href='http://dx.doi.org/10.1063/1.3486801'>download</a>
 *  (for definition of <i>conditional</i> transfer entropy)"
 *  
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 * 
 * TODO Add methods for passing in single time series
 *
 */
public class ConditionalTransferEntropyCalculator extends InfoMeasureCalculator {

	protected int k = 0; // history length k.
	protected int base_power_k = 0;
	protected int base_power_num_others = 0;
	protected int numOtherInfoContributors = 0;
	protected int[][][][] sourceDestPastOthersCount = null;	// count for (i-j[n],i[n+1],i[n]^k,others) tuples
	protected int[][][] sourcePastOthersCount = null;			// count for (i-j[n],i[n]^k,others) tuples
	protected int[][][] destPastOthersCount = null; // Count for (i[n+1], i[n]^k,others) tuples
	protected int[][] pastOthersCount = null; // Count for (i[n]^k,others)
	protected int[] maxShiftedValue = null; // states * (base^(k-1))

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
	 * @param numOtherInfoContributors
	 * 
	 * @return
	 */
	public static ConditionalTransferEntropyCalculator
		newInstance(int base, int history, int numOtherInfoContributors) {
		
		return new ConditionalTransferEntropyCalculator
					(base, history, numOtherInfoContributors);
		
		// Old code for an attempted optimisation:
		/*
		if (isPowerOf2(base)) {
			return new CompleteTransferEntropyCalculatorBase2
						(base, history, numOtherInfoContributors);
		} else {
			return new CompleteTransferEntropyCalculator
					(base, history, numOtherInfoContributors);
		}
		*/
	}
	
	/**
	 * 
	 * 
	 * @param base
	 * @param history
	 * @param numOtherInfoContributors number of information contributors
	 *   (other than the past of the destination, if history < 1,
	 *   of the source) to condition on.
	 */
	public ConditionalTransferEntropyCalculator
		(int base, int history, int numOtherInfoContributors) {

		super(base);
		
		k = history;
		this.numOtherInfoContributors = numOtherInfoContributors;
		base_power_k = MathsUtils.power(base, k);
		base_power_num_others = MathsUtils.power(base, numOtherInfoContributors);
		
		// Relaxing this assumption so we can use this calculation as
		//  a time-lagged conditional MI at will:
		//if (k < 1) {
		//	throw new RuntimeException("History k " + history + " is not >= 1 a ContextOfPastMeasureCalculator");
		//}
		
		// Which time step do we start taking observations from?
		// Normally this is k (to allow k previous time steps)
		//  but if k==0 (becoming a lagged MI), it's 1.
		startObservationTime = Math.max(k, 1);

		// check that we can convert the base tuple into an integer ok
		if (k > Math.log(Integer.MAX_VALUE) / log_base) {
			throw new RuntimeException("Base and history combination too large");
		}
		if (numOtherInfoContributors < 1) {
			throw new RuntimeException("Number of other info contributors < 1 for CompleteTECalculator");
		}
		
		// Create storage for counts of observations
		sourceDestPastOthersCount = new int[base][base][base_power_k][base_power_num_others];
		sourcePastOthersCount = new int[base][base_power_k][base_power_num_others];
		destPastOthersCount = new int [base][base_power_k][base_power_num_others];
		pastOthersCount = new int[base_power_k][base_power_num_others];
		
		// Create constants for tracking prevValues
		maxShiftedValue = new int[base];
		for (int v = 0; v < base; v++) {
			maxShiftedValue[v] = v * MathsUtils.power(base, k-1);
		}
	}

	/**
	 * Initialise calculator, preparing to take observation sets in
	 * Should be called prior to any of the addObservations() methods.
	 * You can reinitialise without needing to create a new object.
	 *
	 */
	public void initialise(){
		super.initialise();
		
		MatrixUtils.fill(sourceDestPastOthersCount, 0);
		MatrixUtils.fill(sourcePastOthersCount, 0);
		MatrixUtils.fill(destPastOthersCount, 0);
		MatrixUtils.fill(pastOthersCount, 0);
	}
	
	/**
 	 * Add observations for a single source-destination-conditionals set
 	 *  to our estimates of the pdfs.
	 *
	 * @param dest destination time series
	 * @param source source time series
	 * @param conditionals conditionals multivariate time series
	 */
	public void addObservations(int[] dest, int[] source, int[][] conditionals)
		throws Exception {

		int rows = dest.length;
		
		if ((source.length != rows) || (conditionals.length != rows)) {
			throw new Exception("Number of observations must match for dest, source and conditionals");
		}
		
		if (rows - startObservationTime <= 0) {
			return;
		}
		
		// increment the count of observations:
		observations += (rows - startObservationTime); 
		
		if (numOtherInfoContributors != conditionals[0].length) {
			throw new Exception(String.format("conditionals does not have the expected number of variables (%d)", numOtherInfoContributors));
		}

		// Initialise and store the current previous value
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += dest[p];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, othersVal;
		for (int r = startObservationTime; r < rows; r++) {
			// Add to the count for this particular transition:
			destVal = dest[r];
			sourceVal = source[r-1];
			othersVal = 0;
			for (int o = 0; o < numOtherInfoContributors; o++) {
				// Include this other contributor
				othersVal *= base;
				othersVal += conditionals[r-1][o];
			}
			sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal]++;
			sourcePastOthersCount[sourceVal][pastVal][othersVal]++;
			destPastOthersCount[destVal][pastVal][othersVal]++;
			pastOthersCount[pastVal][othersVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[dest[r-k]];
				pastVal *= base;
				pastVal += dest[r];
			}
		}
	}

	/**
 	 * <p>Add observations for a single source-destination-conditionals set
 	 *  to our estimates of the pdfs.</p>
 	 *  
 	 * <p>This method takes in a univariate conditionals time series - it is assumed
 	 * that either numOtherInfoContributors == 1 or the user has combined the 
 	 * multivariate tuples into a single value for each observation, e.g. by calling
 	 * {@link MatrixUtils#computeCombinedValues(int[][], int)}. This cannot be
 	 * checked here however, so use at your own risk!
 	 * </p>
	 *
	 * @param dest destination time series
	 * @param source source time series
	 * @param conditionals conditionals univariate time series - it is assumed
	 *  that the user has combined the values of the multivariate conditionals time series
	 */
	public void addObservations(int[] dest, int[] source, int[] conditionals)
		throws Exception {

		int rows = dest.length;
		
		if ((source.length != rows) || (conditionals.length != rows)) {
			throw new Exception("Number of observations must match for dest, source and conditionals");
		}
		
		if (rows - startObservationTime <= 0) {
			return;
		}
		
		// increment the count of observations:
		observations += (rows - startObservationTime); 
		
		// Initialise and store the current previous value
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += dest[p];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, othersVal;
		for (int r = startObservationTime; r < rows; r++) {
			// Add to the count for this particular transition:
			destVal = dest[r];
			sourceVal = source[r-1];
			othersVal = conditionals[r-1];
			sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal]++;
			sourcePastOthersCount[sourceVal][pastVal][othersVal]++;
			destPastOthersCount[destVal][pastVal][othersVal]++;
			pastOthersCount[pastVal][othersVal]++;
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
 	 *  agents will contribute to single pdfs, and all are assumed
 	 *  to have other info contributors at same offsets.
 	 * 
	 * @param states
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. src i-j, dest i: transfer is j cells to the right) 
	 * @param otherSourcesToDestOffsets offsets of the other information contributors.
	 *        (i.e. offsets from each other information source to the destination -
	 *        offset is signed the same way as j!)
	 *        othersOffsets is permitted to include j, it will be ignored.
	 */
	public void addObservations(int states[][], int j, int otherSourcesToDestOffsets[]) {
		addObservations(states, j, otherSourcesToDestOffsets, false);
	}
	private void addObservations(int states[][], int j, int otherSourcesToDestOffsets[], boolean cleanedOthers) {
		
		int[] cleanedOthersOffsets;
		if (cleanedOthers) {
			cleanedOthersOffsets = otherSourcesToDestOffsets;
		} else {
			cleanedOthersOffsets = cleanOffsetOthers(otherSourcesToDestOffsets, j, k > 0);
			// This call made redundant by cleanOffsetOthers:
			// confirmEnoughOffsetOthers(othersOffsets, j);
		}
		
		int rows = states.length;
		int columns = states[0].length;
		// increment the count of observations:
		observations += (rows - startObservationTime)*columns; 
		
		// Initialise and store the current previous value for each column
		int[] pastVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			pastVal[c] = 0;
			for (int p = 0; p < k; p++) {
				pastVal[c] *= base;
				pastVal[c] += states[p][c];
			}
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, othersVal;
		for (int r = startObservationTime; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				destVal = states[r][c];
				sourceVal = states[r-1][(c-j+columns) % columns];
				othersVal = 0;
				for (int o = 0; o < cleanedOthersOffsets.length; o++) {
					// Include this other contributor
					othersVal *= base;
					othersVal += states[r-1][(c-cleanedOthersOffsets[o]+columns) % columns];
				}
				sourceDestPastOthersCount[sourceVal][destVal][pastVal[c]][othersVal]++;
				sourcePastOthersCount[sourceVal][pastVal[c]][othersVal]++;
				destPastOthersCount[destVal][pastVal[c]][othersVal]++;
				pastOthersCount[pastVal[c]][othersVal]++;
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
 	 * Add observations for a single source-destination pair of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states
	 */
	public void addObservations(int states[][], int destCol, int sourceCol, int[] othersAbsolute) {
		addObservations(states, destCol, sourceCol, othersAbsolute, false);
	}
	private void addObservations(int states[][], int destCol, int sourceCol, int[] othersAbsolute, boolean cleanedOthers) {

		int[] cleanedOthersAbsolute;
		if (cleanedOthers) {
			cleanedOthersAbsolute = othersAbsolute;
		} else {
			cleanedOthersAbsolute = cleanAbsoluteOthers(othersAbsolute, destCol,
					sourceCol, k > 0);
			// This call made redundant by cleanAbsoluteOthers:
			// confirmEnoughAbsoluteOthers(othersAbsolute, destCol, sourceCol);
		}
		
		int rows = states.length;
		// increment the count of observations:
		observations += (rows - startObservationTime); 
		
		// Initialise and store the current previous value for each column
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[p][destCol];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, othersVal;
		for (int r = startObservationTime; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = states[r][destCol];
			sourceVal = states[r-1][sourceCol];
			othersVal = 0;
			for (int o = 0; o < cleanedOthersAbsolute.length; o++) {
				// Include this other contributor
				othersVal *= base;
				othersVal += states[r-1][cleanedOthersAbsolute[o]];
			}
			sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal]++;
			sourcePastOthersCount[sourceVal][pastVal][othersVal]++;
			destPastOthersCount[destVal][pastVal][othersVal]++;
			pastOthersCount[pastVal][othersVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[states[r-k][destCol]];
				pastVal *= base;
				pastVal += states[r][destCol];
			}
		}
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
		for (int othersVal = 0; othersVal < this.base_power_num_others; othersVal++) {
			for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
				for (int destVal = 0; destVal < base; destVal++) {
					for (int sourceVal = 0; sourceVal < base; sourceVal++) {
						// Compute TE contribution:
						if (sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal] != 0) {
							/* Double check: should never happen
							if ((sourcePastCount[sourceVal][pastVal][othersVal] == 0) ||
								(destPastCount[destVal][pastVal][othersVal] == 0) ||
								(pastCount[pastVal][othersVal] == 0)) {
								throw new RuntimeException("one subcount was zero!!");
							}
							*/
							// compute p(source,dest,past)
							double p_source_dest_past_others = (double)
								sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal] / (double) observations;
							
							double logTerm = ((double) sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal] / (double) sourcePastOthersCount[sourceVal][pastVal][othersVal]) /
							 	((double) destPastOthersCount[destVal][pastVal][othersVal] / (double) pastOthersCount[pastVal][othersVal]);
							double localValue = Math.log(logTerm) / log_2;
							teCont = p_source_dest_past_others * localValue;
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
		}
		
		average = te;
		std = Math.sqrt(meanSqLocals - average * average);
		return te;
	}
	
	/**
	 * Compute the significance of obtaining the given average TE from the given observations
	 * 
	 * 	This is as per Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128.
	 *  except that we've using conditional/complete TE here.
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
				for (int othersVal = 0; othersVal < this.base_power_num_others; othersVal++) {
					numberOfSamples += sourcePastOthersCount[sourceVal][pastVal][othersVal];
				}
			}
			// Now add all of these as unordered observations:
			MatrixUtils.fill(sourceValues, sourceVal, t_s, numberOfSamples);
			t_s += numberOfSamples;
		}
		
		// And construct unordered (dest,past,others) tuples.
		// It doesn't matter that we've appeared to destroy the ordering here because
		//  the joint distribution destPastOthersCount is actually preserved in
		//  our construction of pastVal and destValues and othersValues together here.
		int[] destValues = new int[observations];
		int[] pastValues = new int[observations];
		int[] othersValues = new int[observations];
		int t_d = 0;
		int t_p = 0;
		int t_o = 0;
		for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
			for (int othersVal = 0; othersVal < this.base_power_num_others; othersVal++) {
				// Add in pastOthersCount[pastVal][othersVal] dummy past values
				MatrixUtils.fill(pastValues, pastVal, t_p,
						pastOthersCount[pastVal][othersVal]);
				t_p += pastOthersCount[pastVal][othersVal];
				// Add in pastOthersCount[pastVal][othersVal] dummy others values
				MatrixUtils.fill(othersValues, othersVal, t_o,
						pastOthersCount[pastVal][othersVal]);
				t_o += pastOthersCount[pastVal][othersVal];
				for (int destVal = 0; destVal < base; destVal++) {
					MatrixUtils.fill(destValues, destVal, t_d,
						destPastOthersCount[destVal][pastVal][othersVal]);
					t_d += destPastOthersCount[destVal][pastVal][othersVal];
				}
			}
		}
		
		// Construct new source orderings based on the source probabilities only
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(observations, numPermutationsToCheck);

		ConditionalTransferEntropyCalculator cte = newInstance(base, k, numOtherInfoContributors);
		cte.initialise();
		cte.observations = observations;
		cte.pastOthersCount = pastOthersCount;
		cte.destPastOthersCount = destPastOthersCount;
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the source
			int[] newSourceData = MatrixUtils.extractSelectedTimePoints(sourceValues, newOrderings[p]);
			// compute the joint probability distributions
			MatrixUtils.fill(cte.sourceDestPastOthersCount, 0);
			MatrixUtils.fill(cte.sourcePastOthersCount, 0);
			for (int t = 0; t < observations; t++) {
				cte.sourcePastOthersCount[newSourceData[t]][pastValues[t]][othersValues[t]]++;
				cte.sourceDestPastOthersCount[newSourceData[t]][destValues[t]]
				                            [pastValues[t]][othersValues[t]]++;
			}
			// And get a TE value for this realisation:
			double newTe = cte.computeAverageLocalOfObservations();
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
	 * Computes local complete transfer entropy for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only
	 *  
	 * @param states 2D multivariate time series of states
	 * @param j - number of columns to compute transfer entropy across
	 *  i.e. offset of destination from the source
	 * 	(i.e. src i-j, dest i: transfer is j cells to the right) 
	 * @param otherSourcesToDestOffsets offsets of the other information contributors.
	 *        (i.e. offsets from each other information source to the destination -
	 *        offset is signed the same way as j!)
	 *        othersOffsets is permitted to include j, it will be ignored.
	 * @return
	 */
	public double[][] computeLocalFromPreviousObservations
		(int states[][], int j, int otherSourcesToDestOffsets[]){
		
		return computeLocalFromPreviousObservations(states, j, otherSourcesToDestOffsets, false);
	}
	private double[][] computeLocalFromPreviousObservations
		(int states[][], int j, int othersOffsets[], boolean cleanedOthers){
		
		int[] cleanedOthersOffsets;
		if (cleanedOthers) {
			cleanedOthersOffsets = othersOffsets;
		} else {
			cleanedOthersOffsets = cleanOffsetOthers(othersOffsets, j, k > 0);
			// This call made redundant by cleanOffsetOthers:
			// confirmEnoughOffsetOthers(othersOffsets, j);
		}

		int rows = states.length;
		int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localTE = new double[rows][columns];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int[] pastVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			pastVal[c] = 0;
			for (int p = 0; p < k; p++) {
				pastVal[c] *= base;
				pastVal[c] += states[p][c];
			}
		}
		int destVal, sourceVal, othersVal;
		double logTerm;
		for (int r = startObservationTime; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				destVal = states[r][c];
				sourceVal = states[r-1][(c-j+columns) % columns];
				othersVal = 0;
				for (int o = 0; o < cleanedOthersOffsets.length; o++) {
					// Include this other contributor
					othersVal *= base;
					othersVal += states[r-1][(c-cleanedOthersOffsets[o]+columns) % columns];
				}

				// Now compute the local value
				logTerm = ((double) sourceDestPastOthersCount[sourceVal][destVal][pastVal[c]][othersVal] / (double) sourcePastOthersCount[sourceVal][pastVal[c]][othersVal]) /
			 		((double) destPastOthersCount[destVal][pastVal[c]][othersVal] / (double) pastOthersCount[pastVal[c]][othersVal]);
				localTE[r][c] = Math.log(logTerm) / log_2;
				average += localTE[r][c];
				if (localTE[r][c] > max) {
					max = localTE[r][c];
				} else if (localTE[r][c] < min) {
					min = localTE[r][c];
				}
				// Update the previous value:
				if (k > 0) {
					pastVal[c] -= maxShiftedValue[states[r-k][c]];
					pastVal[c] *= base;
					pastVal[c] += states[r][c];
				}
			}
		}		

		average = average/(double) (columns * (rows - startObservationTime));
		
		return localTE;
	}
	
	/**
	 * Computes local active information storage for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations
		(int states[][], int destCol, int sourceCol, int[] othersAbsolute){
		
		return computeLocalFromPreviousObservations(states, destCol, sourceCol, othersAbsolute, false);
	}
	private double[] computeLocalFromPreviousObservations
		(int states[][], int destCol, int sourceCol, int[] othersAbsolute, boolean cleanedOthers){

		int[] cleanedOthersAbsolute;
		if (cleanedOthers) {
			cleanedOthersAbsolute = othersAbsolute;
		} else {
			cleanedOthersAbsolute = cleanAbsoluteOthers(othersAbsolute,
					destCol, sourceCol, k > 0);
			// This call made redundant by cleanOffsetOthers:
			// confirmEnoughAbsoluteOthers(othersAbsolute, destCol, sourceCol);
		}
		
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
		int destVal, sourceVal, othersVal;
		double logTerm;
		for (int r = startObservationTime; r < rows; r++) {
			destVal = states[r][destCol];
			sourceVal = states[r-1][sourceCol];
			othersVal = 0;
			for (int o = 0; o < cleanedOthersAbsolute.length; o++) {
				// Include this other contributor
				othersVal *= base;
				othersVal += states[r-1][cleanedOthersAbsolute[o]];
			}
			// Now compute the local value
			logTerm = ((double) sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal] / (double) sourcePastOthersCount[sourceVal][pastVal][othersVal]) /
		 		((double) destPastOthersCount[destVal][pastVal][othersVal] / (double) pastOthersCount[pastVal][othersVal]);
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
	 * Standalone routine to 
	 * compute local transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method to be called for homogeneous agents only
	 * 
	 * @param states - 2D array of states
	 * @param j - TE across j cells to the right
	 * @param otherSourcesToDestOffsets - column offsets from other causal info contributors
	 *  to the destination
	 * @return
	 */
	public double[][] computeLocal(int states[][], int j, int[] otherSourcesToDestOffsets) {
		
		initialise();
		int[] cleanedOthersOffsets = cleanOffsetOthers(otherSourcesToDestOffsets, j, k > 0);
		addObservations(states, j, cleanedOthersOffsets, true);
		return computeLocalFromPreviousObservations(states, j, cleanedOthersOffsets, true);
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
	 * @param otherSourcesToDestOffsets - column offsets from other causal info contributors
	 *  to the destination
	 * @return
	 */
	public double computeAverageLocal(int states[][], int j, int[] otherSourcesToDestOffsets) {
		
		initialise();
		addObservations(states, j, otherSourcesToDestOffsets);
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
	 * @param othersAbsolute - column indices for other causal info contributors
	 * @return
	 */
	public double[] computeLocal(int states[][], int destCol, int sourceCol, int[] othersAbsolute) {
		
		initialise();
		int[] cleanedOthers = cleanAbsoluteOthers(othersAbsolute, destCol,
				sourceCol, k > 0);
		addObservations(states, destCol, sourceCol, cleanedOthers, true);
		return computeLocalFromPreviousObservations(states, destCol, sourceCol, cleanedOthers, true);
	}

	/**
	 * Standalone routine to 
	 * compute average local transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Returns the average
	 * This method suitable for heterogeneous agents
	 * 
	 * @param states - 2D array of states
	 * @param destCol - column index for the destination agent
	 * @param sourceCol - column index for the source agent
	 * @param othersAbsolute - column indices for other causal info contributors
	 * @return
	 */
	public double computeAverageLocal(int states[][], int destCol, int sourceCol, int[] othersAbsolute) {
		
		initialise();
		addObservations(states, destCol, sourceCol, othersAbsolute);
		return computeAverageLocalOfObservations();
	}
	
	/**
	 * Counts the information contributors to this node which
	 * are not equal to the source to dest offset j or the node itself (offset 0,
	 * node itself not included only when removeDest is set to true)
	 * 
	 * @param otherSourcesToDestOffsets array of offsets of the destination from each source
	 * @param j offset of the destination from the source
	 * @param removeDest remove the destination itself from the count
	 *     of offset others.
	 * @return
	 */
	public static int countOfOffsetOthers(int[] otherSourcesToDestOffsets, int j,
			boolean removeDest) {
		int countOfOthers = 0;
		for (int index = 0; index < otherSourcesToDestOffsets.length; index++) {
			if ((otherSourcesToDestOffsets[index] != j) && 
				((otherSourcesToDestOffsets[index] != 0) || !removeDest)) {
				countOfOthers++;
			}
		}
		return countOfOthers;
	}
	
	/**
	 * Counts the information contributors to the dest which
	 * are not equal to src or the node itself (offset 0,
	 * node itself not included only when removeDest is set to true)
	 * 
	 * @param others
	 * @param dest
	 * @param src
	 * @param removeDest remove the destination itself from the count
	 *     of absolute others.
	 * @return
	 */
	public static int countOfAbsoluteOthers(int[] others, int dest, int src,
			boolean removeDest) {
		int countOfOthers = 0;
		for (int index = 0; index < others.length; index++) {
			if ((others[index] != src) &&
				((others[index] != dest) || !removeDest)) {
				countOfOthers++;
			}
		}
		return countOfOthers;
	}

	/**
	 * Check that the supplied array of offsets as other info 
	 * contributors is long enough compared to our expectation
	 * 
	 * @param othersOffsets array of offsets from each source to the destination
	 * @param j offset from the source to the destination
	 * @param removeDest remove the destination itself from the count
	 *     of absolute others.
	 * @return
	 */
	public boolean confirmEnoughOffsetOthers(int[] othersOffsets, int j,
			boolean removeDest) {
		if (countOfOffsetOthers(othersOffsets, j, removeDest) !=
				numOtherInfoContributors) {
			throw new RuntimeException("Incorrect number of others in offsets");
		}
		return true;
	}

	/**
	 * Check that the supplied array of absolutes as other info 
	 * contributors is long enough compared to our expectation
	 * 
	 * @param othersAbsolute
	 * @param dest
	 * @param src
	 * @param removeDest remove the destination itself from the count
	 *     of absolute others.
	 * @return
	 */
	public boolean confirmEnoughAbsoluteOthers(int[] othersAbsolute, int dest,
			int src, boolean removeDest) {
		if (countOfAbsoluteOthers(othersAbsolute, dest, src, removeDest) !=
				numOtherInfoContributors) {
			throw new RuntimeException("Incorrect number of others in absolutes");
		}
		return true;
	}
	
	/**
	 * Returns the information contributors to this node which
	 * are not equal to the offset j or the node itself (offset 0,
	 * removed only if removeDest is set to true).
	 * Checks that there are enough other information contributors.
	 * 
	 * @param othersOffsets array of offsets from each source to the destination
	 * @param j offset from the source to the destination
	 * @param removeDest Remove the destination itself from the cleaned
	 *           other sources (if it is there). Should not be done
	 *           if k == 0 (because then the destination is not included
	 *           in the past history)
	 * @return
	 */
	public int[] cleanOffsetOthers(int[] othersOffsets, int j, boolean removeDest) {
		int[] cleaned = new int[numOtherInfoContributors];
		int countOfOthers = 0;
		for (int index = 0; index < othersOffsets.length; index++) {
			if ((othersOffsets[index] != j) &&
				((othersOffsets[index] != 0) || !removeDest)) {
				// Add this candidate source to the cleaned sources
				if (countOfOthers == numOtherInfoContributors) {
					// We've already taken all the other info 
					//  contributors we expected
					countOfOthers++;
					break;
				}
				cleaned[countOfOthers] = othersOffsets[index];
				countOfOthers++;
			}
		}
		if (countOfOthers < numOtherInfoContributors) {
			throw new RuntimeException("Too few others in offsets");
		} else if (countOfOthers > numOtherInfoContributors) {
			throw new RuntimeException("Too many others in offsets");
		}
		return cleaned;
	}
	
	/**
	 * Returns the information contributors to the dest which
	 * are not equal to src or the node itself (offset 0,
	 * removed only if removeDest is true).
	 * Checks that there are enough other information contributors.
	 * 
	 * @param others
	 * @param dest
	 * @param src
	 * @param removeDest Remove the destination itself from the cleaned
	 *           other sources (if it is there). Should not be done
	 *           if k == 0 (because then the destination is not included
	 *           in the past history)
	 * @return
	 */
	public int[] cleanAbsoluteOthers(int[] others, int dest, int src,
			boolean removeDest) {
		int[] cleaned = new int[numOtherInfoContributors];
		int countOfOthers = 0;
		for (int index = 0; index < others.length; index++) {
			if ((others[index] != src) &&
				((others[index] != dest) || !removeDest)) {
				// Add this candidate source to the cleaned sources
				if (countOfOthers == numOtherInfoContributors) {
					// We've already taken all the other info 
					//  contributors we expected
					countOfOthers++;
					break;
				}
				cleaned[countOfOthers] = others[index];
				countOfOthers++;
			}
		}
		if (countOfOthers < numOtherInfoContributors) {
			throw new RuntimeException("Too few others in absolutes");
		} else if (countOfOthers > numOtherInfoContributors) {
			throw new RuntimeException("Too many others in absolutes");
		}
		return cleaned;
	}
}
