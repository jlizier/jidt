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

import infodynamics.utils.AnalyticMeasurementDistribution;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Implements <b>transfer entropy</b>
 * for univariate discrete time-series data.
 * That is, it is applied to <code>int[]</code> data, indexed
 * by time.
 * See Schreiber below for the definition of transfer entropy,
 * and Lizier et al. for the definition of local transfer entropy.
 * Specifically, this class implements the pairwise or <i>apparent</i>
 * transfer entropy; i.e. we compute the transfer that appears to
 * come from a single source variable, without examining any other
 * potential sources
 * (see Lizier et al, PRE, 2008).</p>
 *  
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * <ol>
 * 		<li>Construct the calculator via {@link #TransferEntropyCalculatorDiscrete(int, int)}
 * 			or {@link #TransferEntropyCalculatorDiscrete(int, int, int)}
 * 			or {@link #TransferEntropyCalculatorDiscrete(int, int, int, int, int, int)};</li>
 *		<li>Initialise the calculator using
 *			{@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			the set of {@link #addObservations(int[], int[])} methods, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average TE: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local TE values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local TE values for a specific set of samples: e.g.
 * 				{@link #computeLocalFromPreviousObservations(int[], int[])} etc.</li>
 * 				<li>the distribution of TE values under the null hypothesis
 * 					of no relationship between source and
 * 					destination values: {@link #computeSignificance(int)} or
 * 					{@link #computeSignificance(int[][])}.</li>
 * 			</ul>
 * 		</li>
 * 		<li>As an alternative to steps 3 and 4, the user may undertake
 * 			standalone computation from a single set of observations, via
 *  		e.g.: {@link #computeLocal(int[], int[])} or
 *  		{@link #computeAverageLocal(int[][], int)}.</li>
 * 		<li>
 * 		Return to step 2 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 * </ul>
 * 
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 */
public class TransferEntropyCalculatorDiscrete extends ContextOfPastMeasureCalculatorDiscrete 
	implements ChannelCalculatorDiscrete, AnalyticNullDistributionComputer {

	/**
	 * Counts of (source,dest_next,dest_embedded_past) tuples
	 */
	protected int[][][] sourceNextPastCount = null;	// count for (source[n],dest[n+1],dest[n]^k) tuples
	/**
	 * Counts of (source,dest_embedded_past) tuples
	 */
	protected int[][] sourcePastCount = null;			// count for (source[n],dest[n]^k) tuples
	/**
	 * Whether to assume periodic boundary conditions for channels across
	 *  the boundary of the multidimensional
	 *  calls supplying observations, e.g.
	 *  {@link #addObservations(int[][], int)} calls
	 */
	protected boolean periodicBoundaryConditions = true;

	/**
	 * Embedding delay for the destination variable,
	 *  i.e. time lag between each sample in the past
	 */
	protected int destEmbeddingDelay = 1;

	/**
	 * Embedding length of the source variable.
	 * This is "l" in Schreiber's notation.
	 */
	protected int sourceHistoryEmbedLength = 1;
	
	/**
	 * Embedding delay for the source variable,
	 *  i.e. time lag between each sample in the past
	 */
	protected int sourceEmbeddingDelay = 1;

	/**
	 * Source-destination delay to consider the information transfer across
	 */
	protected int delay = 1;
	
	/**
	 * A cached value of base^sourceHistoryEmbedLength
	 */
	protected int base_power_l = 1;
	
	/**
	 * A cached value of each discrete value left shifted (in "base" counting) by (sourceHistoryEmbedLength-1).
	 */
	protected int[] maxShiftedSourceValue = null; // states * (base^(sourceHistoryEmbedLength-1))

	/**
	 * First time step at which we can take an observation
	 *  (needs to account for an embedding in the previous steps)
	 */
	protected int startObservationTime = 1;
	
	/**
	 * Tracks whether the measure has been computed since the last initialisation
	 */
	protected boolean estimateComputed = false;
	
	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base
	 * @param destHistoryEmbedLength
	 * 
	 * @return a new TransferEntropyCalculator object
	 * @deprecated
	 */
	public static TransferEntropyCalculatorDiscrete newInstance(int base, int destHistoryEmbedLength) {
		
		return new TransferEntropyCalculatorDiscrete(base, destHistoryEmbedLength);

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
	 * Create a new TE calculator for the given base and destination history embedding length
	 *  (leave the other embedding parameters as default)
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param destHistoryEmbedLength embedded history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 */
	public TransferEntropyCalculatorDiscrete(int base, int destHistoryEmbedLength) {

		this(base, destHistoryEmbedLength, 1, 1, 1, 1);
	}

	/**
	 * Create a new TE calculator for the given base, destination and source history embedding lengths.
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param destHistoryEmbedLength embedded history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 * @param sourceHistoryEmbeddingLength embedded history length of the source to include -
	 *        this is l in Schreiber's notation.
	 */
	public TransferEntropyCalculatorDiscrete(int base, int destHistoryEmbedLength, int sourceHistoryEmbeddingLength) {

		this(base, destHistoryEmbedLength, 1, sourceHistoryEmbeddingLength, 1, 1);
	}

	/**
	 * Create a new TE calculator for the given base, destination and source history embedding lengths
	 *  and delays.
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param destHistoryEmbedLength embedded history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 * @param destEmbeddingDelay embedding delay of the destination for conditioning on -
	 *        this is the delay between each of the k samples from the past history
	 * @param sourceHistoryEmbeddingLength embedded history length of the source to include -
	 *        this is l in Schreiber's notation.
	 * @param sourceEmbeddingDelay embedding delay of the source -
	 *        this is the delay between each of the l samples from the past history
	 * @param delay source-destination delay to consider the information transfer across
	 * 		  (should be >= 0, default is 1)
	 */
	public TransferEntropyCalculatorDiscrete(int base, int destHistoryEmbedLength, int destEmbeddingDelay,
			int sourceHistoryEmbeddingLength, int sourceEmbeddingDelay, int delay) {

		super(base, destHistoryEmbedLength);
		this.destEmbeddingDelay = destEmbeddingDelay;
		if (sourceHistoryEmbeddingLength <= 0) {
			throw new RuntimeException("Cannot have source embedding length of zero or less");
		}
		this.sourceHistoryEmbedLength = sourceHistoryEmbeddingLength;
		this.sourceEmbeddingDelay = sourceEmbeddingDelay;
		this.delay = delay;
		base_power_l = MathsUtils.power(base, sourceHistoryEmbedLength);
		
		// Check that we can convert the history value into an integer ok: 
		if (sourceHistoryEmbedLength > Math.log(Integer.MAX_VALUE) / log_base) {
			throw new RuntimeException("Base and source history combination too large");
		}

		// Create constants for tracking sourceValues
		maxShiftedSourceValue = new int[base];
		for (int v = 0; v < base; v++) {
			maxShiftedSourceValue[v] = v * MathsUtils.power(base, sourceHistoryEmbedLength-1);
		}

		// Create storage for extra counts of observations
		try {
			sourceNextPastCount = new int[base_power_l][base][base_power_k];
			sourcePastCount = new int[base_power_l][base_power_k];
		} catch (OutOfMemoryError e) {
			// Allow any Exceptions to be thrown, but catch and wrap
			//  Error as a RuntimeException
			throw new RuntimeException("Requested memory for the base " +
					base + ", k=" + k + ", l=" + sourceHistoryEmbedLength +
					" is too large for the JVM at this time", e);
		}
		
		// Which time step do we start taking observations from?
		// These two integers represent the earliest next time step, in the cases where the destination
		// embedding itself determines where we can start taking observations, or
		// the case where the source embedding plus delay is longer and so determines
		// where we can start taking observations.
		int startTimeBasedOnDestPast = (k-1)*destEmbeddingDelay + 1;
		int startTimeBasedOnSourcePast = (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay + delay;
		startObservationTime = Math.max(startTimeBasedOnDestPast, startTimeBasedOnSourcePast);

	}

	@Override
	public void initialise(){
		super.initialise();
		estimateComputed = false;
		
		MatrixUtils.fill(sourceNextPastCount, 0);
		MatrixUtils.fill(sourcePastCount, 0);
	}
	
	@Override
	public void addObservations(int[] source, int[] dest) {
		addObservations(source, dest, 0, dest.length-1);
	}
	
	/**
 	 * Add observations for a single source-destination pair
 	 *  to our estimates of the pdfs.
	 * Start and end time are the (inclusive) indices within which to add the observations.
	 * The start time is from the earliest of the k historical values of the destination (inclusive),
	 *  the end time is the last destination time point to add in.
	 *  
	 * @param source source time-series
	 * @param dest destination time-series. 
	 *  Must be same length as source
	 * @param startTime earliest time that we may extract embedded history from
	 * @param endTime last destination (next) time point to add in
	 * 
	 */
	public void addObservations(int[] source, int[] dest, int startTime, int endTime) {
		if ((endTime - startTime) - startObservationTime + 1 <= 0) {
			// No observations to add
			return;
		}
		if ((endTime >= dest.length) || (endTime >= source.length)) {
			throw new ArrayIndexOutOfBoundsException(
					String.format("endTime (%d) must be <= length of input arrays (dest: %d, source: %d)",
							endTime, dest.length, source.length));
		}
		// increment the count of observations:
		observations += (endTime - startTime) - startObservationTime + 1; 
		
		// Initialise and store the current previous values;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[] pastVal = new int[destEmbeddingDelay];
		for (int d = 0; d < destEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			pastVal[d] = 0;
			for (int p = 0; p < k-1; p++) {
				pastVal[d] += dest[startTime + startObservationTime + d - 1
				                   - (k-1)*destEmbeddingDelay
				                   + p*destEmbeddingDelay];
				pastVal[d] *= base;
			}
		}
		// Next for the source:
		int[] sourcePastVal = new int[sourceEmbeddingDelay];
		for (int d = 0; d < sourceEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			sourcePastVal[d] = 0;
			for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
				sourcePastVal[d] += source[startTime + startObservationTime + d - delay
				                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
				                    + p*sourceEmbeddingDelay];
				sourcePastVal[d] *= base;
			}
		}
		
		// 1. Count the tuples observed
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
		for (int r = startTime + startObservationTime; r <= endTime; r++) {
			// First update the embedding values for the current
			//  phases of the embeddings:
			if (k > 0) {
				pastVal[destEmbeddingPhase] += dest[r-1];
			}
			sourcePastVal[sourceEmbeddingPhase] += source[r-delay];
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = dest[r];
			int thisPastVal = pastVal[destEmbeddingPhase];
			int thisSourceVal = sourcePastVal[sourceEmbeddingPhase];
			sourceNextPastCount[thisSourceVal][destVal][thisPastVal]++;
			sourcePastCount[thisSourceVal][thisPastVal]++;
			nextPastCount[destVal][thisPastVal]++;
			pastCount[thisPastVal]++;
			nextCount[destVal]++;
			// Now, update the combined embedding values and phases,
			//  for this phase we back out the oldest value which we'll no longer need:
			if (k > 0) {
				pastVal[destEmbeddingPhase] -= maxShiftedValue[dest[r-1-(k-1)*destEmbeddingDelay]];
				pastVal[destEmbeddingPhase] *= base; // and shift the others up
			}
			sourcePastVal[sourceEmbeddingPhase] -=
					maxShiftedSourceValue[
					    source[r-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay]];
			sourcePastVal[sourceEmbeddingPhase] *= base; // and shift the others up
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
		}
	}

	/**
 	 * Add observations for a single source-destination pair 
 	 *  to our estimates of the pdfs.
 	 *  
	 * @param source source time-series
	 * @param dest destination time-series. 
	 *  Must be same length as source
	 * @param valid time-series of whether the signals
	 *  at the given time should be considered valid 
	 *  and added to our PDFs. We don't include any embedding vectors which
	 *  stretch across any invalid points, even if these invalid points
	 *  are not specifically sampled for the embedding vector.
	 */
	public void addObservations(int[] source, int[] dest, boolean[] valid) {
		
		int rows = dest.length;
		if (dest.length - startObservationTime <= 0) {
			// No observations to add
			return;
		}
		
		// Initialise and store the current previous values;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[] pastVal = new int[destEmbeddingDelay];
		for (int d = 0; d < destEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			pastVal[d] = 0;
			for (int p = 0; p < k-1; p++) {
				pastVal[d] += dest[startObservationTime + d - 1
				                   - (k-1)*destEmbeddingDelay
				                   + p*destEmbeddingDelay];
				pastVal[d] *= base;
			}
		}
		// We can take an observation if timeSinceLastDestInvalid >= minDestLengthRequired
		int minDestLengthRequired = (k>0) ? (k-1)*destEmbeddingDelay + 1 : 0;
		// And initialise the time since last valid dest observation:
		int timeSinceLastDestInvalid = minDestLengthRequired;
		for (int t = startObservationTime - 1; t >= 0; t--) {
			if (!valid[t]) {
				timeSinceLastDestInvalid = startObservationTime - t - 1;
				break;
			}
		}
		
		// Next for the source:
		int[] sourcePastVal = new int[sourceEmbeddingDelay];
		for (int d = 0; d < sourceEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			sourcePastVal[d] = 0;
			for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
				sourcePastVal[d] += source[startObservationTime + d - delay
				                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
				                    + p*sourceEmbeddingDelay];
				sourcePastVal[d] *= base;
			}
		}
		// We can take an observation if timeSinceLastSourceInvalid >= minSourceLengthRequired
		int minSourceLengthRequired = (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay + 1;
		// And initialise the time since last valid source observation:
		int timeSinceLastSourceInvalid = minSourceLengthRequired;
		for (int t = startObservationTime - delay; t >= 0; t--) {
			if (!valid[t]) {
				timeSinceLastSourceInvalid = startObservationTime - t - 1;
				break;
			}
		}

		// 1. Count the tuples observed
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
		for (int r = startObservationTime; r < rows; r++) {
			timeSinceLastDestInvalid++;
			timeSinceLastSourceInvalid++;
			// First update the embedding values for the current
			//  phases of the embeddings:
			// Pre-condition now:
			//  timeSinceLastDestInvalid holds the time from r back to
			//   the last invalid destination point (not including current destination)
			//  timeSinceLastSourceInvalid holds the time from r back to
			//   the last invalid source point (not including new source value)

			// Update embedded values
			if (k > 0) {
				pastVal[destEmbeddingPhase] += dest[r-1];
			}
			sourcePastVal[sourceEmbeddingPhase] += source[r-delay];
			
			// Now check validity of new entrants here:
			if (!valid[r]) {
				timeSinceLastDestInvalid = 0;
			}
			if (!valid[r-delay]) {
				timeSinceLastSourceInvalid = 0;
			}
			if ((timeSinceLastDestInvalid > minDestLengthRequired) &&
				    (timeSinceLastSourceInvalid >= minSourceLengthRequired)) {
				// We have enough of both source and dest valid to continue
				// Note the >= on source only (because dest must have next value as 
				//  valid as well.
				
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				destVal = dest[r];
				int thisPastVal = pastVal[destEmbeddingPhase];
				int thisSourceVal = sourcePastVal[sourceEmbeddingPhase];
				sourceNextPastCount[thisSourceVal][destVal][thisPastVal]++;
				sourcePastCount[thisSourceVal][thisPastVal]++;
				nextPastCount[destVal][thisPastVal]++;
				pastCount[thisPastVal]++;
				nextCount[destVal]++;
				observations++; // Need to increment number of observations explicitly here
			}
			// Now, update the combined embedding values and phases,
			//  for this phase we back out the oldest value which we'll no longer need:
			if (k > 0) {
				pastVal[destEmbeddingPhase] -= maxShiftedValue[dest[r-1-(k-1)*destEmbeddingDelay]];
				pastVal[destEmbeddingPhase] *= base; // and shift the others up
			}
			sourcePastVal[sourceEmbeddingPhase] -=
					maxShiftedSourceValue[
					    source[r-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay]];
			sourcePastVal[sourceEmbeddingPhase] *= base; // and shift the others up
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
		}
	}

	/**
 	 * Add observations in to our estimates of the PDFs,
 	 * from a multivariate time-series.
 	 * This call suitable only for homogeneous variables, as all
 	 *  variable pairs separated by j column will contribute to the PDFs.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i: we
	 *  compute transfer is j cells to the right, using observations
	 *  across all column pairs separated by j) 
	 */
	public void addObservations(int states[][], int j) {

		int timeSteps = states.length;
		if (timeSteps - startObservationTime <= 0) {
			// No observations to add
			return;
		}
		int variables = states[0].length;
		// increment the count of observations:
		if (periodicBoundaryConditions) {
			observations += (timeSteps - startObservationTime)*variables; 			
		} else {
			observations += (timeSteps - startObservationTime)*(variables - Math.abs(j));
		}
		
		// Initialise and store the current previous values for each column;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[][] pastVal = new int[variables][destEmbeddingDelay];
		for (int c = 0; c < variables; c++) {
			for (int d = 0; d < destEmbeddingDelay; d++) {
				// Compute the current previous values for
				//  phase d of the embedding delay, but leave
				//  out the most recent value (we'll add those in
				//  in the main loop)
				pastVal[c][d] = 0;
				for (int p = 0; p < k-1; p++) {
					pastVal[c][d] += states[startObservationTime + d - 1
					                   - (k-1)*destEmbeddingDelay
					                   + p*destEmbeddingDelay][c];
					pastVal[c][d] *= base;
				}
			}
		}
		// Next for the source:
		int[][] sourcePastVal = new int[variables][sourceEmbeddingDelay];
		for (int c = 0; c < variables; c++) {
			int sourceVariable = c-j;
			if ((sourceVariable < 0) || (sourceVariable >= variables)) {
				// Source variable is out of bounds unless we are using periodic boundary conditions
				if (periodicBoundaryConditions) {
					sourceVariable = (sourceVariable+variables) % variables;
				} else {
					// Don't add this to our observations
					continue;
				}
			}
			for (int d = 0; d < sourceEmbeddingDelay; d++) {
				// Compute the current previous values for
				//  phase d of the embedding delay, but leave
				//  out the most recent value (we'll add those in
				//  in the main loop)
				sourcePastVal[c][d] = 0;
				for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
					sourcePastVal[c][d] += states[startObservationTime + d - delay
					                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
					                    + p*sourceEmbeddingDelay][sourceVariable];
					sourcePastVal[c][d] *= base;
				}
			}
		}
		
		// 1. Count the tuples observed
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
		for (int r = startObservationTime; r < timeSteps; r++) {
			for (int c = 0; c < variables; c++) {
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				int sourceVariable = c-j;
				if ((sourceVariable < 0) || (sourceVariable >= variables)) {
					// Source variable is out of bounds unless we are using periodic boundary conditions
					if (periodicBoundaryConditions) {
						sourceVariable = (sourceVariable+variables) % variables;
					} else {
						// Don't add this to our observations
						continue;
					}
				}
				// First update the embedding values for the current
				//  phases of the embeddings:
				if (k > 0) {
					pastVal[c][destEmbeddingPhase] += states[r-1][c];
				}
				sourcePastVal[c][sourceEmbeddingPhase] += states[r-delay][sourceVariable];
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				destVal = states[r][c];
				int thisPastVal = pastVal[c][destEmbeddingPhase];
				int thisSourceVal = sourcePastVal[c][sourceEmbeddingPhase];
				sourceNextPastCount[thisSourceVal][destVal][thisPastVal]++;
				sourcePastCount[thisSourceVal][thisPastVal]++;
				nextPastCount[destVal][thisPastVal]++;
				pastCount[thisPastVal]++;
				nextCount[destVal]++;
				// Now, update the combined embedding values and phases,
				//  for this phase we back out the oldest value which we'll no longer need:
				if (k > 0) {
					pastVal[c][destEmbeddingPhase] -= maxShiftedValue[states[r-1-(k-1)*destEmbeddingDelay][c]];
					pastVal[c][destEmbeddingPhase] *= base; // and shift the others up
				}
				sourcePastVal[c][sourceEmbeddingPhase] -=
						maxShiftedSourceValue[
						    states[r-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay][sourceVariable]];
				sourcePastVal[c][sourceEmbeddingPhase] *= base; // and shift the others up
			}
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
		}
	}

	/**
 	 * Add observations in to our estimates of the PDFs,
 	 * from a multivariate time-series.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by h rows and j columns
 	 *  will contribute to the PDFs.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param h - number of rows to compute transfer entropy across
	 * 	(i.e. source is in row i-h, dest is column i) 
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i) 
	 */
	public void addObservations(int states[][][], int h, int j) {

		int timeSteps = states.length;
		if (timeSteps - startObservationTime <= 0) {
			// No observations to add
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

		// Initialise and store the current previous values for each variable;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[][][] pastVal = new int[agentRows][agentColumns][destEmbeddingDelay];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				for (int d = 0; d < destEmbeddingDelay; d++) {
					// Compute the current previous values for
					//  phase d of the embedding delay, but leave
					//  out the most recent value (we'll add those in
					//  in the main loop)
					pastVal[r][c][d] = 0;
					for (int p = 0; p < k-1; p++) {
						pastVal[r][c][d] += states[startObservationTime + d - 1
						                   - (k-1)*destEmbeddingDelay
						                   + p*destEmbeddingDelay][r][c];
						pastVal[r][c][d] *= base;
					}
				}
			}
		}
		// Next for the source:
		int[][][] sourcePastVal = new int[agentRows][agentColumns][sourceEmbeddingDelay];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
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
				// Now initialise the embedding
				for (int d = 0; d < sourceEmbeddingDelay; d++) {
					// Compute the current previous values for
					//  phase d of the embedding delay, but leave
					//  out the most recent value (we'll add those in
					//  in the main loop)
					sourcePastVal[r][c][d] = 0;
					for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
						sourcePastVal[r][c][d] += states[startObservationTime + d - delay
						                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
						                    + p*sourceEmbeddingDelay][sourceAgentRow][sourceAgentColumn];
						sourcePastVal[r][c][d] *= base;
					}
				}
			}
		}
		
		// 1. Count the tuples observed
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
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
					// First update the embedding values for the current
					//  phases of the embeddings:
					if (k > 0) {
						pastVal[r][c][destEmbeddingPhase] += states[t-1][r][c];
					}
					sourcePastVal[r][c][sourceEmbeddingPhase] += states[t-delay][sourceAgentRow][sourceAgentColumn];
					// Add to the count for this particular transition:
					// (cell's assigned as above)
					destVal = states[t][r][c];
					int thisPastVal = pastVal[r][c][destEmbeddingPhase];
					int thisSourceVal = sourcePastVal[r][c][sourceEmbeddingPhase];
					sourceNextPastCount[thisSourceVal][destVal][thisPastVal]++;
					sourcePastCount[thisSourceVal][thisPastVal]++;
					nextPastCount[destVal][thisPastVal]++;
					pastCount[thisPastVal]++;
					nextCount[destVal]++;
					// Now, update the combined embedding values and phases,
					//  for this phase we back out the oldest value which we'll no longer need:
					if (k > 0) {
						pastVal[r][c][destEmbeddingPhase] -= maxShiftedValue[states[t-1-(k-1)*destEmbeddingDelay][r][c]];
						pastVal[r][c][destEmbeddingPhase] *= base; // and shift the others up
					}
					sourcePastVal[r][c][sourceEmbeddingPhase] -=
							maxShiftedSourceValue[
							    states[t-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay][sourceAgentRow][sourceAgentColumn]];
					sourcePastVal[r][c][sourceEmbeddingPhase] *= base; // and shift the others up
				}
			}
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
		}
	}

	/**
 	 * Add observations for a single source-destination pair of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to {@link #addObservations(int[][], int)}
 	 *  for computing TE for heterogeneous agents.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param sourceIndex source variable index in states
	 * @param destIndex destination variable index in states
	 */
	public void addObservations(int states[][], int sourceIndex, int destIndex) {

		int timeSteps = states.length;
		if (timeSteps - startObservationTime <= 0) {
			// No observations to add
			return;
		}
		// increment the count of observations:
		observations += timeSteps - startObservationTime; 
		
		// Initialise and store the current previous values;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[] pastVal = new int[destEmbeddingDelay];
		for (int d = 0; d < destEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			pastVal[d] = 0;
			for (int p = 0; p < k-1; p++) {
				pastVal[d] += states[startObservationTime + d - 1
				                   - (k-1)*destEmbeddingDelay
				                   + p*destEmbeddingDelay][destIndex];
				pastVal[d] *= base;
			}
		}
		// Next for the source:
		int[] sourcePastVal = new int[sourceEmbeddingDelay];
		for (int d = 0; d < sourceEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			sourcePastVal[d] = 0;
			for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
				sourcePastVal[d] += states[startObservationTime + d - delay
				                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
				                    + p*sourceEmbeddingDelay][sourceIndex];
				sourcePastVal[d] *= base;
			}
		}
		
		// 1. Count the tuples observed
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
		for (int r = startObservationTime; r < timeSteps; r++) {
			// First update the embedding values for the current
			//  phases of the embeddings:
			if (k > 0) {
				pastVal[destEmbeddingPhase] += states[r-1][destIndex];
			}
			sourcePastVal[sourceEmbeddingPhase] += states[r-delay][sourceIndex];
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = states[r][destIndex];
			int thisPastVal = pastVal[destEmbeddingPhase];
			int thisSourceVal = sourcePastVal[sourceEmbeddingPhase];
			sourceNextPastCount[thisSourceVal][destVal][thisPastVal]++;
			sourcePastCount[thisSourceVal][thisPastVal]++;
			nextPastCount[destVal][thisPastVal]++;
			pastCount[thisPastVal]++;
			nextCount[destVal]++;
			// Now, update the combined embedding values and phases,
			//  for this phase we back out the oldest value which we'll no longer need:
			if (k > 0) {
				pastVal[destEmbeddingPhase] -= maxShiftedValue[states[r-1-(k-1)*destEmbeddingDelay][destIndex]];
				pastVal[destEmbeddingPhase] *= base; // and shift the others up
			}
			sourcePastVal[sourceEmbeddingPhase] -=
					maxShiftedSourceValue[
					    states[r-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay][sourceIndex]];
			sourcePastVal[sourceEmbeddingPhase] *= base; // and shift the others up
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
		}
	}

	/**
 	 * Add observations for a single source-destination pair of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to {@link #addObservations(int[][][], int, int)}
 	 *  for computing TE for heterogeneous agents.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param sourceRowIndex source variable row index in states
	 * @param sourceColumnIndex source variable column index in states
	 * @param destRowIndex destination variable row index in states
	 * @param destColumnIndex destination variable column index in states
	 */
	public void addObservations(int states[][][], int sourceRowIndex, int sourceColumnIndex,
												  int destRowIndex, int destColumnIndex) {

		int timeSteps = states.length;
		if (timeSteps - startObservationTime <= 0) {
			// No observations to add
			return;
		}
		// increment the count of observations:
		observations += timeSteps - startObservationTime; 
		
		// Initialise and store the current previous values;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[] pastVal = new int[destEmbeddingDelay];
		for (int d = 0; d < destEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			pastVal[d] = 0;
			for (int p = 0; p < k-1; p++) {
				pastVal[d] += states[startObservationTime + d - 1
				                   - (k-1)*destEmbeddingDelay
				                   + p*destEmbeddingDelay][destRowIndex][destColumnIndex];
				pastVal[d] *= base;
			}
		}
		// Next for the source:
		int[] sourcePastVal = new int[sourceEmbeddingDelay];
		for (int d = 0; d < sourceEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			sourcePastVal[d] = 0;
			for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
				sourcePastVal[d] += states[startObservationTime + d - delay
				                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
				                    + p*sourceEmbeddingDelay][sourceRowIndex][sourceColumnIndex];
				sourcePastVal[d] *= base;
			}
		}
		
		// 1. Count the tuples observed
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
		for (int r = startObservationTime; r < timeSteps; r++) {
			// First update the embedding values for the current
			//  phases of the embeddings:
			if (k > 0) {
				pastVal[destEmbeddingPhase] += states[r-1][destRowIndex][destColumnIndex];
			}
			sourcePastVal[sourceEmbeddingPhase] += states[r-delay][sourceRowIndex][sourceColumnIndex];
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = states[r][destRowIndex][destColumnIndex];
			int thisPastVal = pastVal[destEmbeddingPhase];
			int thisSourceVal = sourcePastVal[sourceEmbeddingPhase];
			sourceNextPastCount[thisSourceVal][destVal][thisPastVal]++;
			sourcePastCount[thisSourceVal][thisPastVal]++;
			nextPastCount[destVal][thisPastVal]++;
			pastCount[thisPastVal]++;
			nextCount[destVal]++;
			// Now, update the combined embedding values and phases,
			//  for this phase we back out the oldest value which we'll no longer need:
			if (k > 0) {
				pastVal[destEmbeddingPhase] -= maxShiftedValue[states[r-1-(k-1)*destEmbeddingDelay][destRowIndex][destColumnIndex]];
				pastVal[destEmbeddingPhase] *= base; // and shift the others up
			}
			sourcePastVal[sourceEmbeddingPhase] -=
					maxShiftedSourceValue[
					    states[r-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay][sourceRowIndex][sourceColumnIndex]];
			sourcePastVal[sourceEmbeddingPhase] *= base; // and shift the others up
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
		}
	}

	/**
	 * 
	 * Returns the count of observations of the supplied past state
	 *  pastVal.
	 * The past state is indicated by a unique discrete integer representing the joint variable
	 *  of the k past states: (dest[n-k+1],dest[n-k+2],...,dest[n-1],dest[n]).
	 * The integer is computed as:<br/>
	 * pastVal = dest[n-k+1] * base^(k-1) + dest[n-k+2] * base^(k-2) + ... + dest[n-1] * base + dest[n]
	 * 
	 * 
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return count of observations of this given past state
	 */
	public int getPastCount(int pastVal) {
		return pastCount[pastVal];
	}
	
	/**
	 * Returns the probability of the supplied past state
	 *  pastVal.
	 * See {@link #getPastCount(int)} for how the joint value representing the past is calculated.
	 * 
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state
	 */
	public double getPastProbability(int pastVal) {
		return (double) pastCount[pastVal] / (double) observations;
	}

	/**
	 * Returns the count of observations of the past given past state and next value.
	 * 
	 * See {@link #getPastCount(int)} for how the joint value representing the past is calculated.
	 * 
	 * @param destVal next state of the destination dest[n+1]
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return count of observations of the given past state and next state
	 */
	public int getNextPastCount(int destVal, int pastVal) {
		return nextPastCount[destVal][pastVal];
	}
	
	/**
	 * Returns the probability of the past given past state and next value.
	 * 
	 * See {@link #getPastCount(int)} for how the joint value representing the past is calculated.
	 * 
	 * @param destVal next state of the destination dest[n+1]
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state and next state
	 */
	public double getNextPastProbability(int destVal, int pastVal) {
		return (double) nextPastCount[destVal][pastVal] / (double) observations;
	}

	/**
	 * Returns the count of observations of the past given state dest[n]^k and the source state source[n]^l.
	 * 
	 * See {@link #getPastCount(int)} for how the joint values representing the past states are calculated.
	 * 
	 * @param sourceVal int representing the joint state of the source source[n]^l
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return count of observations of the given past state and the source state
	 */
	public int getSourcePastCount(int sourceVal, int pastVal) {
		return sourcePastCount[sourceVal][pastVal];
	}
	
	/**
	 * Returns the probability of the past given state dest[n]^k and the source state source[n]^l.
	 * 
	 * See {@link #getPastCount(int)} for how the joint values representing the past states are calculated.
	 * 
	 * @param sourceVal int representing the joint state of the source source[n]^l
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state and the source state
	 */
	public double getSourcePastProbability(int sourceVal, int pastVal) {
		return (double) sourcePastCount[sourceVal][pastVal] / (double) observations;
	}

	/**
	 * Returns the count of observations of the past given state dest[n]^k,
	 *  the next state of the destination dest[n+1] and the source state source[n]^l.
	 *  
	 * See {@link #getPastCount(int)} for how the joint values representing the past states are calculated.
	 * 
	 * @param sourceVal int representing the joint state of the source source[n]^l
	 * @param nextVal next state of the destination dest[n+1]
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return count of observations of the given past state, next state of destination and the source state
	 */
	public int getSourceNextPastCount(int sourceVal, int destVal, int pastVal) {
		return sourceNextPastCount[sourceVal][destVal][pastVal];
	}
	
	/**
	 * Returns the probability of the past given state dest[n]^k,
	 *  the next state of the destination dest[n+1] and the source state source[n]^l.
	 *  
	 * See {@link #getPastCount(int)} for how the joint values representing the past states are calculated.
	 * 
	 * @param sourceVal int representing the joint state of the source source[n]^l
	 * @param nextVal next state of the destination dest[n+1]
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state, next state of destination and the source state
	 */
	public double getSourceNextPastProbability(int sourceVal, int destVal, int pastVal) {
		return (double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) observations;
	}

	/**
	 * Returns the count of observations of the next state dest[n+1].
	 * 
	 * @param nextVal next state of the destination dest[n+1]
	 * @return count of observations of the given next state
	 */
	public int getNextCount(int destVal) {
		return nextCount[destVal];
	}
	
	/**
	 * Returns the probability of the next state dest[n+1].
	 * 
	 * @param nextVal state of the next destination dest[n+1]
	 * @return probability of the given next state
	 */
	public double getNextProbability(int destVal) {
		return (double) nextCount[destVal] / (double) observations;
	}

	@Override
	public double computeAverageLocalOfObservations() {
		double te = 0.0;
		double teCont = 0.0;

		max = 0;
		min = 0;
		double meanSqLocals = 0;
		for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
			// compute p(past)
			// double p_past = (double) pastCount[pastVal] / (double) observations;
			if (pastCount[pastVal] == 0) {
				continue;
			}
			for (int destVal = 0; destVal < base; destVal++) {
				// compute p(dest,past)
				// double p_dest_past = (double) destPastCount[destVal][pastVal] / (double) observations;
				if (nextPastCount[destVal][pastVal] == 0) {
					continue;
				}
				double denom = (double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal];
				for (int sourceVal = 0; sourceVal < base_power_l; sourceVal++) {
					// Compute TE contribution:
					if (sourceNextPastCount[sourceVal][destVal][pastVal] != 0) {
						/* Double check: should never happen
						if ((sourcePastCount[sourceVal][pastVal] == 0) ||
							(destPastCount[destVal][pastVal] == 0) ||
							(pastCount[pastVal] == 0)) {
							throw new RuntimeException("one subcount was zero!!");
						}
						*/
						
						// compute p(source,dest,past)
						double p_source_dest_past = (double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) observations;
						// compute p(source,past)
						// double p_source_past = (double) sourcePastCount[sourceVal][pastVal] / (double) observations;
						double logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) sourcePastCount[sourceVal][pastVal]) /
						 	(denom);
						double localValue = Math.log(logTerm); // We'll / log_2 later, to save a floating pt op;
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
		
		te /= log_2;
		max /= log_2;
		min /= log_2;
		meanSqLocals /= (log_2 * log_2);
		
		average = te;
		std = Math.sqrt(meanSqLocals - average * average);
		estimateComputed = true;
		return te;
	}
	
	/**
	 * Returns the average active information storage from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @see ActiveInformationCalculatorDiscrete
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
					double localValue = Math.log(logTerm) / log_2;
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
				for (int sourceVal = 0; sourceVal < base_power_l; sourceVal++) {
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

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) {
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(observations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) {
		double actualTE = computeAverageLocalOfObservations();
		
		int numPermutationsToCheck = newOrderings.length;
		
		// Reconstruct the *joint* source values (not necessarily in order, but using joint values retains their l-tuples)
		int[] sourceValues = new int[observations];
		int t_s = 0;
		for (int sourceVal = 0; sourceVal < base_power_l; sourceVal++) {
			// Count up the number of times this joint source value was observed:
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
		
		// If we want a calculator just like this one, we should provide all of 
		//  the same parameters:
		TransferEntropyCalculatorDiscrete ate2 =
				new TransferEntropyCalculatorDiscrete(base, k, destEmbeddingDelay,
						sourceHistoryEmbedLength, sourceEmbeddingDelay, delay);
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
				// This looks like we're doing no source embedding --
				//  but we're already using joint source values in newSourceData
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
	
	@Override
	public AnalyticMeasurementDistribution computeSignificance()
			throws Exception {
		if (!estimateComputed) {
			computeAverageLocalOfObservations();
		}
		return new ChiSquareMeasurementDistribution(average,
				observations,
				(base_power_l - 1)*(base - 1)*(base_power_k));
	}

	/**
	 * Computes local transfer entropy for the given values
	 * 
	 * See {@link #getPastCount(int)} for how the joint values representing the past are calculated.
	 * 
	 * @param sourceCurrent int representing the joint state of the source source[n]^l
	 * @param destNext next state of the destination dest[n+1]
	 * @param destPast int representing the joint state of the past of the destination dest[n]^k
	 * 
	 * @return local TE for the given observation
	 */
	public double computeLocalFromPreviousObservations(int sourceCurrent, int destNext, int destPast){

		double logTerm = ((double) sourceNextPastCount[sourceCurrent][destNext][destPast] / (double) sourcePastCount[sourceCurrent][destPast]) /
	 		((double) nextPastCount[destNext][destPast] / (double) pastCount[destPast]);
		return Math.log(logTerm) / log_2;
	}

	/**
	 * Computes local apparent transfer entropy for the given
	 *  states, using PDFs built up from observations previously
	 *  sent in via the addObservations method.
	 *  
 	 * @param source source time-series
	 * @param dest destination time-series. 
	 *  Must be same length as source
	 * @return time-series of local TE values
	 */
	public double[] computeLocalFromPreviousObservations(int source[], int dest[]){
		int timeSteps = dest.length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localTE = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;

		if (timeSteps - startObservationTime <= 0) {
			// No observations to compute locals for
			return localTE;
		}
		
		// Initialise and store the current previous values;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[] pastVal = new int[destEmbeddingDelay];
		for (int d = 0; d < destEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			pastVal[d] = 0;
			for (int p = 0; p < k-1; p++) {
				pastVal[d] += dest[startObservationTime + d - 1
				                   - (k-1)*destEmbeddingDelay
				                   + p*destEmbeddingDelay];
				pastVal[d] *= base;
			}
		}
		// Next for the source:
		int[] sourcePastVal = new int[sourceEmbeddingDelay];
		for (int d = 0; d < sourceEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			sourcePastVal[d] = 0;
			for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
				sourcePastVal[d] += source[startObservationTime + d - delay
				                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
				                    + p*sourceEmbeddingDelay];
				sourcePastVal[d] *= base;
			}
		}
		
		// now compute the local values
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
		double logTerm;
		for (int t = startObservationTime; t < timeSteps; t++) {
			// First update the embedding values for the current
			//  phases of the embeddings:
			if (k > 0) {
				pastVal[destEmbeddingPhase] += dest[t-1];
			}
			sourcePastVal[sourceEmbeddingPhase] += source[t-delay];
			destVal = dest[t];
			int thisPastVal = pastVal[destEmbeddingPhase];
			int thisSourceVal = sourcePastVal[sourceEmbeddingPhase];
			// Now compute the local value
			logTerm = ((double) sourceNextPastCount[thisSourceVal][destVal][thisPastVal] /
							(double) sourcePastCount[thisSourceVal][thisPastVal]) /
			 		((double) nextPastCount[destVal][thisPastVal] / (double) pastCount[thisPastVal]);
			localTE[t] = Math.log(logTerm) / log_2;
			average += localTE[t];
			if (localTE[t] > max) {
				max = localTE[t];
			} else if (localTE[t] < min) {
				min = localTE[t];
			}
			// Now, update the combined embedding values and phases,
			//  for this phase we back out the oldest value which we'll no longer need:
			if (k > 0) {
				pastVal[destEmbeddingPhase] -= maxShiftedValue[dest[t-1-(k-1)*destEmbeddingDelay]];
				pastVal[destEmbeddingPhase] *= base; // and shift the others up
			}
			sourcePastVal[sourceEmbeddingPhase] -=
					maxShiftedSourceValue[
					    source[t-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay]];
			sourcePastVal[sourceEmbeddingPhase] *= base; // and shift the others up
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
		}

		average = average/(double) (timeSteps - startObservationTime);
		
		return localTE;
	}

	/**
	 * Computes local transfer for the given
	 *  multivariate states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by j column will
 	 *  have their local TE computed.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i: we
	 *  compute transfer is j cells to the right, using observations
	 *  across all column pairs separated by j) 
	 * @return multivariate time series of local TE values
	 *  (first index is time, second index is destination variable)
	 */
	public double[][] computeLocalFromPreviousObservations(int states[][], int j){
		int timeSteps = states.length;
		if ((timeSteps == 0) || (states[0] == null)) {
			// No variables supplied
			return new double[timeSteps][];
		}
		int variables = states[0].length;
		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localTE = new double[timeSteps][variables];
		average = 0;
		max = 0;
		min = 0;
		if (timeSteps - startObservationTime <= 0) {
			// No observations to add
			return localTE;
		}

		// Initialise and store the current previous values for each column;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[][] pastVal = new int[variables][destEmbeddingDelay];
		for (int c = 0; c < variables; c++) {
			for (int d = 0; d < destEmbeddingDelay; d++) {
				// Compute the current previous values for
				//  phase d of the embedding delay, but leave
				//  out the most recent value (we'll add those in
				//  in the main loop)
				pastVal[c][d] = 0;
				for (int p = 0; p < k-1; p++) {
					pastVal[c][d] += states[startObservationTime + d - 1
					                   - (k-1)*destEmbeddingDelay
					                   + p*destEmbeddingDelay][c];
					pastVal[c][d] *= base;
				}
			}
		}
		// Next for the source:
		int[][] sourcePastVal = new int[variables][sourceEmbeddingDelay];
		for (int c = 0; c < variables; c++) {
			int sourceVariable = c-j;
			if ((sourceVariable < 0) || (sourceVariable >= variables)) {
				// Source variable is out of bounds unless we are using periodic boundary conditions
				if (periodicBoundaryConditions) {
					sourceVariable = (sourceVariable+variables) % variables;
				} else {
					// Don't add this to our observations
					continue;
				}
			}
			for (int d = 0; d < sourceEmbeddingDelay; d++) {
				// Compute the current previous values for
				//  phase d of the embedding delay, but leave
				//  out the most recent value (we'll add those in
				//  in the main loop)
				sourcePastVal[c][d] = 0;
				for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
					sourcePastVal[c][d] += states[startObservationTime + d - delay
					                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
					                    + p*sourceEmbeddingDelay][sourceVariable];
					sourcePastVal[c][d] *= base;
				}
			}
		}
		
		// now compute the local values
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
		double logTerm;
		for (int t = startObservationTime; t < timeSteps; t++) {
			for (int c = 0; c < variables; c++) {
				int sourceVariable = c-j;
				if ((sourceVariable < 0) || (sourceVariable >= variables)) {
					// Source variable is out of bounds unless we are using periodic boundary conditions
					if (periodicBoundaryConditions) {
						sourceVariable = (sourceVariable+variables) % variables;
					} else {
						// Don't compute a local value for this one
						continue;
					}
				}
				// First update the embedding values for the current
				//  phases of the embeddings:
				if (k > 0) {
					pastVal[c][destEmbeddingPhase] += states[t-1][c];
				}
				sourcePastVal[c][sourceEmbeddingPhase] += states[t-delay][sourceVariable];
				destVal = states[t][c];
				int thisPastVal = pastVal[c][destEmbeddingPhase];
				int thisSourceVal = sourcePastVal[c][sourceEmbeddingPhase];
				// Now compute the local value
				logTerm = ((double) sourceNextPastCount[thisSourceVal][destVal][thisPastVal] /
							(double) sourcePastCount[thisSourceVal][thisPastVal]) /
			 		((double) nextPastCount[destVal][thisPastVal] / (double) pastCount[thisPastVal]);
				localTE[t][c] = Math.log(logTerm) / log_2;
				average += localTE[t][c];
				if (localTE[t][c] > max) {
					max = localTE[t][c];
				} else if (localTE[t][c] < min) {
					min = localTE[t][c];
				}
				// Now, update the combined embedding values and phases,
				//  for this phase we back out the oldest value which we'll no longer need:
				if (k > 0) {
					pastVal[c][destEmbeddingPhase] -= maxShiftedValue[states[t-1-(k-1)*destEmbeddingDelay][c]];
					pastVal[c][destEmbeddingPhase] *= base; // and shift the others up
				}
				sourcePastVal[c][sourceEmbeddingPhase] -=
						maxShiftedSourceValue[
						    states[t-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay][sourceVariable]];
				sourcePastVal[c][sourceEmbeddingPhase] *= base; // and shift the others up
			}
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
		}
		
		if (periodicBoundaryConditions) {
			average = average/(double) ((timeSteps - startObservationTime) * variables);
		} else {
			average = average/(double) ((timeSteps - startObservationTime) * (variables - Math.abs(j)));
		}
		
		return localTE;
	}
	
	/**
	 * Computes local transfer for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by h rows and j columns
 	 *  will have their local TE computed.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param h - number of rows to compute transfer entropy across
	 * 	(i.e. source is in row i-h, dest is column i) 
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i) 
	 * @return multivariate time series of local TE values
	 *  (first index is time, second index is destination variable
	 *  row number, third is destination variable column number)
	 */
	public double[][][] computeLocalFromPreviousObservations(int states[][][], int h, int j){
		
		int timeSteps = states.length;
		if ((timeSteps == 0) || (states[0] == null)) {
			// No variables supplied
			return new double[timeSteps][][];
		}
		int agentRows = states[0].length;
		if (agentRows == 0) {
			return new double[timeSteps][agentRows][];
		}
		int agentColumns = states[0][0].length;
		if (agentRows == 0) {
			return new double[timeSteps][agentRows][agentColumns];
		}
		if (timeSteps - startObservationTime <= 0) {
			// No observations to add
			return new double[timeSteps][][];
		}
		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][][] localTE = new double[timeSteps][agentRows][agentColumns];
		average = 0;
		max = 0;
		min = 0;		

		// Initialise and store the current previous values for each variable;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[][][] pastVal = new int[agentRows][agentColumns][destEmbeddingDelay];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				for (int d = 0; d < destEmbeddingDelay; d++) {
					// Compute the current previous values for
					//  phase d of the embedding delay, but leave
					//  out the most recent value (we'll add those in
					//  in the main loop)
					pastVal[r][c][d] = 0;
					for (int p = 0; p < k-1; p++) {
						pastVal[r][c][d] += states[startObservationTime + d - 1
						                   - (k-1)*destEmbeddingDelay
						                   + p*destEmbeddingDelay][r][c];
						pastVal[r][c][d] *= base;
					}
				}
			}
		}
		// Next for the source:
		int[][][] sourcePastVal = new int[agentRows][agentColumns][sourceEmbeddingDelay];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
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
				// Now initialise the embedding
				for (int d = 0; d < sourceEmbeddingDelay; d++) {
					// Compute the current previous values for
					//  phase d of the embedding delay, but leave
					//  out the most recent value (we'll add those in
					//  in the main loop)
					sourcePastVal[r][c][d] = 0;
					for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
						sourcePastVal[r][c][d] += states[startObservationTime + d - delay
						                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
						                    + p*sourceEmbeddingDelay][sourceAgentRow][sourceAgentColumn];
						sourcePastVal[r][c][d] *= base;
					}
				}
			}
		}
		
		// Now compute the local values
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
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
							// Don't compute a local value here
							continue;
						}
					}
					int sourceAgentColumn = c-j;
					if ((sourceAgentColumn < 0) || (sourceAgentColumn >= agentColumns)) {
						// Source agent is out of bounds unless we are using periodic boundary conditions
						if (periodicBoundaryConditions) {
							sourceAgentColumn = (sourceAgentColumn+agentColumns) % agentColumns;
						} else {
							// Don't compute a local value here
							continue;
						}
					}
					// First update the embedding values for the current
					//  phases of the embeddings:
					if (k > 0) {
						pastVal[r][c][destEmbeddingPhase] += states[t-1][r][c];
					}
					sourcePastVal[r][c][sourceEmbeddingPhase] += states[t-delay][sourceAgentRow][sourceAgentColumn];
					// Add to the count for this particular transition:
					// (cell's assigned as above)
					destVal = states[t][r][c];
					int thisPastVal = pastVal[r][c][destEmbeddingPhase];
					int thisSourceVal = sourcePastVal[r][c][sourceEmbeddingPhase];
					// Now compute the local value
					logTerm = ((double) sourceNextPastCount[thisSourceVal][destVal][thisPastVal] /
							    (double) sourcePastCount[thisSourceVal][thisPastVal]) /
				 		((double) nextPastCount[destVal][thisPastVal] / (double) pastCount[thisPastVal]);
					localTE[t][r][c] = Math.log(logTerm) / log_2;
					average += localTE[t][r][c];
					if (localTE[t][r][c] > max) {
						max = localTE[t][r][c];
					} else if (localTE[t][r][c] < min) {
						min = localTE[t][r][c];
					}
					// Now, update the combined embedding values and phases,
					//  for this phase we back out the oldest value which we'll no longer need:
					if (k > 0) {
						pastVal[r][c][destEmbeddingPhase] -= maxShiftedValue[states[t-1-(k-1)*destEmbeddingDelay][r][c]];
						pastVal[r][c][destEmbeddingPhase] *= base; // and shift the others up
					}
					sourcePastVal[r][c][sourceEmbeddingPhase] -=
							maxShiftedSourceValue[
							    states[t-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay][sourceAgentRow][sourceAgentColumn]];
					sourcePastVal[r][c][sourceEmbeddingPhase] *= base; // and shift the others up
				}
			}
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
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
	 *  single source-destination pair of the multi-agent system,
	 *  using pdfs built up from observations previously
	 *  sent in via the addObservations methods.
 	 * This call should be made as opposed to {@link #addObservations(int[][], int)}
 	 *  for computing local TE for heterogeneous agents.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param sourceIndex source variable index in states
	 * @param destIndex destination variable index in states
	 * @return time-series of local TE values between the series
	 */
	public double[] computeLocalFromPreviousObservations(int states[][], int sourceIndex, int destIndex){

		int timeSteps = states.length;
		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localTE = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;

		if (timeSteps - startObservationTime <= 0) {
			// No observations to add
			return localTE;
		}
		
		// Initialise and store the current previous values;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[] pastVal = new int[destEmbeddingDelay];
		for (int d = 0; d < destEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			pastVal[d] = 0;
			for (int p = 0; p < k-1; p++) {
				pastVal[d] += states[startObservationTime + d - 1
				                   - (k-1)*destEmbeddingDelay
				                   + p*destEmbeddingDelay][destIndex];
				pastVal[d] *= base;
			}
		}
		// Next for the source:
		int[] sourcePastVal = new int[sourceEmbeddingDelay];
		for (int d = 0; d < sourceEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			sourcePastVal[d] = 0;
			for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
				sourcePastVal[d] += states[startObservationTime + d - delay
				                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
				                    + p*sourceEmbeddingDelay][sourceIndex];
				sourcePastVal[d] *= base;
			}
		}
		
		// Now compute the local values
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
		double logTerm;
		for (int r = startObservationTime; r < timeSteps; r++) {
			// First update the embedding values for the current
			//  phases of the embeddings:
			if (k > 0) {
				pastVal[destEmbeddingPhase] += states[r-1][destIndex];
			}
			sourcePastVal[sourceEmbeddingPhase] += states[r-delay][sourceIndex];
			destVal = states[r][destIndex];
			int thisPastVal = pastVal[destEmbeddingPhase];
			int thisSourceVal = sourcePastVal[sourceEmbeddingPhase];
			// Now compute the local value
			logTerm = ((double) sourceNextPastCount[thisSourceVal][destVal][thisPastVal] /
						(double) sourcePastCount[thisSourceVal][thisPastVal]) /
		 		((double) nextPastCount[destVal][thisPastVal] / (double) pastCount[thisPastVal]);
			localTE[r] = Math.log(logTerm) / log_2;
			average += localTE[r];
			if (localTE[r] > max) {
				max = localTE[r];
			} else if (localTE[r] < min) {
				min = localTE[r];
			}
			// Now, update the combined embedding values and phases,
			//  for this phase we back out the oldest value which we'll no longer need:
			if (k > 0) {
				pastVal[destEmbeddingPhase] -= maxShiftedValue[states[r-1-(k-1)*destEmbeddingDelay][destIndex]];
				pastVal[destEmbeddingPhase] *= base; // and shift the others up
			}
			sourcePastVal[sourceEmbeddingPhase] -=
					maxShiftedSourceValue[
					    states[r-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay][sourceIndex]];
			sourcePastVal[sourceEmbeddingPhase] *= base; // and shift the others up
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
		}
		average = average/(double) (timeSteps - startObservationTime);
		
		return localTE;
	}

	/**
	 * Computes local transfer for the given
	 *  single source-destination pair of the multi-agent system,
	 *  using pdfs built up from observations previously
	 *  sent in via the addObservations method.
 	 * This call should be made as opposed to {@link #addObservations(int[][][], int, int)}
 	 *  for computing local TE for heterogeneous agents.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param sourceRowIndex source variable row index in states
	 * @param sourceColumnIndex source variable column index in states
	 * @param destRowIndex destination variable row index in states
	 * @param destColumnIndex destination variable column index in states
	 * @return time-series of local TE values between the series
	 */
	public double[] computeLocalFromPreviousObservations(int states[][][],
			int sourceRowIndex, int sourceColumnIndex, int destRowIndex, int destColumnIndex){
		int timeSteps = states.length;
		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localTE = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;
		if (timeSteps - startObservationTime <= 0) {
			// No observations to add
			return localTE;
		}
		
		// Initialise and store the current previous values;
		//  one for each phase of the embedding delay.
		// First for the destination:
		int[] pastVal = new int[destEmbeddingDelay];
		for (int d = 0; d < destEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			pastVal[d] = 0;
			for (int p = 0; p < k-1; p++) {
				pastVal[d] += states[startObservationTime + d - 1
				                   - (k-1)*destEmbeddingDelay
				                   + p*destEmbeddingDelay][destRowIndex][destColumnIndex];
				pastVal[d] *= base;
			}
		}
		// Next for the source:
		int[] sourcePastVal = new int[sourceEmbeddingDelay];
		for (int d = 0; d < sourceEmbeddingDelay; d++) {
			// Compute the current previous values for
			//  phase d of the embedding delay, but leave
			//  out the most recent value (we'll add those in
			//  in the main loop)
			sourcePastVal[d] = 0;
			for (int p = 0; p < sourceHistoryEmbedLength - 1; p++) {
				sourcePastVal[d] += states[startObservationTime + d - delay
				                    - (sourceHistoryEmbedLength-1)*sourceEmbeddingDelay
				                    + p*sourceEmbeddingDelay][sourceRowIndex][sourceColumnIndex];
				sourcePastVal[d] *= base;
			}
		}
		
		// Now compute the local values
		int destVal, destEmbeddingPhase = 0, sourceEmbeddingPhase = 0;
		double logTerm;
		for (int r = startObservationTime; r < timeSteps; r++) {
			// First update the embedding values for the current
			//  phases of the embeddings:
			if (k > 0) {
				pastVal[destEmbeddingPhase] += states[r-1][destRowIndex][destColumnIndex];
			}
			sourcePastVal[sourceEmbeddingPhase] += states[r-delay][sourceRowIndex][sourceColumnIndex];
			destVal = states[r][destRowIndex][destColumnIndex];
			int thisPastVal = pastVal[destEmbeddingPhase];
			int thisSourceVal = sourcePastVal[sourceEmbeddingPhase];
			// Now compute the local value
			logTerm = ((double) sourceNextPastCount[thisSourceVal][destVal][thisPastVal] /
						(double) sourcePastCount[thisSourceVal][thisPastVal]) /
		 		((double) nextPastCount[destVal][thisPastVal] / (double) pastCount[thisPastVal]);
			localTE[r] = Math.log(logTerm) / log_2;
			average += localTE[r];
			if (localTE[r] > max) {
				max = localTE[r];
			} else if (localTE[r] < min) {
				min = localTE[r];
			}
			// Now, update the combined embedding values and phases,
			//  for this phase we back out the oldest value which we'll no longer need:
			if (k > 0) {
				pastVal[destEmbeddingPhase] -= maxShiftedValue[states[r-1-(k-1)*destEmbeddingDelay][destRowIndex][destColumnIndex]];
				pastVal[destEmbeddingPhase] *= base; // and shift the others up
			}
			sourcePastVal[sourceEmbeddingPhase] -=
					maxShiftedSourceValue[
					    states[r-delay-(sourceHistoryEmbedLength-1)*sourceEmbeddingDelay][sourceRowIndex][sourceColumnIndex]];
			sourcePastVal[sourceEmbeddingPhase] *= base; // and shift the others up
			// then update the phase
			destEmbeddingPhase = (destEmbeddingPhase + 1) % destEmbeddingDelay; 
			sourceEmbeddingPhase = (sourceEmbeddingPhase + 1) % sourceEmbeddingDelay; 
		}

		average = average/(double) (timeSteps - startObservationTime);
		
		return localTE;
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy between two time series
	 * Return a time series of local values.
	 * First max(k,l) values are zeros since TE is not defined there
	 * 
 	 * @param sourceStates source time-series
	 * @param destStates destination time-series. 
	 *  Must be same length as sourceStates
	 * @return time-series of local TE values
	 */
	public double[] computeLocal(int sourceStates[], int destStates[]) {
		
		initialise();
		addObservations(sourceStates, destStates);
		return computeLocalFromPreviousObservations(sourceStates, destStates);
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents.
	 * Return a 2D spatiotemporal array of local values.
	 * First max(k,l) values are zeros since TE is not defined there.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by j column will
 	 *  have their local TE computed.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i: we
	 *  compute transfer is j cells to the right, using observations
	 *  across all column pairs separated by j) 
	 * @return multivariate time series of local TE values
	 *  (first index is time, second index is destination variable)
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
	 * First max(k,l) values are zeros since TE is not defined there.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by h rows and j columns
 	 *  will have their local TE computed.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param h - number of rows to compute transfer entropy across
	 * 	(i.e. source is in row i-h, dest is column i) 
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i) 
	 * @return multivariate time series of local TE values
	 *  (first index is time, second index is destination variable
	 *  row number, third is destination variable column number)
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
	 * Return the average TE.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by j column will
 	 *  have their local TE computed.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i: we
	 *  compute transfer is j cells to the right, using observations
	 *  across all column pairs separated by j) 
	 * @return average TE across j variables to the right
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
	 * Return the average.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by h rows and j columns
 	 *  will have their PDFs combined.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param h - number of rows to compute transfer entropy across
	 * 	(i.e. source is in row i-h, dest is column i) 
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i) 
	 * @return
	 */
	public double computeAverageLocal(int states[][][], int h, int j) {
		
		initialise();
		addObservations(states, h, j);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy between specific variables in
	 * a 2D spatiotemporal multivariate time-series.
	 * First max(k,l) values are zeros since TE is not defined there.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param sourceCol source variable index in states
	 * @param destCol destination variable index in states
	 * @return time-series of local TE values between the series
	 */
	public double[] computeLocal(int states[][], int sourceCol, int destCol) {
		
		initialise();
		addObservations(states, sourceCol, destCol);
		return computeLocalFromPreviousObservations(states, sourceCol, destCol);
	}

	/**
	 * Standalone routine to 
	 * computes local transfer for the given
	 *  single source-destination pair of the 3D multi-agent system.
	 * This method suitable for heterogeneous variables.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param sourceRowIndex source variable row index in states
	 * @param sourceColumnIndex source variable column index in states
	 * @param destRowIndex destination variable row index in states
	 * @param destColumnIndex destination variable column index in states
	 * @return time-series of local TE values between the series
	 */
	public double[] computeLocal(int states[][][], int sourceRowIndex, int sourceColumnIndex,
											int destRowIndex, int destColumnIndex) {
		
		initialise();
		addObservations(states, sourceRowIndex, sourceColumnIndex, destRowIndex, destColumnIndex);
		return computeLocalFromPreviousObservations(states, sourceRowIndex, sourceColumnIndex,
				destRowIndex, destColumnIndex);
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy between specific variables in
	 * a 2D spatiotemporal multivariate time-series.
	 * Returns the average.
	 * This method suitable for heterogeneous agents.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param sourceCol source variable index in states
	 * @param destCol destination variable index in states
	 * @return average TE for the given pair
	 */
	public double computeAverageLocal(int states[][], int sourceCol, int destCol) {
		
		initialise();
		addObservations(states, sourceCol, destCol);
		return computeAverageLocalOfObservations();
	}
	
	/**
	 * Standalone routine to 
	 * compute local transfer entropy between specific variables in
	 * a 3D spatiotemporal multivariate time-series.
	 * Returns the average.
	 * This method suitable for heterogeneous agents.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param sourceRowIndex source variable row index in states
	 * @param sourceColumnIndex source variable column index in states
	 * @param destRowIndex destination variable row index in states
	 * @param destColumnIndex destination variable column index in states
	 * @return average TE for the given pair
	 */
	public double computeAverageLocal(int states[][][], int sourceRowIndex, int sourceColumnIndex,
			int destRowIndex, int destColumnIndex) {
		
		initialise();
		addObservations(states, sourceRowIndex, sourceColumnIndex, destRowIndex, destColumnIndex);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Whether we assume periodic boundary conditions in the calls
	 *  for homogeneous variables.
	 *  
	 * @return as above
	 */
	public boolean isPeriodicBoundaryConditions() {
		return periodicBoundaryConditions;
	}
	/**
	 * set whether we assume periodic boundary conditions in the calls
	 *  for homogeneous variables.
	 *  
	 * @param periodicBoundaryConditions as above
	 */
	public void setPeriodicBoundaryConditions(boolean periodicBoundaryConditions) {
		this.periodicBoundaryConditions = periodicBoundaryConditions;
	}
}
