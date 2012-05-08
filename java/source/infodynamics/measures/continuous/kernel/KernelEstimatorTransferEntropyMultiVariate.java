/**
 * 
 */
package infodynamics.measures.continuous.kernel;

import infodynamics.utils.MatrixUtils;

/**
 * <p>Kernel estimator for use with the transfer entropy.</p>
 * 
 * <p>Extends KernelEstimatorMultiVariate, using the super class to manage the history of the destination
 * variable, and adds the next state and source on top of this. Any calls to the super class methods will only
 * function on the joint history.
 * </p> 
 * 
 * @author Joseph Lizier joseph.lizier at gmail.com
 *
 */
public class KernelEstimatorTransferEntropyMultiVariate extends KernelEstimatorMultiVariate {

	private int sourceDimensions;
	private double[] epsilonSource;
	private double[] epsilonSourceInUse;

	// Keep a separate epsilon for the destination next state, just in case
	//  we're doing something funky and using different variables to track
	//  the destination's past and next state (e.g. in swarm analysis).
	private double epsilonDestNextFixed;
	private double[] epsilonDestNext;
	private double[] epsilonDestNextInUse;

	private double[][] destNext;
	private double[][] source;
	
	// Store the current observations passed in the synchronized method getCount
	//  waiting for callbacks from the underlying kernel estimator
	private double[] destNextObs;
	private double[] sourceObs;
	// Counts of destPastNext, destPastSource, destPastNextSource to be filled in 
	//  with the callbacks
	private int countNextPast;
	private int countPastSource;
	private int countNextPastSource;
	
	public KernelEstimatorTransferEntropyMultiVariate() {
		super();
		// Make sure when get a callbacl when correlated points are found
		makeCorrelatedPointAddedCallback = true;
	}

	public void initialise(int dimensions, double epsilon) {
		initialise(dimensions, dimensions, epsilon, epsilon);
	}
	
	public void initialise(int destDimensionsWithPast, int sourceDimensions, 
			double epsilonDest, double epsilonSource) {
		super.initialise(destDimensionsWithPast, epsilonDest);
		this.sourceDimensions = sourceDimensions;
		this.epsilonSource = new double[sourceDimensions];
		for (int d = 0; d < sourceDimensions; d++) {
			this.epsilonSource[d] = epsilonSource;
		}
		epsilonSourceInUse = new double[sourceDimensions];
		// Track that we're using a fixed epsilon for destination next
		epsilonDestNext = null;
		epsilonDestNextFixed = epsilonDest;
	}
	
	public void initialise(double[] epsilonDest, 
			double[] epsilonSource) {
		super.initialise(epsilonDest);
		this.epsilonSource = epsilonSource;
		epsilonSourceInUse = new double[sourceDimensions];
		// Assume that we are doing an ordinary TE computation (dest past
		//  and next state are the same variable)
		epsilonDestNext = epsilonDest;
	}

	/**
	 * 
	 * @param destPastVectors is the joint vector of the 
	 *   past k states of the <dimensions> joint variables.
	 *   Each vector has values x_c,t: [x_0,0; x_0,0; .. ; x_d,0; x_0,1 .. ] 
	 *   i.e. it puts all values for a given time step in at once.
	 * @param destNext
	 * @param source
	 */
	public void setObservations(double[][] destPastVectors,
			double[][] destNext, double[][] source) {
		setObservations(destPastVectors);
		// epsilonInUse has been computed for the destination.
		// TODO We could compute and set it directly here so
		//  we don't have a mismatch between any of the vector
		//  variables.
		
		if (normalise) {
			for (int d = 0; d < sourceDimensions; d++) {
				double std = MatrixUtils.stdDev(source, d);
				epsilonSourceInUse[d] = epsilonSource[d] * std;
			}
		} else {
			for (int d = 0; d < sourceDimensions; d++) {
				epsilonSourceInUse[d] = epsilonSource[d];
			}
		}

		// Set epsilon for the destination next state here, now
		//  that we can be sure of it's number of dimensions (in case
		//  we're doing something funky with different dest past and next
		//  state variables).
		int destNextDimensions = destNext[0].length;
		if (epsilonDestNext == null) {
			// We must be using a fixed epsilon
			epsilonDestNext = new double[destNextDimensions];
			for (int d = 0; d < destNextDimensions; d++) {
				epsilonDestNext[d] = epsilonDestNextFixed;
			}
		}
		epsilonDestNextInUse = new double[destNextDimensions];
		if (normalise) {
			for (int d = 0; d < destNextDimensions; d++) {
				double std = MatrixUtils.stdDev(destNext, d);
				epsilonDestNextInUse[d] = epsilonDestNext[d] * std;
			}
		} else {
			for (int d = 0; d < destNextDimensions; d++) {
				epsilonDestNextInUse[d] = epsilonDestNext[d];
			}
		}

		this.source = source;
		this.destNext = destNext;
	}
	
	/**
	 * Compute the required counts for Transfer Entropy using kernel estimation on the destination's past.
	 * Use callbacks to check if the joint counts need to be incremented.
	 * 
	 * If observationTimeStep &lt; 0, then no dynamic correlation exclusion will be attempted
	 * 
	 * @param destPast
	 * @param destNextObs
	 * @param sourceObs
	 * @param observationTimeStep
	 * @return
	 */
	public synchronized TransferEntropyKernelCounts getCount(
			double[] destPast, double[] destNextObs,
			double[] sourceObs, int observationTimeStep) {
		
		// Prepare for any callbacks
		countNextPast = 0;
		countPastSource = 0;
		countNextPastSource = 0;
		this.destNextObs = destNextObs;
		this.sourceObs = sourceObs;

		// Get the count, and have the joint counts filled in via callbacks
		int countPast;
		if (observationTimeStep < 0) {
			countPast = super.getCount(destPast);
		} else {
			countPast = super.getCount(destPast, observationTimeStep);
		}
		
		TransferEntropyKernelCounts teKernelCount = 
			new TransferEntropyKernelCounts(countPast, countNextPast, countPastSource, countNextPastSource);
		return teKernelCount;
	}

	public void setEpsSource(double[] epsilonSource) {
		this.epsilonSource = epsilonSource;
	}

	/**
	 * A callback for where a correlated point is found at
	 * correlatedTimeStep in the destination's past.
	 * Now check whether we need to increment the joint counts.
	 *
	 */
	protected void correlatedPointAddedCallback(int correlatedTimeStep) {
		boolean sourceMatches = false;
		if (stepKernel(sourceObs, source[correlatedTimeStep], epsilonSourceInUse) > 0) {
			countPastSource++;
			sourceMatches = true;
		}
		// The epsilons across the history of each of the joint
		//  destination variables should all be approximately 
		//  equal, so just use the first ones for each dimension.
		//  (i.e. use the first dimensions epsilons, which
		//  is achieved simply by passing in epsilon, since
		//  stepKernel just uses the first dimensions elements
		//  of the widths argument).
		if (stepKernel(destNextObs, destNext[correlatedTimeStep], epsilonDestNextInUse) > 0) {
			countNextPast++;
			if (sourceMatches) {
				countNextPastSource++;
			}
		}
	}

	/**
	 * A callback for where a correlated point is removed due to dynamic correlated exclusion.
	 * The removal is for the point at correlatedTimeStep in the destination's past.
	 * Now check whether we need to decrement the joint counts.
	 */
	protected void correlatedPointRemovedCallback(int removedCorrelatedTimeStep) {
		boolean sourceMatches = false;
		if (stepKernel(sourceObs, source[removedCorrelatedTimeStep], epsilonSourceInUse) > 0) {
			countPastSource--;
			sourceMatches = true;
		}
		// The epsilons across the history of each of the joint
		//  destination variables should all be approximately 
		//  equal, so just use the first ones for each dimension.
		//  (i.e. use the first dimensions epsilons, which
		//  is achieved simply by passing in epsilon, since
		//  stepKernel just uses the first dimensions elements
		//  of the widths argument).
		if (stepKernel(destNextObs, destNext[removedCorrelatedTimeStep], epsilonDestNextInUse) > 0) {
			countNextPast--;
			if (sourceMatches) {
				countNextPastSource--;
			}
		}
	}
	
	protected int stepKernel(double[] vector1, double[] vector2, double[] widths) {
		for (int d = 0; d < vector1.length; d++) {
			if (Math.abs(vector1[d] - vector2[d]) > widths[d]) {
				return 0;
			}
		}
		// The candidate is within epsilon of the observation
		return 1;
	}
}
