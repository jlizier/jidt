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

package infodynamics.measures.continuous.kernel;

import infodynamics.measures.continuous.TransferEntropyCalculator;
import infodynamics.measures.continuous.TransferEntropyCommon;
import infodynamics.measures.continuous.kernel.TransferEntropyKernelCounts;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

import java.util.Iterator;


/**
 * 
 * <p>
 * Implements a transfer entropy calculator using kernel estimation.
 * (see Schreiber, PRL 85 (2) pp.461-464, 2000)</p> 
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct</li>
 * 		<li>SetProperty() for each property</li>
 *		<li>intialise()</li>
 * 		<li>setObservations(), or [startAddObservations(), addObservations()*, finaliseAddObservations()]
 *   Note: If not using setObservations(), the results from computeLocal or getSignificance
 *    are not likely to be particularly sensible.</li> 
 * 		<li>computeAverageLocalOfObservations() or ComputeLocalOfPreviousObservations()</li>
 * 	</ol>
 * </p>
 * 
 * <p>
 * TODO Implement dynamic correlation exclusion with multiple observation sets. (see the
 *  way this is done in Plain calculator).
 * </p>
 * 
 * @author Joseph Lizier
 * @see For transfer entropy: Schreiber, PRL 85 (2) pp.461-464, 2000; http://dx.doi.org/10.1103/PhysRevLett.85.461
 * @see For local transfer entropy: Lizier et al, PRE 77, 026110, 2008; http://dx.doi.org/10.1103/PhysRevE.77.026110
 *
 */
public class TransferEntropyCalculatorKernel
	extends TransferEntropyCommon implements TransferEntropyCalculator {

	protected KernelEstimatorTransferEntropy teKernelEstimator = null;

	// Keep joint vectors so we don't need to regenerate them
	protected double[][] destPastVectors;
	protected double[]   destNextValues;
	protected double[]   sourceValues;
	
	private boolean normalise = true;
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	private boolean dynCorrExcl = false;
	private int dynCorrExclTime = 100;
	public static final String DYN_CORR_EXCL_TIME_NAME = "DYN_CORR_EXCL";
	
	private boolean forceCompareToAll = false;
	public static final String FORCE_KERNEL_COMPARE_TO_ALL = "FORCE_KERNEL_COMPARE_TO_ALL";
	
	/**
	 * Default value for epsilon (kernel width)
	 */
	public static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	private double epsilon = DEFAULT_EPSILON;
	public static final String EPSILON_PROP_NAME = "EPSILON";
	
	/**
	 * Creates a new instance of the kernel-estimate style transfer entropy calculator
	 *
	 */
	public TransferEntropyCalculatorKernel() {
		super();
		teKernelEstimator = new KernelEstimatorTransferEntropy();
		teKernelEstimator.setNormalise(normalise);
	}

	/**
	 * Initialises the calculator with the existing value for epsilon
	 * 
	 * @param k history length
	 */
	public void initialise(int k) throws Exception {
		initialise(k, epsilon);
	}
	
	/**
	 * Initialises the calculator
	 * 
	 * @param k history length
	 * @param epsilon kernel width
	 */
	public void initialise(int k, double epsilon) throws Exception {
		this.epsilon = epsilon;
		super.initialise(k); // calls initialise();
	}

	/**
	 * Initialise using default or existing values for k and epsilon
	 */
	public void initialise() {
		teKernelEstimator.initialise(k, epsilon);
		destPastVectors = null;
		destNextValues = null;
		sourceValues = null;
	}
	
	/**
	 * Set properties for the transfer entropy calculator.
	 * These can include:
	 * <ul>
	 * 		<li>K_PROP_NAME</li>
	 * 		<li>EPSILON_PROP_NAME</li>
	 * 		<li>NORMALISE_PROP_NAME</li>
	 * 		<li>DYN_CORR_EXCL_TIME_NAME</li>
	 * 		<li>FORCE_KERNEL_COMPARE_TO_ALL</li>
	 * </ul> 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception 
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		super.setProperty(propertyName, propertyValue);
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			epsilon = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
			teKernelEstimator.setNormalise(normalise);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
			if (dynCorrExcl) {
				teKernelEstimator.setDynamicCorrelationExclusion(dynCorrExclTime);
			} else {
				teKernelEstimator.clearDynamicCorrelationExclusion();
			}
		} else if (propertyName.equalsIgnoreCase(FORCE_KERNEL_COMPARE_TO_ALL)) {
			forceCompareToAll = Boolean.parseBoolean(propertyValue);
			teKernelEstimator.setForceCompareToAll(forceCompareToAll);
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/**
	 * Flag that the observations are complete, probability distribution functions can now be built.
	 *
	 */
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[] destination : vectorOfDestinationObservations) {
			totalObservations += destination.length - k;
		}
		destPastVectors = new double[totalObservations][k];
		destNextValues = new double[totalObservations];
		sourceValues = new double[totalObservations];
		
		// Construct the joint vectors from the given observations
		int startObservation = 0;
		Iterator<double[]> iterator = vectorOfDestinationObservations.iterator();
		for (double[] source : vectorOfSourceObservations) {
			double[] destination = iterator.next();
			double[][] currentDestPastVectors = makeJointVectorForPast(destination);
			MatrixUtils.arrayCopy(currentDestPastVectors, 0, 0,
					destPastVectors, startObservation, 0, currentDestPastVectors.length, k);
			System.arraycopy(destination, k, destNextValues, startObservation, destination.length - k);
			System.arraycopy(source, k - 1, sourceValues, startObservation, source.length - k);
			startObservation += destination.length - k;
		}
		
		// Now set the joint vectors in the kernel estimators
		teKernelEstimator.setObservations(destPastVectors, destNextValues, sourceValues);

		// Store whether there was more than one observation set:
		addedMoreThanOneObservationSet = vectorOfDestinationObservations.size() > 1;
		if (addedMoreThanOneObservationSet && dynCorrExcl) {
			// We have not properly implemented dynamic correlation exclusion for
			//  multiple observation sets, so throw an error
			throw new RuntimeException("Addition of multiple observation sets is not currently " +
					"supported with property DYN_CORR_EXCL set");
		}

		// And clear the vector of observations
		vectorOfSourceObservations = null;
		vectorOfDestinationObservations = null;
	}
	
	/**
	 * <p>Computes the average Transfer Entropy for the previously supplied observations</p> 
	 * 
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		
		double te = 0.0;
		if (debug) {
			MatrixUtils.printMatrix(System.out, destPastVectors);
		}
		for (int b = 0; b < totalObservations; b++) {
			TransferEntropyKernelCounts kernelCounts = teKernelEstimator.getCount(destPastVectors[b],
					destNextValues[b], sourceValues[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (kernelCounts.countNextPastSource > 0) {
				logTerm = ((double) kernelCounts.countNextPastSource / (double) kernelCounts.countPastSource) /
							((double) kernelCounts.countNextPast / (double) kernelCounts.countPast);
				cont = Math.log(logTerm);
			}
			te += cont;
			if (debug) {
				System.out.println(b + ": " + destPastVectors[b][0] + " (" +
						kernelCounts.countNextPastSource + " / " + kernelCounts.countPastSource + ") / (" +
						kernelCounts.countNextPast + " / " + kernelCounts.countPast + ") = " + 
						logTerm + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
		}
		lastAverage = te / (double) totalObservations / Math.log(2.0);
		return lastAverage;
	}

	/**
	 * <p>Computes the average Transfer Entropy for the previously supplied observations,
	 *  using the Grassberger correction for the point count k: log_e(k) ~= digamma(k).</p>
	 * <p>Kaiser and Schreiber, Physica D 166 (2002) pp. 43-62 suggest (on p. 57) that for the TE
	 *   though the adverse correction of the bias correction is worse than the correction
	 *   itself (because the probabilities being multiplied/divided are not independent), 
	 *   so recommend not to use this method.
	 *   (Bias correction is implemented properly in the Kraskov et al. 
	 *   estimators, see {@link infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov}
	 * </p>
	 * <p>It is implemented here for testing purposes only.</p>
	 * 
	 */
	public double computeAverageLocalOfObservationsWithCorrection() throws Exception {
		
		double te = 0.0;
		int numNoNeighbours = 0;
		double contributionsNoNeighbours = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			TransferEntropyKernelCounts kernelCounts = teKernelEstimator.getCount(destPastVectors[b],
					destNextValues[b], sourceValues[b], b);
			double cont = 0.0;
			// Original code:
			/* if (kernelCounts.countNextPastSource > 0) {
				cont = MathsUtils.digamma(kernelCounts.countNextPastSource) - 
						MathsUtils.digamma(kernelCounts.countPastSource) -
						MathsUtils.digamma(kernelCounts.countNextPast) +
						MathsUtils.digamma(kernelCounts.countPast);
			} */
			// But Schreiber confirmed to me that with dynamic correlation 
			//  exclusion, when you may have no nearest neighbours,
			//  he was allowing contributions from other groups (e.g. countPastSource)
			//  to be added even if the full joint count was zero.
			// Implement it like this:
			//  TODO Do we need to correct the totalObservations
			//   divisor for each digamma sum to account for this?
			//   (Schreiber does not do that; for the moment we'll accept that
			//    as the right approach)
			if (kernelCounts.countPastSource > 0) {
				cont -= MathsUtils.digamma(kernelCounts.countPastSource);
			}
			if (kernelCounts.countNextPast > 0) {
				cont -= MathsUtils.digamma(kernelCounts.countNextPast);
			}
			if (kernelCounts.countPast > 0) {
				cont += MathsUtils.digamma(kernelCounts.countPast);
			}
			if (kernelCounts.countNextPastSource > 0) {
				cont += MathsUtils.digamma(kernelCounts.countNextPastSource);
			} else {
				// These contributions are from a set with no neighbours in 
				//  the full joint space.
				contributionsNoNeighbours += cont;
				numNoNeighbours++;
			}
			te += cont;
			/*
			if (debug) {
				System.out.println(b + ": " + cont + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
			*/
		}
		// Average it, and convert results to bytes
		lastAverage = te / (double) totalObservations / Math.log(2.0);
		if (debug) {
			System.out.printf("TE=%.4f, with %d contributions from 0 neighbour sets being %.4f\n",
					lastAverage,
					numNoNeighbours,
					contributionsNoNeighbours / (double) totalObservations / Math.log(2.0));
		}
		return lastAverage;
	}

	/**
	 * Computes the local transfer entropies for the previous supplied observations.
	 * 
	 * Where more than one time series has been added, the array
	 *  contains the local values for each tuple in the order in
	 *  which they were added.
	 * 
	 * If there was only a single time series added, the array
	 *  contains k zero values before the local values.
	 *  (This means the length of the return array is the same
	 *  as the length of the input time series).
	 * 
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		return computeLocalUsingPreviousObservations(null, null, true);
	}

	/**
	 * Comptues local transfer entropies for the given observations, using the previously supplied
	 * observations to compute the PDFs.
	 * I don't think it's such a good idea to do this for continuous variables (e.g. where
	 * one can get kernel estimates for probabilities of zero now) but I've implemented
	 * it anyway. I guess getting kernel estimates of zero here is no different than what
	 * can occur with dynamic correlation exclusion.
	 * 
	 * @param source
	 * @param destination
	 * @return
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double[] source, double[] destination) throws Exception {
		return computeLocalUsingPreviousObservations(source, destination, false);
	}
	
	/**
	 * Returns the local TE at every time point.
	 * 
	 * @param source
	 * @param destination
	 * @param isPreviousObservations
	 * @return
	 * @throws Exception
	 */
	private double[] computeLocalUsingPreviousObservations(double[] source, double[] destination, boolean isPreviousObservations) throws Exception {
		
		double[][] newDestPastVectors;
		double[] newDestNextValues;
		double[] newSourceValues;

		if (isPreviousObservations) {
			// We've already computed the joint vectors for these observations
			newDestPastVectors = destPastVectors;
			newDestNextValues = destNextValues;
			newSourceValues = sourceValues;
		} else {
			// We need to compute a new set of joint vectors
			newDestPastVectors = makeJointVectorForPast(destination);
			newDestNextValues = MatrixUtils.select(destination, k, destination.length - k);
			newSourceValues = MatrixUtils.select(source, k - 1, source.length - k);
		}

		double te = 0.0;
		int numLocalObservations = newDestPastVectors.length;
		double[] localTE;
		int offset = 0;
		if (isPreviousObservations && addedMoreThanOneObservationSet) {
			// We're returning the local values for a set of disjoint
			//  observations. So we don't add k zeros to the start
			localTE = new double[numLocalObservations];
			offset = 0;
		} else {
			localTE = new double[numLocalObservations + k];
			offset = k;
		}
		TransferEntropyKernelCounts kernelCounts;
		for (int b = 0; b < numLocalObservations; b++) {
			if (isPreviousObservations) {
				kernelCounts = teKernelEstimator.getCount(newDestPastVectors[b],
						newDestNextValues[b], newSourceValues[b], b);
			} else {
				kernelCounts = teKernelEstimator.getCount(newDestPastVectors[b],
						newDestNextValues[b], newSourceValues[b], -1);
			}
			double logTerm = 0.0;
			double local = 0.0;
			if (kernelCounts.countNextPastSource > 0) {
				logTerm = ((double) kernelCounts.countNextPastSource / (double) kernelCounts.countPastSource) /
							((double) kernelCounts.countNextPast / (double) kernelCounts.countPast);
				local = Math.log(logTerm);
			}
			localTE[offset + b] = local;
			te += local;
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (local/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
		}
		lastAverage = te / (double) numLocalObservations / Math.log(2.0);
		return localTE;
	}

	/**
	 * Compute the significance of obtaining the given average TE from the given observations
	 * 
	 * 	This is as per Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128.
	 *
	 * Basically, we shuffle the source observations against the destination tuples.
	 * This keeps the marginal PDFs the same (including the entropy rate of the destination)
	 *  but destroys any correlation between the source and state change of the destination.
	 * 
	 * @param numPermutationsToCheck number of new orderings of the source values to compare against
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(totalObservations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * As per {@link computeSignificance(int) computeSignificance()} but supplies
	 *  the re-orderings of the observations of the source variables.
	 *  
	 * 
	 * @param newOrderings first index is permutation number, i.e. newOrderings[i]
	 * 		is an array of 1 permutation of 0..n-1, where there were n observations.
	 * @return
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		
		double actualTE = computeAverageLocalOfObservations();
		
		// Space for the source observations:
		double[] oldSourceValues = sourceValues;
		
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the source in the destPastSourceVectors 
			//  and destNextPastSourceVectors vectors
			sourceValues = MatrixUtils.extractSelectedTimePoints(oldSourceValues, newOrderings[p]);
			
			// Make the equivalent operations of intialise
			teKernelEstimator.initialise(k, epsilon);
			// Make the equivalent operations of setObservations:
			teKernelEstimator.setObservations(destPastVectors, destNextValues, sourceValues);
			// And get a TE value for this realisation:
			double newTe = computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newTe;
			if (newTe >= actualTE) {
				countWhereTeIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the local variables:
		lastAverage = actualTE;
		sourceValues = oldSourceValues;
		// And set the kernel estimator back to their previous state
		teKernelEstimator.initialise(k, epsilon);
		teKernelEstimator.setObservations(destPastVectors, destNextValues, sourceValues);
		
		// And return the significance
		measDistribution.pValue = (double) countWhereTeIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualTE;
		return measDistribution;
	}

	public void setDebug(boolean debug) {
		super.setDebug(debug);
		teKernelEstimator.setDebug(debug);
	}
}
