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

package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.MultiInfoCalculator;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

import java.util.Vector;


/**
 * <p>Compute the Multi-Information (or integration) using the Kraskov estimation method.
 * Two child classes actually implement the two algorithms in the Kraskov paper.</p>
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * @see http://dx.doi.org/10.1103/PhysRevE.69.066138
 * 
 * @author Joseph Lizier
 */
public abstract class MultiInfoCalculatorKraskov implements
		MultiInfoCalculator {

	/**
	 * we compute distances to the kth neighbour
	 */
	protected int k;
	protected double[][] data;
	protected boolean debug;
	protected double mi;
	protected boolean miComputed;
	private Vector<double[]> individualObservations;

	protected int N; // number of observations
	protected int V; // number of variables
	
	protected EuclideanUtils normCalculator;
	// Storage for the norms for each marginal variable from each observation to each other one
	protected double[][][] norms;
	// Keep the norms each time (making reordering very quick)
	//  (Should only be set to false for testing)
	protected boolean tryKeepAllPairsNorms = true;
	public static int MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM = 4000;
	
	public final static String PROP_K = "k";
	public final static String PROP_NORM_TYPE = "NORM_TYPE";
	public final static String PROP_TRY_TO_KEEP_ALL_PAIRS_NORM = "TRY_KEEP_ALL_PAIRS_NORM";
	
	public MultiInfoCalculatorKraskov() {
		super();
		k = 1; // by default
		normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
	}

	public void initialise(int dimensions) {
		V = dimensions;
		mi = 0.0;
		miComputed = false;
		norms = null;
		data = null;
	}

	/**
	 * Sets properties for the calculator.
	 * Valid properties include:
	 * <ul>
	 *  <li>{@link #PROP_K} - number of neighbouring points in joint kernel space</li>
	 * 	<li>{@link #PROP_NORM_TYPE}</li> - normalization type to apply to 
	 * 		working out the norms between the points in each marginal space.
	 * 		Options are defined by {@link EuclideanUtils#setNormToUse(String)} -
	 * 		default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 *  <li>{@link #PROP_TRY_TO_KEEP_ALL_PAIRS_NORM})</li>
	 * </ul>
	 * 
	 * @param propertyName
	 * @param propertyValue
	 */
	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normCalculator.setNormToUse(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_TRY_TO_KEEP_ALL_PAIRS_NORM)) {
			tryKeepAllPairsNorms = Boolean.parseBoolean(propertyValue);
		}
	}

	public void setObservations(double[][] observations) throws Exception {
		if ((observations == null) || (observations[0].length == 0)) {
			throw new Exception("Computing MI with a null set of data");
		}
		if (observations[0].length != V) {
			throw new Exception("Incorrect number of dimensions " + observations[0].length +
					" in supplied observations (expected " + V + ")");
		}
		data = observations;
		N = data.length;
	}

	/**
	 * Set observations from two separate time series (join the rows at each time step
	 *  together to make a joint vector)
	 * 
	 * @param observations1
	 * @param observations2
	 */
	public void setObservations(double[][] observations1, double[][] observations2) throws Exception {
		if ((observations1 == null) || (observations1[0].length == 0) ||
				(observations2 == null) || (observations2[0].length == 0)) {
			throw new Exception("Computing MI with a null set of data");
		}
		if (observations1.length != observations2.length) {
			throw new Exception("Length of the time series to be joined to not match");
		}
		if (observations1[0].length + observations2[0].length != V) {
			throw new Exception("Incorrect number of dimensions " + 
					(observations1[0].length + observations2[0].length) +
					" in supplied observations (expected " + V + ")");
		}
		N = observations1.length;
		data = new double[N][V];
		for (int t = 0; t < N; t++) {
			int v = 0;
			for (int i = 0; i < observations1[t].length; i++) {
				data[t][v++] = observations1[t][i];
			}
			for (int i = 0; i < observations2[t].length; i++) {
				data[t][v++] = observations2[t][i];
			}
		}
		return;
	}
	
	/**
	 * User elects to set observations one by one rather than in one go.
	 * Will need to call endIndividualObservations before calling any of the
	 *  compute functions, otherwise the previous observations will be used.
	 */
	public void startIndividualObservations() {
		individualObservations = new Vector<double[]>();
	}
	
	public void addObservation(double observation[]) {
		individualObservations.add(observation);
	}
	
	public void endIndividualObservations() throws Exception {
		double[][] data = new double[individualObservations.size()][];
		for (int t = 0; t < data.length; t++) {
			data[t] = individualObservations.elementAt(t);
		}
		// Allow vector to be reclaimed
		individualObservations = null;
		setObservations(data);
	}


	/**
	 * Compute the norms for each marginal time series
	 *
	 */
	protected void computeNorms() {
		
		norms = new double[V][N][N];
		for (int t = 0; t < N; t++) {
			// Compute the norms from t to all other time points
			double[][] normsForT = EuclideanUtils.computeNorms(data, t);
			for (int t2 = 0; t2 < N; t2++) {
				for (int v = 0; v < V; v++) {
					norms[v][t][t2] = normsForT[t2][v];
				}
			}
		}
	}
	
	public abstract double computeAverageLocalOfObservations() throws Exception;

	/**
	 * Compute what the average MI would look like were all time series reordered
	 *  as per the array of time indices in reordering.
	 * The reordering array contains the reordering for each marginal variable (first index).
	 * The user should ensure that all values 0..N-1 are represented exactly once in the
	 *  array reordering and that no other values are included here.
	 * 
	 * @param reordering the specific new orderings to use. First index is the variable number
	 *  (minus 1, since we don't reorder the first variable),
	 *  second index is the time step, the value is the reordered time step to use
	 *  for that variable at the given time step.
	 * @return
	 * @throws Exception
	 */
	public abstract double computeAverageLocalOfObservations(int[][] reordering) throws Exception;

	/**
	 * Compute the significance of the multi-information of the previously supplied observations.
	 * We destroy the p(x,y,z,..) correlations, while retaining the p(x), p(y),.. marginals, to check how
	 *  significant this multi-information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param numPermutationsToCheck
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public synchronized EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		int[][][] newOrderings = new int[numPermutationsToCheck][][];
		// Generate numPermutationsToCheck * V permutations of 0 .. data.length-1
		for (int n = 0; n < numPermutationsToCheck; n++) {
			// (Not necessary to check for distinct random perturbations)
			newOrderings[n] = rg.generateRandomPerturbations(data.length, V-1);
		}
		return computeSignificance(newOrderings);
	}
	
	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y,z,..) correlations, while retaining the p(x), p(y),.. marginals, to check how
	 *  significant this mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param newOrderings the specific new orderings to use. First index is the reordering index,
	 *  second index is the variable number (minus 1, since we don't reorder the first variable),
	 *  third index is the reordered variable number for that position.
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][][] newOrderings) throws Exception {
		int numPermutationsToCheck = newOrderings.length;
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		// Store the real observations and their MI:
		double actualMI = mi;
		
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		
		int countWhereMiIsMoreSignificantThanOriginal = 0;
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Compute the MI under this reordering
			double newMI = computeAverageLocalOfObservations(newOrderings[i]);
			measDistribution.distribution[i] = newMI;
			if (debug){
				System.out.println("New MI was " + newMI);
			}
			if (newMI >= actualMI) {
				countWhereMiIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the actual MI and the observations
		mi = actualMI;

		// And return the significance
		measDistribution.pValue = (double) countWhereMiIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualMI;
		return measDistribution;
	}

	public abstract double[] computeLocalOfPreviousObservations() throws Exception;

	public double[] computeLocalUsingPreviousObservations(double[][] states) throws Exception {
		// TODO If this is implemented, will need to normalise the incoming
		//  observations the same way that previously supplied ones were
		//  normalised (if they were normalised, that is)
		throw new Exception("Local method not implemented yet");
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return mi;
	}
	
	public abstract String printConstants(int N) throws Exception;

	/**
	 * Utility to take a reordering matrix and return the array of reordered time indices from
	 *  which to find the reordered data to be inserted at timeStep.
	 * 
	 * @param reordering the specific new orderings to use. First index is the variable number
	 *  (can be for all variables, or one less than all if the first is not to be reordered),
	 *  second index is the time step, the value is the reordered time step to use
	 *  for that variable at the given time step.
	 *  If null, no reordering is performed.
	 * @param timeStep
	 * @return
	 */
	protected int[] reorderedTimeStepsForEachMarginal(int[][] reordering, int timeStep) {
		// Create storage for the reordered time steps for the variables
		int[] tForEachMarginal = new int[V];
		if (reordering == null) {
			// We're not reordering
			for (int v = 0; v < V; v++) {
				tForEachMarginal[v] = timeStep;
			}
		} else {
			boolean reorderingFirstColumn = (reordering.length == V);
			int reorderIndex = 0;
			// Handle the first column
			if (reorderingFirstColumn) {
				tForEachMarginal[0] = reordering[reorderIndex++][timeStep];
			} else {
				tForEachMarginal[0] = timeStep;
			}
			// Handle subsequent columns
			for (int v = 1; v < V; v++) {
				tForEachMarginal[v] = reordering[reorderIndex++][timeStep];
			}
		}
		return tForEachMarginal;
	}
}
