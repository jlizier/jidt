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
 * <p>Computes the differential multi-information of a given multivariate set of
 *  observations (implementing {@link MultiInfoCalculator}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below).
 *  This is an abstract class to gather common functionality between the two
 *  algorithms defined by Kraskov et al.
 *  Two child classes {@link MultiInfoCalculatorKraskov1} and
 *  {@link MultiInfoCalculatorKraskov2} then
 *  actually implement the two algorithms in the Kraskov et al. paper</p>
 *    
 * <p>Usage is as per the paradigm outlined for {@link MultiInfoCalculator},
 * with:
 * <ul>
 * 	<li>For constructors see the child classes.</li>
 *  <li>Further properties are defined in {@link #setProperty(String, String)}.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 * 
 * <p>
 * TODO Add fast nearest neighbour searches to the child classes
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class MultiInfoCalculatorKraskov implements
		MultiInfoCalculator {

	/**
	 * we compute distances to the kth nearest neighbour
	 */
	protected int k = 4;
	/**
	 * Cached observations
	 */
	protected double[][] data;
	/**
	 * Whether we are in debug mode
	 */
	protected boolean debug;
	/**
	 * Last average multi-info computed
	 */
	protected double mi;
	protected boolean miComputed;
	
	/**
	 * Set of individually supplied observations
	 */
	private Vector<double[]> individualObservations;

	/**
	 * number of observations supplied
	 */
	protected int N;
	/**
	 * number of joint variables
	 */
	protected int V; 
	
	/**
	 * Calculator for the norm between data points
	 */
	protected EuclideanUtils normCalculator;
	/**
	 * Cached norms for each marginal variable from each observation to each other one
	 */
	protected double[][][] norms;
	/**
	 * Whether to keep the norms each time (making reordering very quick)
	 * (Should only be set to false for testing)
	 */
	protected boolean tryKeepAllPairsNorms = true;
	public static int MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM = 4000;
	
	/**
	 * Property name for the number of K nearest neighbours used in
	 * the KSG algorithm (default 4).
	 */
	public final static String PROP_K = "k";
	/**
	 * Property name for what type of norm to use between data points
	 *  for each marginal variable -- Options are defined by 
	 *  {@link EuclideanUtils#setNormToUse(String)} and the
	 *  default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 */
	public final static String PROP_NORM_TYPE = "NORM_TYPE";
	/**
	 * Property name for whether to keep the norms
	 * each time, making reordering very quick
	 * (default is true, Should only be set to false for testing)
	 */
	public final static String PROP_TRY_TO_KEEP_ALL_PAIRS_NORM = "TRY_KEEP_ALL_PAIRS_NORM";
	
	// TODO Add a NORMALISE property
	
	/**
	 * Construct an instance
	 */
	public MultiInfoCalculatorKraskov() {
		super();
		normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
	}

	@Override
	public void initialise(int dimensions) {
		V = dimensions;
		mi = 0.0;
		miComputed = false;
		norms = null;
		data = null;
	}

	/**
	 * Sets properties for the KSG multi-info calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 *  <li>{@link #PROP_K} -- number of k nearest neighbours to use in joint kernel space
	 *      in the KSG algorithm (default is 4).</li>
	 * 	<li>{@link #PROP_NORM_TYPE}</li> -- normalization type to apply to 
	 * 		working out the norms between the points in each marginal space.
	 * 		Options are defined by {@link EuclideanUtils#setNormToUse(String)} -
	 * 		default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 *  <li>{@link #PROP_TRY_TO_KEEP_ALL_PAIRS_NORM} -- whether to keep the norms
	 * 		each time, making reordering very quick
	 * 		(default is true, Should only be set to false for testing)</li>
	 * </ul>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normCalculator.setNormToUse(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_TRY_TO_KEEP_ALL_PAIRS_NORM)) {
			tryKeepAllPairsNorms = Boolean.parseBoolean(propertyValue);
		}
	}

	@Override
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
	 * Set observations from two separate time series: join the rows at each time step
	 *  together to make a joint vector, then effectively call {@link #setObservations(double[][])}
	 * 
	 * @param observations1 observations for first few variables
	 * @param observations2 observations for the other variables
	 * @see #setObservations(double[][])
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
	
	@Override
	public void startAddObservations() {
		individualObservations = new Vector<double[]>();
	}
	
	@Override
	public void addObservation(double observation[]) {
		individualObservations.add(observation);
	}
	
	@Override
	public void addObservations(double[][] observations) {
		// This implementation is not particularly efficient,
		//  however for the little use this calculator will 
		//  attract, it will suffice.
		for (int s = 0; s < observations.length; s++) {
			addObservation(observations[s]);
		}
	}

	@Override
	public void finaliseAddObservations() throws Exception {
		double[][] data = new double[individualObservations.size()][];
		for (int t = 0; t < data.length; t++) {
			data[t] = individualObservations.elementAt(t);
		}
		// Allow vector to be reclaimed
		individualObservations = null;
		setObservations(data);
	}


	/**
	 * Compute the norms between each observation
	 * for each marginal time series
	 * and cache them
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
	
	/**
	 * Compute what the average multi-info would look like were all time series
	 *  (bar the first) reordered
	 *  as per the array of time indices in reordering.
	 * The reordering array contains the reordering for each marginal variable (first index).
	 * The user should ensure that all values 0..N-1 are represented exactly once in the
	 *  array reordering and that no other values are included here.
	 * 
	 * @param reordering the specific new orderings to use. First index is the variable number
	 *  (minus 1, since we don't reorder the first variable),
	 *  second index is the time step, the value is the reordered time step to use
	 *  for that variable at the given time step.
	 *  If null, no reordering is performed.
	 * @return what the average multi-info would look like under this reordering
	 * @throws Exception
	 */
	public abstract double computeAverageLocalOfObservations(int[][] reordering) throws Exception;

	/**
	 * Generate a bootstrapped distribution of what the multi-information would look like,
	 * under a null hypothesis that the individual values of each
	 * variable in the 
	 * samples have no relation to eachother.
	 * That is, we destroy the p(x,y,z,..) correlations, while
	 * retaining the p(x), p(y),.. marginals, to check how
	 *  significant this multi-information actually was.
	 *  
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * we are extending that here.
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p>This method (in contrast to {@link #computeSignificance(int[][][])})
	 * creates <i>random</i> shufflings of the next values for the surrogate AIS
	 * calculations.</p>
	 * 
	 * @param numPermutationsToCheck number of surrogate samples to bootstrap
	 *  to generate the distribution.
	 * @return the distribution of surrogate multi-info values under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception
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
	 * Generate a bootstrapped distribution of what the multi-information would look like,
	 * under a null hypothesis that the individual values of each
	 * variable in the 
	 * samples have no relation to eachother.
	 * That is, we destroy the p(x,y,z,..) correlations, while
	 * retaining the p(x), p(y),.. marginals, to check how
	 *  significant this multi-information actually was.
	 *  
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * we are extending that here.
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p>This method (in contrast to {@link #computeSignificance(int)})
	 * allows the user to specify how to construct the surrogates,
	 * such that repeatable results may be obtained.</p>
	 * 
	 * @param newOrderings a specification of how to shuffle the values
	 *  to create the surrogates to generate the distribution with. The first
	 *  index is the permutation number (i.e. newOrderings.length is the number
	 *  of surrogate samples we use to bootstrap to generate the distribution here.)
	 *  The second index is the variable number (minus 1, since we don't reorder
	 *  the first variable),
	 *  Each array newOrderings[i][v] should be an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 * @return the distribution of surrogate multi-info values under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception where the length of each permutation in newOrderings
	 *   is not equal to the number N samples that were previously supplied.
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
	
	/**
	 * Utility function used for debugging, printing digamma constants
	 * 
	 * @param N
	 * @return
	 * @throws Exception
	 */
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
	 * @return array of reordered time indices from
	 *  which to find the reordered data to be inserted at timeStep
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
