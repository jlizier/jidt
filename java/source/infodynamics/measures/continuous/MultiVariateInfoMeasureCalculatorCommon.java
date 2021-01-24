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

package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

import java.util.Arrays;
import java.util.Random;
import java.util.Vector;

/**
 * Implements a base class with common functionality for child class
 * implementations of multivariate information measures via various estimators. 
 *
 * <p>Multivariate information measures are functionals of probability
 * distributions over <code>R^n</code>, and typical examples include multi-information
 * (a.k.a. total correlation), dual total correlation, O-information, and connected
 * information.</p>
 *
 * <p>These measures can be computed via different kinds of estimators, such as
 * linear-gaussian, KSG estimators, etc (see the child classes linked above).
 * </p>
 * 
 * <p>
 * Usage of the child classes is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator using {@link #initialise(int)};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[][])}
 * 					for calculations based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to {@link #addObservations(double[][])} or
 * 							{@link #addObservation(double[])}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average measure:
 *					{@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local values for these samples:
 * 					{@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local values for a specific set of samples:
 * 					{@link #computeLocalUsingPreviousObservations(double[][])}</li>
 * 			</ul>
 * 		</li>
 * 		<li>
 * 		Return to step 2 or 3 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Rosas, F., Mediano, P., Gastpar, M, Jensen, H., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.100.032305">"Quantifying high-order
 *   interdependencies via multivariate extensions of the mutual information"</a>,
 *   Physical Review E 100, (2019) 032305.</li>
 * </ul>
 * 
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public abstract class MultiVariateInfoMeasureCalculatorCommon
  implements InfoMeasureCalculatorContinuous {

	/**
	 * Number of joint variables to consider
	 */
	protected int dimensions = 1;
	/**
	 * Number of samples supplied
	 */
	protected int totalObservations = 0;
	/**
	 * Whether we are in debug mode
	 */
	protected boolean debug = false;
	/**
	 * Cached supplied observations
	 */
	protected double[][] observations;
	/**
	 * Set of individually supplied observations
	 */
	protected Vector<double[]> individualObservations;
	/**
	 * Whether the user has supplied more than one (disjoint) set of samples
	 */
	protected boolean addedMoreThanOneObservationSet;
	/**
	 * Whether the measure has been computed for the latest supplied data
	 */
	protected boolean isComputed = false;
	/**
	 * Cached last the measure value calculated
	 */
	protected double lastAverage;

	/**
	 * Whether to normalise incoming values
	 */
	protected boolean normalise = true;

	/**
	 * Property name for whether to normalise incoming values to mean 0,
	 * standard deviation 1 (default true)
	 */
	public static final String PROP_NORMALISE = "NORMALISE";


	/**
	 * Initialise the calculator for (re-)use, with the existing
	 * (or default) values of parameters, with number of 
	 * joint variables specified.
	 * Clears any PDFs of previously supplied observations.
	 *
	 * @param dimensions the number of joint variables to consider
	 */
	public void initialise() {
		initialise(dimensions);
	}
	
	/**
	 * Initialise the calculator for (re-)use, with the existing
	 * (or default) values of parameters, with number of 
	 * joint variables specified.
	 * Clears an PDFs of previously supplied observations.
	 *
	 * @param dimensions the number of joint variables to consider
	 */
	public void initialise(int dimensions) {
		this.dimensions = dimensions;
		lastAverage = 0.0;
		totalObservations = 0;
		isComputed = false;
		observations = null;
		addedMoreThanOneObservationSet = false;
	}

	/**
	 * Set properties for the calculator.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #PROP_NORMALISE} -- whether to normalise the incoming variables 
	 * 			to mean 0, standard deviation 1, or not (default false).</li>
	 * </ul>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		} else {
			// No property was set here
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			return Boolean.toString(normalise);
		} else {
			return null;
		}
	}

	/**
	 * Sets a single series from which to compute the PDF for the instantiated measure..
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * <p>The supplied series may be a time-series, or may be simply
	 * a set of separate observations
	 * without a time interpretation.</p>
	 * 
	 * <p>Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated).</p> 
	 * 
	 * @param observations series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 * @throws Exception
	 */
	public void setObservations(double[][] observations) throws Exception {
		startAddObservations();
		addObservations(observations);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	/**
	 * Signal that we will add in the samples for computing the PDF 
	 * from several disjoint time-series or trials via calls to
	 * {@link #addObservation(double[])} or {@link #addObservations(double[][])}
	 * rather than {@link #setDebug(boolean)}.
	 */
	public void startAddObservations() {
		individualObservations = new Vector<double[]>();
	}

	/**
	 * <p>Adds a new (single) observation to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p>Note that the arrays must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param observation a single multivariate observation
	 *  (index is variable number)
	 */
	public void addObservation(double[] observation) {
		if (individualObservations.size() > 0) {
			addedMoreThanOneObservationSet = true;
		}
		individualObservations.add(observation);
	}

	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p>Note that the arrays must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param observations series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 */
	public void addObservations(double[][] observations) {
		// This implementation is not particularly efficient,
		//  however it will suffice for now.
		for (int s = 0; s < observations.length; s++) {
			addObservation(observations[s]);
		}
	}

	/**
	 * {@inheritDoc} 
	 * 
	 * This class provides a basic implementation, generating
	 * the internal set of samples in observations; child classes
	 * should then process these observations as required.
	 * 
	 */
	public void finaliseAddObservations() throws Exception {
		observations = new double[individualObservations.size()][];
		for (int t = 0; t < observations.length; t++) {
			observations[t] = individualObservations.elementAt(t);
		}
		// Allow vector to be reclaimed
		individualObservations = null;
		
		if (observations[0].length != dimensions) {
			throw new Exception("Incorrect number of dimensions " + observations[0].length +
					" in supplied observations (expected " + dimensions + ")");
		}
		totalObservations = observations.length;
	}
	
	/**
	 * Generate a resampled distribution of what the measure would look like,
	 * under a null hypothesis that the individual values of each
	 * variable in the 
	 * samples have no relation to each other.
	 * That is, we destroy the p(x,y,z,..) correlations, while
	 * retaining the p(x), p(y),.. marginals, to check how
	 *  significant this measure actually was.
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
	 * creates <i>random</i> shufflings of the next values for the surrogate MultiInfo
	 * calculations.</p>
	 * 
	 * @param numPermutationsToCheck number of surrogate samples to permute
	 *  to generate the distribution.
	 * @return the distribution of surrogate measure values under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		int[][][] newOrderings = new int[numPermutationsToCheck][][];
		// Generate numPermutationsToCheck * (dimensions-1) permutations of 0 .. data.length-1
		for (int n = 0; n < numPermutationsToCheck; n++) {
			// (Not necessary to check for distinct random perturbations)
			newOrderings[n] = rg.generateRandomPerturbations(totalObservations, dimensions-1);
		}
		return computeSignificance(newOrderings);
	}

	/**
	 * Generate a resampled distribution of what the measure would look like,
	 * under a null hypothesis that the individual values of each
	 * variable in the 
	 * samples have no relation to eachother.
	 * That is, we destroy the p(x,y,z,..) correlations, while
	 * retaining the p(x), p(y),.. marginals, to check how
	 *  significant this measure actually was.
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
	 * @return the distribution of surrogate measure values under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception where the length of each permutation in newOrderings
	 *   is not equal to the number N samples that were previously supplied.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		if (!isComputed) {
			computeAverageLocalOfObservations();
		}
		
		// Store the real observations and their measure value:
		double actualMeasure = lastAverage;
		
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		
		int countWhereSurrogateIsMoreSignificantThanOriginal = 0;
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Compute the measure under this reordering
			double newMeasure = computeAverageLocalOfObservations(newOrderings[i]);
			measDistribution.distribution[i] = newMeasure;
			if (debug){
				System.out.println("New measure value was " + newMeasure);
			}
			if (newMeasure >= actualMeasure) {
				countWhereSurrogateIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the actual measure and the observations
		lastAverage = actualMeasure;

		// And return the significance
		measDistribution.pValue = (double) countWhereSurrogateIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualMeasure;
		return measDistribution;
	}

  /**
   * Compute what the measure would look like were all time series (bar the first)
   * reordered as per the array of time indices in newOrdering.
	 * 
   * <p>The reordering array contains the reordering for each marginal variable
   * (first index). The user should ensure that all values 0..N-1 are
   * represented exactly once in the array reordering and that no other values
   * are included here.</p>
	 *  
   * <p>Note that if several disjoint time-series have been added as
   * observations using {@link #addObservations(double[])} etc., then these
   * separate "trials" will be mixed up in the generation of a shuffled source
   * series here.</p>
	 * 
   * <p>This method is primarily intended for use in {@link
   * #computeSignificance(int[][])} however has been made public in case users
   * wish to access it.</p>
	 * 
   * @param newOrdering the specific permuted new orderings to use. First index
   * is the variable number (minus 1, since we don't reorder the first
   * variable), second index is the time step, the value is the reordered time
   * step to use for that variable at the given time step.  The values must be
   * an array of length N (where
   *  would be the value returned by {@link #getNumObservations()}), containing
   *  a permutation of the values in 0..(N-1).  If null, no reordering is
   *  performed.
   * @return what the average measure would look like under this reordering
   * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[][] newOrdering)
			throws Exception {
		
		if (newOrdering == null) {
			return computeAverageLocalOfObservations();
		}
		
		// Take a clone of the object to compute the measure of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays)
		MultiVariateInfoMeasureCalculatorCommon surrogateCalculator =
				(MultiVariateInfoMeasureCalculatorCommon) this.clone();
		
		// Generate a new re-ordered source data
		double[][] shuffledData =
				MatrixUtils.reorderDataForVariables(
					observations, newOrdering);
		// Perform new initialisations
		surrogateCalculator.initialise(dimensions);
		// Set new observations
		surrogateCalculator.setObservations(shuffledData);
		// Compute the MI
		return surrogateCalculator.computeAverageLocalOfObservations();
	}

	/**
	 * Calculates the local measure at every sample provided since the last time the
   * calculator was initialised.
	 * 
	 * @return the "time-series" of local measure values in nats (not bits!)
	 * @throws Exception
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		// Cannot do if observations haven't been set
		if (observations == null) {
			throw new Exception("Cannot compute local values of previous observations " +
					"if they have not been set!");
		}
		
		return computeLocalUsingPreviousObservations(observations);
	}

	/**
	 * Compute the local measure values for each of the
	 * supplied samples in <code>states</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>states</code>
	 * (unless they were
	 * some of the previously supplied samples).</p>
	 * 
	 * @param states series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 * @return the series of local measure values.
	 * @throws Exception
	 */
	public abstract double[] computeLocalUsingPreviousObservations(double states[][])
		throws Exception;

  /**
   * Shortcut method to initialise the calculator, set observations and compute
   * the average measure in one line.
   *
	 * @param new_observations series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
   */
  public double compute(double[][] new_observations) throws Exception {
    initialise(new_observations[0].length);
    setObservations(new_observations);
    return computeAverageLocalOfObservations();
  }

  /**
   * Shortcut method to initialise the calculator, set observations and compute
   * the local measure in one line.
   *
	 * @param new_observations series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
   */
  public double[] computeLocals(double[][] new_observations) throws Exception {
    initialise(new_observations[0].length);
    setObservations(new_observations);
    return computeLocalOfPreviousObservations();
  }
	
	public int getNumObservations() throws Exception {
		return totalObservations;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return lastAverage;
	}

  /**
   * Returns an <code>int[]</code> array with all integers from 0 to
   * <code>N-1</code>, except <code>idx</code>.
   *
   * <p>This method is primarily intended for internal use (to extract blocks
   * of covariance matrices excluding one variable).</p>
	 * 
   * @param idx index of integer to omit
   * @param N upper limit of the integer array
   */
  protected int[] allExcept(int idx, int N) {
    boolean[] v = new boolean[N];
    Arrays.fill(v, true);
    v[idx] = false;

    int[] v2 = new int[N - 1];
    int counter = 0;
    for (int i = 0; i < N; i++) {
      if (v[i]) {
        v2[counter] = i;
        counter++;
      }
    }

    return v2;
  }

}

