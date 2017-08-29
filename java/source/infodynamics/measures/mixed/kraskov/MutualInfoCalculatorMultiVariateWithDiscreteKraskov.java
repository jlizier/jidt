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

package infodynamics.measures.mixed.kraskov;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.mixed.MutualInfoCalculatorMultiVariateWithDiscrete;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;
import infodynamics.utils.KdTree;

import java.util.Iterator;
import java.util.Vector;
import java.util.Arrays;

/**
 * <p>Compute the Mutual Information between a vector of continuous variables and discrete
 *  variable using the Kraskov estimation method.</p>
 * <p>Uses Kraskov method type 2, since type 1 only looks at points with
 * distances strictly less than the kth variable, which won't work for one marginal
 * being discrete.</p>
 * <p>I have noticed that there are quite large bias negative values here
 * where small K is used (e.g. for binary data splitting continuous into
 * two distinct groups, 1400 observations, K=4 has bias ~ -.35)</p>
 * 
 * <p>These calculators are <b>EXPERIMENTAL</b> -- not properly tested,
 * and not well documented. The intended calling pattern is similar to
 * {@link MutualInfoCalculatorMultiVariate}
 * </p>
 * 
 * <p>This calculator effectively implements the nearest neighbour method
 * laid out <a href="http://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0087357">here</a>
 * (see Ross reference below),
 * but we have subtract out the sample point from the count in the discrete space,
 * following the original KSG method.
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li> B.C. Ross, "Mutual Information between Discrete and Continuous Data Sets",
 * 	PLoS ONE 9(2): e87357. doi:
 * <a href="http://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0087357">10.1371/journal.pone.0087357</a></li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MutualInfoCalculatorMultiVariateWithDiscreteKraskov implements MutualInfoCalculatorMultiVariateWithDiscrete, Cloneable {

	/**
	 * we compute distances to the kth neighbour
	 */
	protected int k = 4;
	/**
	 * The set of continuous data observations, retained in case the user wants
   * to retrieve the local entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link addObservations(double[][], int[])} functions.
	 */
	protected double[][] continuousData;
	/**
	 * The set of discrete data observations, retained in case the user wants
   * to retrieve the local entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link addObservations(double[][], int[])} functions.
	 */
	protected int[] discreteData;
  /**
   * Counts of how many observations belong to each discrete bin.
   */
	protected int[] counts;
	/**
	 * Number of possible states of the discrete variable
	 */
	protected int base;
	/**
	 * Number of dimenions of the joint continuous variable
	 */
	protected int dimensions;
  /**
   * Number of observations added to the calculator
   */
  protected int totalObservations;
  /**
   * Whether to print debug information
   */
	protected boolean debug;
  /**
   * Last computed value of MI
   */
	protected double mi;
  /**
   * Whether MI has been computed with the current observations
   */
	protected boolean miComputed;
  /**
   * Saved digammaN for average and local computations
   */
  protected double digammaN;
  /**
   * Saved digammaK for average and local computations
   */
  protected double digammaK;
	/**
	 * Storage for continuous observations supplied via {@link #addObservations(double[][], int[])}
	 * type calls
	 */
	protected Vector<double[][]> vectorOfContinuousObservations;
	/**
	 * Storage for discrete observations supplied via {@link #addObservations(double[][], int[])}
	 * type calls
	 */
	protected Vector<int[]> vectorOfDiscreteObservations;
  /**
   * KdTree for range searches in full continuous dataset
   */
  protected KdTree kdTreeJoint;
  /**
   * KdTrees for neighbour searches in binned continuous datasets
   */
  protected KdTree[] kdTreeBins;
  /**
   * The norm type in use (see {@link #PROP_NORM_TYPE})
   */
  protected int normType = EuclideanUtils.NORM_MAX_NORM;
  /**
   * Property name for the number of K nearest neighbours used in
   * the KSG algorithm in the full joint space (default 4).
   */
	public final static String PROP_K = "k";
  /**
   * Property name for what type of norm to use between data points
   *  for each marginal variable -- Options are defined by 
   *  {@link KdTree#setNormType(String)} and the
   *  default is {@link EuclideanUtils#NORM_MAX_NORM}.
   */
	public final static String PROP_NORM_TYPE = "NORM_TYPE";	
  /**
   * Property name for whether to normalise the incoming data to 
   * mean 0, standard deviation 1 (default true)
   */
	public static final String PROP_NORMALISE = "NORMALISE";
	/**
	 * Property name for time difference between source and destination (0 by default,
	 * must be >= 0)
	 */ 
	public static final String PROP_TIME_DIFF = "TIME_DIFF";
	/**
	 * Track whether we're going to normalise the joint variables individually
	 */
	protected boolean normalise = true;
	/**
	 * Track the means of the joint variables if we are normalising them
	 */
	protected double[] means;
	/**
	 * Track the std devs of the joint variables if we are normalising them
	 */
	protected double[] stds;
	/**
	 * Time difference from the source to the destination observations
	 *  (ie destination lags the source by this time:
	 *  	we compute I(source_{n}; dest_{n+timeDiff}).
	 * (Note that our internal sourceObservations and destObservations
	 *  are adjusted so that there is no timeDiff between them).
	 */
	protected int timeDiff = 0;

	public MutualInfoCalculatorMultiVariateWithDiscreteKraskov() {
		super();
	}

	/**
	 * Initialise the calculator.
	 * 
	 * @param dimensions number of joint continuous variables
	 * @param base number of discrete states
	 */
	public void initialise(int dimensions, int base) {
		mi = 0.0;
		miComputed = false;
    totalObservations = 0;
		continuousData = null;
		means = null;
		stds = null;
		discreteData = null;
		this.dimensions = dimensions;
		this.base = base;
    kdTreeJoint = null;
    kdTreeBins = null;
	}

	/**
	 * Sets properties for the calculator.
	 * Valid properties include:
	 * <ul>
	 *  <li>{@link #PROP_K} - number of neighbouring points in joint kernel space (default 4)</li>
	 * 	<li>{@link #PROP_NORM_TYPE}</li> - normalization type to apply to 
	 * 		working out the norms between the points in each marginal space.
	 * 		Options are defined by {@link EuclideanUtils#setNormToUse(String)} -
	 * 		default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 *  <li>{@link #PROP_NORMALISE} - whether to normalise the individual
	 *      variables (true by default)</li>
	 * 	<li>{@link #PROP_TIME_DIFF} - Time difference between source and
   *     	destination (0 by default). Must be >= 0.</li>
	 * </ul>
	 * 
	 * @param propertyName name of the property to set
	 * @param propertyValue value to set on that property
	 */
	public void setProperty(String propertyName, String propertyValue)
      throws Exception {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
      normType = KdTree.validateNormType(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_TIME_DIFF)) {
			int val = Integer.parseInt(propertyValue);
			if (val < 0) {
				throw new Exception("Time difference must be >= 0. Flip data1 and data2 around if required.");
			}
      timeDiff = val;
		}
	}

	public void addObservations(double[][] continuousObservations,
      int[] discreteObservations) throws Exception {
		if (vectorOfContinuousObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (continuousObservations.length != discreteObservations.length) {
			throw new Exception("Time steps for observations2 " +
					discreteObservations.length + " does not match the length " +
					"of observations1 " + continuousObservations.length);
		}
		if (continuousObservations[0].length != dimensions) {
			throw new Exception("The continuous observations do not have the expected number of variables (" + dimensions + ")");
		}
    if (continuousObservations.length > timeDiff) {
      vectorOfContinuousObservations.add(continuousObservations);
      vectorOfDiscreteObservations.add(discreteObservations);
    }
	}

	public void addObservations(double[][] source, double[][] destination, int startTime, int numTimeSteps) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination, boolean[] sourceValid, boolean[] destValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination, boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void startAddObservations() {
		vectorOfContinuousObservations = new Vector<double[][]>();
		vectorOfDiscreteObservations   = new Vector<int[]>();
	}

	public void finaliseAddObservations() throws Exception {
		if (vectorOfContinuousObservations.size() < 1) {
			throw new Exception("Cannot compute MI with a null set of data");
		}

		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[][] destination : vectorOfContinuousObservations) {
			totalObservations += destination.length - timeDiff;
		}
		continuousData = new double[totalObservations][dimensions];
		discreteData   = new int[totalObservations];
		
		// Construct the joint vectors from the given observations
		//  (removing redundant data which is outside any timeDiff)
		int startObservation = 0;
		Iterator<double[][]> iterator = vectorOfContinuousObservations.iterator();
		for (int[] dct : vectorOfDiscreteObservations) {
			double[][] cnt = iterator.next();
			// Copy the data from these given observations into our master 
			//  array, aligning them incorporating the timeDiff:
			MatrixUtils.arrayCopy(cnt, 0, 0,
					continuousData, startObservation, 0,
					cnt.length - timeDiff, dimensions);
      System.arraycopy(dct, timeDiff, discreteData, startObservation, dct.length - timeDiff);
			startObservation += cnt.length - timeDiff;
		}
		// We don't need to keep the vectors of observation sets anymore:
		vectorOfContinuousObservations = null;
		vectorOfDiscreteObservations   = null;

		if (normalise) {
			// We need to keep the means/stds ready to normalise local values
			//  that are supplied later:
			means = MatrixUtils.means(continuousData);
			stds = MatrixUtils.stdDevs(continuousData, means);
			MatrixUtils.normalise(continuousData, means, stds);
		}

		// count the discrete states:
		counts = new int[base];
    try {
      for (int t = 0; t < discreteData.length; t++) {
        counts[discreteData[t]]++;
      }
    } catch (ArrayIndexOutOfBoundsException e) {
      totalObservations = 0;
      continuousData = null;
      discreteData   = null;
      throw new RuntimeException("Values of the discrete variable must range from 0 to base-1");
    }

		for (int b = 0; b < counts.length; b++) {
			if (counts[b] < k) {
				throw new RuntimeException("This implementation assumes there are at least k items in each discrete bin");
			}
		}

    digammaN = MathsUtils.digamma(totalObservations);
    digammaK = MathsUtils.digamma(k);

    ensureKdTreesConstructed();
	}


	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations) throws Exception {
    startAddObservations();
    addObservations(continuousObservations, discreteObservations);
    finaliseAddObservations();
	}

  /**
   * Internal method to ensure that the Kd-tree data structures to represent the
   * observational data have been constructed (should be called prior to attempting
   * to use these data structures)
   */
  public void ensureKdTreesConstructed() {
    if (kdTreeJoint == null) {
      kdTreeJoint = new KdTree(continuousData);
      kdTreeJoint.setNormType(normType);
    }
    if (kdTreeBins == null) {
      kdTreeBins = new KdTree[base];
      for (int b = 0; b < base; b++) {
        kdTreeBins[b] = new KdTree(
            MatrixUtils.extractSelectedPointsMatchingCondition(
              continuousData, discreteData, b));
        kdTreeBins[b].setNormType(normType);
      }
    }
  }
	
	/**
	 * Compute what the average MI would look like were the second time series reordered
	 *  as per the array of time indices in reordering.
	 * The user should ensure that all values 0..N-1 are represented exactly once in the
	 *  array reordering and that no other values are included here. 
	 * 
	 * @param reordering
	 * @return
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[] newOrdering)
			throws Exception {
		
		if (newOrdering == null) {
			return computeAverageLocalOfObservations();
		}

		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays)
		MutualInfoCalculatorMultiVariateWithDiscreteKraskov miSurrogateCalculator =
				(MutualInfoCalculatorMultiVariateWithDiscreteKraskov) this.clone();
		
		// Generate a new re-ordered data set
		int[] shuffledDiscreteData = MatrixUtils.extractSelectedTimePoints(discreteData, newOrdering);
		// Perform new initialisations
		miSurrogateCalculator.initialise(dimensions, base);
		// Set new observations
		miSurrogateCalculator.setObservations(continuousData, shuffledDiscreteData);
		// Compute the MI
		return miSurrogateCalculator.computeAverageLocalOfObservations();
	}

  public double computeAverageLocalOfObservations() throws Exception {
    return computeFromObservations(false)[0];
  }

  public double[] computeLocalOfPreviousObservations() throws Exception {
    return computeFromObservations(true);
  }

	public double[] computeFromObservations(boolean returnLocals) throws Exception {

    // FIXME
    int dynCorrExclTime = 0;

    int[] cumcount = new int[base];
    Arrays.fill(cumcount, 0);


		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double avNx = 0;
		double avNy = 0;
		
		double testSum = 0.0; // Used for debugging prints
    int N = totalObservations;
    double[] locals = new double[N];

		for (int t = 0; t < N; t++) {

      // Compute eps_x for this time step:
      // find neighbours in the same discrete bin, then
      // count points within eps_x in the whole dataset.
      int b = discreteData[t];
      double eps_x = kdTreeBins[b].findKNearestNeighbours(k, cumcount[b]++).poll().distance;
      int n_x = kdTreeJoint.countPointsWithinOrOnR(t, eps_x, dynCorrExclTime);
      int n_y = counts[b] - 1;

			avNx += n_x;
			avNy += n_y;
			// And take the digamma before adding into the 
			//  average:
			double localSum = MathsUtils.digamma(n_x) + MathsUtils.digamma(n_y);
			averageDiGammas += localSum;

      // Don't need the 1/k correction here because the conditional entropy term
      //  is taken over the continuous space only. The correction is (m-1)/k
      //  for an entropy over m subspaces.
      // double localValue = MathsUtils.digamma(k) - 1.0/(double)k - localSum + MathsUtils.digamma(N);
      // Instead do:
      double localValue = digammaK + digammaN - localSum;
      locals[t] = localValue;

			if (debug) {
				testSum += localValue;
				if (dimensions == 1) {
					System.out.printf("t=%d: x=%.3f, eps_x=%.3f, n_x=%d, n_y=%d, local=%.3f, running total = %.5f\n",
							t, continuousData[t][0], eps_x, n_x, n_y, localValue, testSum);
				} else {
					System.out.printf("t=%d: eps_x=%.3f, n_x=%d, n_y=%d, local=%.3f, running total = %.5f\n",
							t, eps_x, n_x, n_y, localValue, testSum);
				}
			}
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f (-> digam=%.3f %.3f), Average n_y=%.3f (-> digam=%.3f)",
					avNx, MathsUtils.digamma((int) avNx), MathsUtils.digamma((int) avNx - 1), avNy, MathsUtils.digamma((int) avNy)));
			System.out.printf("Independent average num in joint box is %.3f\n", (avNx * avNy / (double) N));
			System.out.println(String.format("digamma(k)=%.3f - averageDiGammas=%.3f + digamma(N)=%.3f\n",
					digammaK, averageDiGammas, digammaN));
		}
		
		// Don't need the 1/k correction here because the conditional entropy term
		//  is taken over the continuous space only. The correction is (m-1)/k
		//  for an entropy over m subspaces.
		// mi = MathsUtils.digamma(k) - 1.0/(double)k - averageDiGammas + MathsUtils.digamma(N);
		// Instead, do:
		mi = digammaK + digammaN - averageDiGammas;
		miComputed = true;

    double[] returnValues;
    if (returnLocals) {
      returnValues = locals;
    } else {
      returnValues = new double[] {mi};
    }
		return returnValues;
	}


	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param numPermutationsToCheck
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 *  (i.e. 1 - CDF of our score)
	 */
	public synchronized EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(continuousData.length, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param newOrderings the specific new orderings to use
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception {
		
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
		measDistribution.actualValue = mi;
		return measDistribution;
	}

	/**
	 * <p>Compute the local MI values for the given observations,
	 *  using the previously set observations to compute the PDFs.
	 *  That is to say, we will evaluate the required counts for each
	 *  observation here based on the state space constructed from the
	 *  previous observations only (not these ones). This is non-standard,
	 *  and is not considered in the Kraskov paper. It is slightly unclear
	 *  how to account for the fact that this observation itself is not
	 *  in the data set used for computing the PDFs - I think though that
	 *  it should be ignored, since the data point was not counted in its
	 *  own counts in the standard version anyway.</p>
	 *  
	 *  <p><b>Importantly</b>, the supplied observations are intended to be new observations,
	 *  not those fed in to compute the PDFs from. There would be
	 *  some subtle changes necessary to accomodate computing locals on
	 *  the same data set (e.g. not counting the current point as one
	 *  of those within eps_x etc.). 
	 *  </p>
	 *  
	 *  @param continuousStates multivariate observations of the continuous variable
	 *    (1st index is time, 2nd is variable number)
	 *  @param discreteStates unvariate observations of the discrete variable
	 * 
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] continuousNewStates,
			int[] discreteNewStates) throws Exception {
		
		if (normalise) {
			// The stored observations continuousData have been normalised
			//  according to their stored means and stds; we need to 
			//  normalise the incoming observations the same way before
			//  comparing them
			continuousNewStates = MatrixUtils.normaliseIntoNewArray(
					continuousNewStates, means, stds);
		}

		double fixedPartOfLocals = digammaK + MathsUtils.digamma(totalObservations); // Use N_samplesForPdfs here because that's what would be in denominator of probability functions
		double testSum = 0.0, avNx = 0.0, avNy = 0.0;
    double[] locals = new double[discreteNewStates.length];
    for (int t = 0; t < discreteNewStates.length; t++) {

      int b = discreteNewStates[t];

      double[][] x = new double[][] {continuousNewStates[t]};
      double eps_x = kdTreeBins[b].findKNearestNeighbours(k, x).poll().distance;
      int n_x = kdTreeJoint.countPointsWithinR(x, eps_x, true);
			int n_y = counts[b];

			// Now compute the local value:
			locals[t] = fixedPartOfLocals -
					MathsUtils.digamma(n_x) - MathsUtils.digamma(n_y);

      if (debug) {
				testSum += locals[t];
        avNx += n_x;
        avNy += n_y;
				if (dimensions == 1) {
					System.out.printf("t=%d: x=%.3f, eps_x=%.3f, n_x=%d, n_y=%d, local=%.3f, running total = %.5f\n",
							t, continuousNewStates[t][0], eps_x, n_x, n_y, locals[t], testSum);
				} else {
					System.out.printf("t=%d: eps_x=%.3f, n_x=%d, n_y=%d, local=%.3f, running total = %.5f\n",
							t, eps_x, n_x, n_y, locals[t], testSum);
				}
      }
    }
		if (debug) {
			avNx /= (double) discreteNewStates.length;
			avNy /= (double) discreteNewStates.length;
			System.out.printf("Average n_x=%.3f, Average n_y=%.3f\n", avNx, avNy);
		}
		
		return locals;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return mi;
	}
	
	public int getNumObservations() {
		return totalObservations;
	}
}
