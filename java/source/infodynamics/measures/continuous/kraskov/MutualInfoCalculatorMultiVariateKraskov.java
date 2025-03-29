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

import java.util.Calendar;
import java.util.PriorityQueue;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.MutualInfoMultiVariateCommon;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.KdTree;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.NearestNeighbourSearcher;
import infodynamics.utils.NeighbourNodeData;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.NativeUtils;

/**
 * <p>Computes the differential mutual information of two given multivariate sets of
 *  observations (implementing {@link MutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below).
 *  The implementation is made using fast-neighbour searches with an
 *  underlying k-d tree algorithm.
 *  This is an abstract class to gather common functionality between the two
 *  algorithms defined by Kraskov et al.
 *  Two child classes {@link MutualInfoCalculatorMultiVariateKraskov1} and
 *  {@link MutualInfoCalculatorMultiVariateKraskov2} then
 *  actually implement the two algorithms in the Kraskov et al. paper</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link MutualInfoCalculatorMultiVariate},
 * with:
 * <ul>
 *  <li>For constructors see the child classes.</li>
 *  <li>Further properties are defined in {@link #setProperty(String, String)}.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 *  
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @author Ipek Özdemir
 */
public abstract class MutualInfoCalculatorMultiVariateKraskov
  extends MutualInfoMultiVariateCommon
  implements MutualInfoCalculatorMultiVariate {

  /**
   * we compute distances to the kth nearest neighbour
   */
  protected int k = 4;
  
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
   * Property name for the number of parallel threads to use in the
   *  computation (default is to use all available)
   */
  public static final String PROP_NUM_THREADS = "NUM_THREADS";
  /**
   * Valid property value for {@link #PROP_NUM_THREADS} to indicate
   *  that all available processors should be used. 
   */
  public static final String USE_ALL_THREADS = "USE_ALL";
  /**
   * Property name for the flag to enable or disable the GPU module.
   */
  public static final String PROP_USE_GPU = "USE_GPU";
  /**
   * Property name for the path to JIDT GPU library.
   *
   * Path must be full and contain the library filename.
   * Example: /home/johndoe/myfolder/libKraskov.so
   */
  public static final String PROP_GPU_LIBRARY_PATH = "GPU_LIBRARY_PATH";
  /**
   * Number of parallel threads to use in the computation;
   *  defaults to use all available.
   */
  protected int numThreads = Runtime.getRuntime().availableProcessors();
  /**
   * Private variable to record which KSG algorithm number
   *  this instance is implementing
   */
  protected boolean isAlgorithm1 = false;
  /**
   * Whether to enable the GPU module
   */
  protected boolean useGPU = false;
  /**
   * Path to JIDT GPU library
   */
  protected String gpuLibraryPath = "";
  /**
   * Check whether C native code has been loaded
   */
  protected boolean cudaLibraryLoaded = false;
  
  /**
   * protected k-d tree data structure (for fast nearest neighbour searches)
   *  representing the joint source-dest space
   */
  protected KdTree kdTreeJoint;
  /**
   * protected data structure (for fast nearest neighbour searches)
   *  representing the source space
   */
  protected NearestNeighbourSearcher nnSearcherSource;
  /**
   * protected data structure (for fast nearest neighbour searches)
   *  representing the dest space
   */
  protected NearestNeighbourSearcher nnSearcherDest;
  
  /**
   * Constant for digamma(k), with k the number of nearest neighbours selected
   */
  protected double digammaK;
  /**
   * Constant for digamma(N), with N the number of samples.
   */
  protected double digammaN;
  
  /**
   * Construct an instance of the KSG MI calculator
   */
  public MutualInfoCalculatorMultiVariateKraskov() {
    super();
    // Switch on adding noise to the data by default for the KSG estimator
    addNoise = true;
    noiseLevel = (double) 1e-8;
  }

  @Override
  public void initialise(int sourceDimensions, int destDimensions) {
    kdTreeJoint = null;
    nnSearcherSource = null;
    nnSearcherDest = null;
    super.initialise(sourceDimensions, destDimensions);
  }

  /**
   * Sets properties for the KSG MI calculator.
   *  New property values are not guaranteed to take effect until the next call
   *  to an initialise method. 
   *  
   * <p>Valid property names, and what their
   * values should represent, include:</p>
   * <ul>
   *  <li>{@link #PROP_K} -- number of k nearest neighbours to use in joint kernel space
   *      in the KSG algorithm (default is 4).</li>
   *  <li>{@link #PROP_NORM_TYPE} -- normalization type to apply to 
   *    working out the norms between the points in each marginal space.
   *    Options are defined by {@link KdTree#setNormType(String)} -
   *    default is {@link EuclideanUtils#NORM_MAX_NORM}.</li>
   *  <li>{@link #PROP_DYN_CORR_EXCL_TIME} -- a dynamics exclusion time window,
   *      also known as Theiler window (see Kantz and Schreiber);
   *      default is 0 which means no dynamic exclusion window.</li>
   *  <li>{@link #PROP_NUM_THREADS} -- the integer number of parallel threads
   *    to use in the computation. Can be passed as a string "USE_ALL"
   *      to use all available processors on the machine.
   *      Default is "USE_ALL".
   *  <li>any valid properties for {@link MutualInfoMultiVariateCommon#setProperty(String, String)}.</li>
   * </ul>
   * 
   * <p>Unknown property values are ignored.</p>
   * 
   * @param propertyName name of the property
   * @param propertyValue value of the property
   * @throws Exception for invalid property values
   */
  public void setProperty(String propertyName, String propertyValue) throws Exception {
    boolean propertySet = true;
    if (propertyName.equalsIgnoreCase(PROP_K)) {
      k = Integer.parseInt(propertyValue);
    } else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
      normType = KdTree.validateNormType(propertyValue);
    } else if (propertyName.equalsIgnoreCase(PROP_NUM_THREADS)) {
      if (propertyValue.equalsIgnoreCase(USE_ALL_THREADS)) {
        numThreads = Runtime.getRuntime().availableProcessors();
      } else { // otherwise the user has passed in an integer:
        numThreads = Integer.parseInt(propertyValue);
      }
    } else if (propertyName.equalsIgnoreCase(PROP_USE_GPU)) {
      useGPU = Boolean.parseBoolean(propertyValue);
    } else if (propertyName.equalsIgnoreCase(PROP_GPU_LIBRARY_PATH)) {
      gpuLibraryPath = propertyValue;
    } else {
      // No property was set here
      propertySet = false;
      // try the superclass:
      super.setProperty(propertyName, propertyValue);
    }
    if (debug && propertySet) {
      System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
          " to " + propertyValue);
    }
  }

  /**
   * Get property values for the calculator.
   * 
   * <p>Valid property names, and what their
   * values should represent, are the same as those for
   * {@link #setProperty(String, String)}</p>
   * 
   * <p>Unknown property values are responded to with a null return value.</p>
   * 
   * @param propertyName name of the property
   * @return current value of the property
   * @throws Exception for invalid property values
   */
  public String getProperty(String propertyName)
      throws Exception {
    
    if (propertyName.equalsIgnoreCase(PROP_K)) {
      return Integer.toString(k);
    } else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
      return KdTree.convertNormTypeToString(normType);
    } else if (propertyName.equalsIgnoreCase(PROP_NUM_THREADS)) {
      return Integer.toString(numThreads);
    } else if (propertyName.equalsIgnoreCase(PROP_USE_GPU)) {
        return Boolean.toString(useGPU);
    } else if (propertyName.equalsIgnoreCase(PROP_GPU_LIBRARY_PATH)) {
        return gpuLibraryPath;
    } else {
      // try the superclass:
      return super.getProperty(propertyName);
    }
  }

  /* (non-Javadoc)
   * @see infodynamics.measures.continuous.MutualInfoMultiVariateCommon#finaliseAddObservations()
   */
  @Override
  public void finaliseAddObservations() throws Exception {
    // Allow the parent to generate the data for us first
    super.finaliseAddObservations();
    
    if (totalObservations <= k + 2*dynCorrExclTime) {
      throw new Exception("There are less observations provided (" +
          totalObservations +
          ") than required for the number of nearest neighbours parameter (" +
          k + ") and any dynamic correlation exclusion (" + dynCorrExclTime + ")");
    }
    
    // Set the constants:
    digammaK = MathsUtils.digamma(k);
    digammaN = MathsUtils.digamma(totalObservations);
  }

  /**
   * {@inheritDoc} 
   * 
   * @return the average MI in nats (not bits!)
   */
  public double computeAverageLocalOfObservations() throws Exception {
    // Compute the MI
    double startTime = Calendar.getInstance().getTimeInMillis();
    lastAverage = computeFromObservations(false, null)[0];
    miComputed = true;
    if (debug) {
      Calendar rightNow2 = Calendar.getInstance();
      long endTime = rightNow2.getTimeInMillis();
      System.out.println("Calculation time: " + ((endTime - startTime)/1000.0) + " sec" );
    }
    return lastAverage;
  }

  /**
   * {@inheritDoc} 
   * 
   * @return the MI under the new ordering, in nats (not bits!).
   */
  public double computeAverageLocalOfObservations(int[] reordering) throws Exception {
    if (reordering == null) {
      return computeAverageLocalOfObservations();
    }
    
    // Save internal variables pertaining to the original order for
    //  later reinstatement:
    // (don't need to save and reinstate kdTreeSource, since we're not
    //  altering the source data order).
    KdTree originalKdTreeJoint = kdTreeJoint;
    kdTreeJoint = null; // So that it is rebuilt for the new ordering
    NearestNeighbourSearcher originalKdTreeDest = nnSearcherDest;
    nnSearcherDest = null; // So that it is rebuilt for the new ordering
    double[][] originalData2 = destObservations;
    
    // Generate a new re-ordered data2
    destObservations = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData2, reordering);
    // Compute the MI
    double newMI = computeFromObservations(false, null)[0];
    
    // restore original variables:
    destObservations = originalData2;
    kdTreeJoint = originalKdTreeJoint;
    nnSearcherDest = originalKdTreeDest;
    
    return newMI;
  }

  /**
   * <p>Computes the local values of the MI,
   *  for each valid observation in the previously supplied observations
   *  (with PDFs computed using all of the previously supplied observation sets).</p>
   *  
   * <p>If the samples were supplied via a single call such as
   * {@link #setObservations(double[][], double[][])},
   * then the return value is a single time-series of local
   * channel measure values corresponding to these samples.</p>
   * 
   * <p>Otherwise where disjoint time-series observations were supplied using several 
   *  calls such as {@link #addObservations(double[][], double[][])}
   *  then the local values for each disjoint observation set will be appended here
   *  to create a single "time-series" return array.</p>
   * 
   * @return the "time-series" of local MIs in nats
   * @throws Exception
   */
  public double[] computeLocalOfPreviousObservations() throws Exception {
    double[] localValues = computeFromObservations(true, null);
    lastAverage = MatrixUtils.mean(localValues);
    miComputed = true;
    return localValues;
  }

  @Override
  public double[] computeLocalUsingPreviousObservations(double[][] states1, double[][] states2) throws Exception {
	// Do normalisation of the incoming data if required:
	double[][] states1ToUse, states2ToUse;
	if (normalise) {
		states1ToUse = MatrixUtils.normaliseIntoNewArray(states1, sourceMeansBeforeNorm, sourceStdsBeforeNorm, 0, states1.length-timeDiff);
		states2ToUse = MatrixUtils.normaliseIntoNewArray(states2, destMeansBeforeNorm, destStdsBeforeNorm, timeDiff, states2.length-timeDiff);
	} else {
		if (timeDiff > 0) {
			states1ToUse = MatrixUtils.selectRows(states1, 0, states1.length-timeDiff);
			states2ToUse = MatrixUtils.selectRows(states2, 0, states2.length-timeDiff);
		} else {
			states1ToUse = states1;
			states2ToUse = states2;
		}
	}
	// And call the algorithm:
	double[] localValues = computeFromObservations(true,
			new double[][][]{states1ToUse, states2ToUse});
	return localValues;
  }
  
  /**
   * This protected method handles the multiple threads which
   *  computes either the average or local MI (over parts of the total
   *  observations), computing the x and y 
   *  distances between all tuples in time.
   * 
   * <p>The method returns:<ol>
   *  <li>for (returnLocals == false), an array of size 1,
   *      containing the average MI </li>
   *  <li>for local MIs (returnLocals == true), the array of local MI values</li>
   *  </ol>
   * 
   * @param returnLocals whether to return an array or local values, or else
   *  sums of these values
   * @param newObservations set to null for computing for the observation set for the PDF, or pass in a new set
   *  of observations to compute the average/locals for (using the existing observations to construct the PDF)
   * @return either the average MI, or array of local MI value, in nats not bits
   * @throws Exception
   */
  protected double[] computeFromObservations(boolean returnLocals, double[][][] newObservations) throws Exception {
    int N = sourceObservations.length; // number of observations for the PDFs
    
    double[] returnValues = null;
    
    // How many time points are we averaging over?
    int numTimePointsToComputeFor = (newObservations == null) ?
        N : newObservations[0].length;
    
    if (useGPU && (newObservations != null)) {
    	System.out.println("Cannot use GPU for estimation based on new observations -- falling back to CPU calculation...");
    }
    
    if (useGPU && (newObservations == null)) {
      returnValues = gpuComputeFromObservations(0, N, returnLocals);
    } else if (numThreads == 1) {
      // Single-threaded implementation:
      ensureKdTreesConstructed();

      if (newObservations == null) {
          returnValues = partialComputeFromObservations(0, numTimePointsToComputeFor, returnLocals);
        } else {
          returnValues = partialComputeFromNewObservations(0, numTimePointsToComputeFor,
              newObservations[0], newObservations[1], returnLocals);
        }

    } else {
      // We're going multithreaded:
      ensureKdTreesConstructed();

      if (returnLocals) {
        // We're computing local MI
        returnValues = new double[numTimePointsToComputeFor];
      } else {
        // We're computing average MI
        returnValues = new double[MiKraskovThreadRunner.RETURN_ARRAY_LENGTH_MI];
      }
      
      // Distribute the observations to the threads for the parallel processing
      int lTimesteps = numTimePointsToComputeFor / numThreads; // each thread gets the same amount of data
      int res = numTimePointsToComputeFor % numThreads; // the first thread gets the residual data
      if (debug) {
        System.out.printf("Computing Kraskov MI with %d threads (%d timesteps each, plus %d residual)\n",
            numThreads, lTimesteps, res);
      }
      Thread[] tCalculators = new Thread[numThreads];
      MiKraskovThreadRunner[] runners = new MiKraskovThreadRunner[numThreads];
      for (int t = 0; t < numThreads; t++) {
        int startTime = (t == 0) ? 0 : lTimesteps * t + res;
        int numTimesteps = (t == 0) ? lTimesteps + res : lTimesteps;
        if (debug) {
          System.out.println(t + ".Thread: from " + startTime +
              " to " + (startTime + numTimesteps)); // Trace Message
        }
        runners[t] = new MiKraskovThreadRunner(this, startTime, numTimesteps, newObservations, returnLocals);
        tCalculators[t] = new Thread(runners[t]);
        tCalculators[t].start();
      }
      
      // Here, we should wait for the termination of the all threads
      //  and collect their results
      for (int t = 0; t < numThreads; t++) {
        if (tCalculators[t] != null) { // TODO Ipek: can you comment on why we're checking for null here?
          tCalculators[t].join(); 
        }
        // Now we add in the data from this completed thread:
        if (returnLocals) {
          // We're computing local MI; copy these local values
          //  into the full array of locals
          System.arraycopy(runners[t].getReturnValues(), 0, 
              returnValues, runners[t].myStartTimePoint, runners[t].numberOfTimePoints);
        } else {
          // We're computing the average MI, keep the running sums of digammas and counts
          MatrixUtils.addInPlace(returnValues, runners[t].getReturnValues(), MiKraskovThreadRunner.RETURN_ARRAY_LENGTH_MI);
        }
      }
    }
    
    // Finalise the results:
    if (returnLocals) {
      return returnValues;
    } else {
      // Compute the average number of points within eps_x and eps_y
      double averageDiGammas = returnValues[MiKraskovThreadRunner.INDEX_SUM_DIGAMMAS] / (double) numTimePointsToComputeFor;
      double avNx = returnValues[MiKraskovThreadRunner.INDEX_SUM_NX] / (double) numTimePointsToComputeFor;
      double avNy = returnValues[MiKraskovThreadRunner.INDEX_SUM_NY] / (double) numTimePointsToComputeFor;
      if (debug) {
        System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
      }
      
      // Use digamma(N) normally, unless we're looking at new observations:
      double digammaNToUse = (newObservations == null) ? digammaN : MathsUtils.digamma(totalObservations+1);
      
      // Finalise the average result, depending on which algorithm we are implementing:
      if (isAlgorithm1) {
        return new double[] { digammaK - averageDiGammas + digammaNToUse };
      } else {
        return new double[] { digammaK - (1.0 / (double)k) - averageDiGammas + digammaNToUse };
      }
    }
  }
  
  /**
   * Use the same nearest neighbour searches to return a KL-based 
   * conditional entropy, sharing the neighbour search radius.
   * With reference to the KSG paper, this takes the KL-based estimator for H(X) (eqn 22) and adds it to the negative of MI.
   * It should also be noted that if the variables are normalised, this changes the absolute values of the entropies,
   *  and they become relative to normalised variables.
   * Since this is a side-method for this estimator, it is currently is only implemented single threaded on CPU.
   * It should also be considered experimental - this will eventually be moved to its own estimator.
   * 
   * @return the average conditional entropy of variable1 given variable2 in nats (not bits!)
   */
  public double computeAverageConditionalEntropy() throws Exception {
    // Compute the MI
    double startClockTime = Calendar.getInstance().getTimeInMillis();
    double conditionalEntropy = 0.0;
    conditionalEntropy = computeFromObservations(false, null)[0];
    
    int N = sourceObservations.length; // number of observations for the PDFs
    
    // Single-threaded implementation:
    ensureKdTreesConstructed();
    double avDigammaCond = 0;
    double avLog2xkNNDist = 0;
    for (int t = 0; t < N; t++) {
    	// For this sample only:
    	double[] returnValues = partialComputeFromObservations(t, 1, false);
    	int n_y = (int) returnValues[MiKraskovThreadRunner.INDEX_SUM_NY];
    	if (isAlgorithm1) {
    		avDigammaCond += MathsUtils.digamma(n_y+1);
    	} else {
    		avDigammaCond += MathsUtils.digamma(n_y);
    	}
    	double doublekNNDist = returnValues[MiKraskovThreadRunner.INDEX_SUM_2X_NEIGH_DIST];
    	avLog2xkNNDist += Math.log(doublekNNDist);
    }
    avDigammaCond /= (double) N;
    avLog2xkNNDist *= (double) dimensionsSource / (double) N;
    if (isAlgorithm1) {
    	conditionalEntropy = -digammaK + avDigammaCond + avLog2xkNNDist;
    } else {
    	conditionalEntropy = -digammaK + (1.0 / (double)k) + avDigammaCond + avLog2xkNNDist;
    }
    
    if (debug) {
      Calendar rightNow2 = Calendar.getInstance();
      long endTime = rightNow2.getTimeInMillis();
      System.out.println("Calculation time: " + ((endTime - startClockTime)/1000.0) + " sec" );
    }
    return conditionalEntropy;
  }

  /**
   * Protected method to be used internally for threaded implementations.
   * This method implements the guts of each Kraskov algorithm, computing the number of 
   *  nearest neighbours in each dimension for a sub-set of the data points.
   *  It is intended to be called by one thread to work on that specific
   *  sub-set of the data.
   * 
   * <p>The method returns:<ol>
   *  <li>for average MIs (returnLocals == false), the relevant sums of digamma(n_x+1) and digamma(n_y+1)
   *     for a partial set of the observations</li>
   *  <li>for local MIs (returnLocals == true), the array of local MI values</li>
   *  </ol>
   * 
   * @param startTimePoint start time for the partial set we examine
   * @param numTimePoints number of time points (including startTimePoint to examine)
   * @param returnLocals whether to return an array or local values, or else
   *  sums of these values
   * @return an array of sum of digamma(n_x+1) and digamma(n_y+1), then
   *  sum of n_x and finally sum of n_y (these latter two are for debugging purposes).
   * @throws Exception
   */
  protected abstract double[] partialComputeFromObservations(
      int startTimePoint, int numTimePoints, boolean returnLocals) throws Exception;

	/**
	 * Protected method to be used internally for threaded implementations.
	 * This method implements the guts of each Kraskov algorithm, computing the number of 
	 *  nearest neighbours in each dimension for a sub-set of the data points.
	 *  It is intended to be called by one thread to work on that specific
	 *  sub-set of the data.
	 * In particular, this method differs from {@link #partialComputeFromObservations(int, int, boolean)}
	 *  because it operates on a new set of observations (using the old set of observations for
	 *  constructing the search spaces and PDFs)
	 * 
	 * <p>The method returns:<ol>
	 *  <li>for average MIs (returnLocals == false), the relevant sums of digamma(n_x+1), digamma(n_y+1)
	 *     for a partial set of the observations</li>
	 *  <li>for local MIs (returnLocals == true), the array of local MI values</li>
	 *  </ol>
	 * 
	 * @param startTimePoint start time for the partial set we examine
	 * @param numTimePoints number of time points (including startTimePoint to examine)
	 * @param newVar1Observations new time series of observations for variable 1
	 * @param newVar2Observations new time series of observations for variable 2
	 * @param returnLocals whether to return an array or local values, or else
	 *  sums of these values
     * @return an array of sum of digamma(n_x+1) and digamma(n_y+1), then
     *  sum of n_x and finally sum of n_y (these latter two are for debugging purposes).
	 * @throws Exception
	 */
	protected abstract double[] partialComputeFromNewObservations(
			int startTimePoint, int numTimePoints,
			double[][] newVar1Observations, double[][] newVar2Observations, 
			boolean returnLocals) throws Exception;

  /**
   * Protected method to be used internally for GPU implementations.
   * This method serves the same purpose as partialComputeFromObservations,
   *   but for GPU computation. Each algorithm must override this method
   *   and implement a GPU routine to calculate all values in a single call
   *   to the GPU code.
   */
  protected double[] gpuComputeFromObservations(int startTimePoint,
      int numTimePoints, boolean returnLocals, int nb_surrogates,
      int[][] newOrderings) throws Exception {

    if (debug) {
      System.out.println("Start GPU calculation");
    }

    ensureCudaLibraryLoaded();

    boolean useMaxNorm;

    if ( normType == EuclideanUtils.NORM_MAX_NORM) {
      useMaxNorm = true;
    } else if ( normType == EuclideanUtils.NORM_EUCLIDEAN || normType == EuclideanUtils.NORM_EUCLIDEAN_SQUARED) {
      useMaxNorm = false;
    } else {
      throw new Exception("Only max and square norms are implemented. Abort.");
    }

    double[] res;

    try {
      if (debug) {
        System.out.printf("Calling GPU calculation with returnLocals=%b and nb_surrogates=%d\n", returnLocals, nb_surrogates);
      }
      res = MIKraskov(totalObservations, sourceObservations, dimensionsSource,
          destObservations, dimensionsDest, k, dynCorrExclTime, returnLocals, useMaxNorm,
          isAlgorithm1, nb_surrogates, null!=newOrderings, newOrderings);
      if (debug) {
        System.out.println("GPU calculation finished successfully. Returning results");
      }
    } catch (Throwable e) {
      System.out.println("WARNING. Error in GPU code. Reverting back to CPU.");
      e.printStackTrace();
      res = partialComputeFromObservations(0, totalObservations, returnLocals);
    }

    return res;
  }


  /**
   * FIXME
   */
  protected double[] gpuComputeFromObservations(int startTimePoint,
      int numTimePoints, boolean returnLocals) throws Exception {
    return gpuComputeFromObservations(startTimePoint, numTimePoints, returnLocals, 0, null);
  }

  /**
   * FIXME
   */
  protected double[] gpuComputeFromObservations(int startTimePoint,
      int numTimePoints, boolean returnLocals, int[][] newOrderings) throws Exception {
    return gpuComputeFromObservations(startTimePoint, numTimePoints,
        returnLocals, newOrderings.length, newOrderings);
  }

  /**
   * Native method to calculate MI in GPU.
   */
  private native double[] MIKraskov(
      int N, double[][] source, int dimx, double[][] dest, int dimy,
      int k, int theiler, boolean returnLocals, boolean useMaxNorm,
      boolean isAlgorithm1, int nb_surrogates, boolean orderingsGiven,
      int[][] newOrderings);

  /**
   * Internal method to ensure that the Kd-tree data structures to represent the
   * observational data have been constructed (should be called prior to attempting
   * to use these data structures)
   */
  protected void ensureKdTreesConstructed() throws Exception {
    
    // We need to construct the k-d trees for use by the child
    //  classes. We check each tree for existence separately
    //  since source can be used across original and surrogate data
    // TODO can parallelise these -- best done within the kdTree --
    //  though it's unclear if there's much point given that
    //  the tree construction itself afterwards can't really be well parallelised.
    if (kdTreeJoint == null) {
      kdTreeJoint = new KdTree(new int[] {dimensionsSource, dimensionsDest},
            new double[][][] {sourceObservations, destObservations},
            observationSetIndices, observationTimePoints);
      kdTreeJoint.setNormType(normType);
    }
    if (nnSearcherSource == null) {
      nnSearcherSource = NearestNeighbourSearcher.create(sourceObservations,
    		  observationSetIndices, observationTimePoints);
      nnSearcherSource.setNormType(normType);
    }
    if (nnSearcherDest == null) {
      nnSearcherDest = NearestNeighbourSearcher.create(destObservations,
    		  observationSetIndices, observationTimePoints);
      nnSearcherDest.setNormType(normType);
    }
  }

  /**
   * Internal method to ensure that the Cuda native library has been correctly
   * loaded.
   *
   * @throws Exception if library not found or unable to load
   */
  protected void ensureCudaLibraryLoaded() throws Exception {

    if (!cudaLibraryLoaded) {

      try {
        if (gpuLibraryPath.length() < 1) {
          NativeUtils.loadLibraryFromJar("/cuda/libKraskov.so");
        } else {
          System.load(gpuLibraryPath);
        }
      } catch (Throwable e) {
        String errmsg = "GPU library not found. To compile GPU code set the enablegpu flag to true in build.xml, or run `ant gpu jar`.";
        errmsg += "\nFor more information see the JIDT GPU wiki page: https://github.com/jlizier/jidt/wiki/GPU";
        if (gpuLibraryPath.length() > 0) {
          errmsg += "\n\nGPU library was not found in the path provided. Provide full path including library file name.";
          errmsg += "\nExample: /home/johndoe/myfolder/libKraskov.so";
        }
        throw new Exception(errmsg);
      }

      cudaLibraryLoaded = true;

    }
  }

  /**
   * {@inheritDoc} 
   */
  @Override
  public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck)
      throws Exception {
    if (useGPU) {
      double[] res = gpuComputeFromObservations(0, totalObservations, false, numPermutationsToCheck, null);
      return new EmpiricalMeasurementDistribution(
          MatrixUtils.select(res, 1, res.length - 1), res[0]);
    } else {
      return super.computeSignificance(numPermutationsToCheck);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings)
      throws Exception {
    if (useGPU) {
      double[] res = gpuComputeFromObservations(0, totalObservations, false, newOrderings.length, newOrderings);
      return new EmpiricalMeasurementDistribution(
          MatrixUtils.select(res, 1, res.length - 1), res[0]);
    } else {
      return super.computeSignificance(newOrderings);
    }
  }
  
  /**
   * Compute the prediction error in one variable from the k nearest neighbours (kNNs) of the observation
   *  of the other variable.
   * The kNNs of the variable to predict from are formed from the supplied norm type within that variable.
   * The number of kNNs to use here is the current property value set for {@link #PROP_K}.
   * The prediction error is a sum of absolute errors for each dimension within the variable to predict
   * 
   * @param predictFirstVariable true for predicting the first variable (Source) or
   * false for predicting the second variable (destination)
   * @return array of prediction errors for each dimension of the predicted variable
   * @throws Exception
   */
  public double[] computePredictionErrorsFromObservations(boolean predictFirstVariable) throws Exception {
    return computePredictionErrorsFromObservations(predictFirstVariable, k);
  }
  
  /**
   * Compute the prediction error in one variable from the k nearest neighbours (kNNs) of the observation
   *  of the other variable.
   * The kNNs of the variable to predict from are formed from the supplied norm type within that variable.
   * The prediction error is a sum of absolute errors for each dimension within the variable to predict
   * 
   * @param predictFirstVariable true for predicting the first variable (Source) or
   * false for predicting the second variable (destination)
   * @param kNNs number of nearest neighbours to use
   * @return array of prediction errors for each dimension of the predicted variable
   * @throws Exception
   */
  public double[] computePredictionErrorsFromObservations(boolean predictFirstVariable, int kNNs) throws Exception {
    int N = sourceObservations.length; // number of observations
    
    double[] totalErrors = null;
    
    ensureKdTreesConstructed();
    
    if (numThreads == 1) {
      // Single-threaded implementation:
      totalErrors = partialComputePredictionErrorFromObservations(0, N, kNNs, predictFirstVariable);
    } else {
      // We're going multithreaded:
      totalErrors = new double[predictFirstVariable ? dimensionsSource : dimensionsDest];
      
      // Distribute the observations to the threads for the parallel processing
      int lTimesteps = N / numThreads; // each thread gets the same amount of data
      int res = N % numThreads; // the first thread gets the residual data
      if (debug) {
        System.out.printf("Computing prediction errors for variable %d from variable %d with %d threads (%d timesteps each, plus %d residual)\n",
            predictFirstVariable ? 1 : 2, predictFirstVariable ? 2 : 1,
            numThreads, lTimesteps, res);
      }
      Thread[] tCalculators = new Thread[numThreads];
      MiKraskovPredictionThreadRunner[] runners = new MiKraskovPredictionThreadRunner[numThreads];
      for (int t = 0; t < numThreads; t++) {
        int startTime = (t == 0) ? 0 : lTimesteps * t + res;
        int numTimesteps = (t == 0) ? lTimesteps + res : lTimesteps;
        if (debug) {
          System.out.println(t + ".Thread: from " + startTime +
              " to " + (startTime + numTimesteps)); // Trace Message
        }
        runners[t] = new MiKraskovPredictionThreadRunner(this, startTime,
            numTimesteps, kNNs, predictFirstVariable);
        tCalculators[t] = new Thread(runners[t]);
        tCalculators[t].start();
      }
      
      // Here, we should wait for the termination of the all threads
      //  and collect their results
      for (int t = 0; t < numThreads; t++) {
        if (tCalculators[t] != null) {
          tCalculators[t].join();
        }
        // Now we add in the data from this completed thread:
        MatrixUtils.addInPlace(totalErrors, runners[t].getReturnValues());
      }
    }
    
    // Finalise the results:
    if (debug) {
      System.out.printf("Total prediction error from variable %d to variable %d=",
          predictFirstVariable ? 2 : 1, predictFirstVariable ? 1 : 2);
      MatrixUtils.printArray(System.out, 3, totalErrors);
    }
    return totalErrors;
  }
  
  /**
   * Private class to handle multi-threading of the Kraskov algorithms.
   * Each instance calls partialComputeFromObservations()
   * to compute nearest neighbours for a part of the data.
   * 
   * 
   * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
   * <a href="http://lizier.me/joseph/">www</a>)
   * @author Ipek Özdemir
   */
  private class MiKraskovThreadRunner implements Runnable {
    protected MutualInfoCalculatorMultiVariateKraskov miCalc;
    protected int myStartTimePoint;
    protected int numberOfTimePoints;
	protected double[][][] newObservations;
    protected boolean computeLocals;
    
    protected double[] returnValues = null;
    protected Exception problem = null;
    
    public static final int INDEX_SUM_DIGAMMAS = 0;
    public static final int INDEX_SUM_NX = 1;
    public static final int INDEX_SUM_NY = 2;
    public static final int INDEX_SUM_2X_NEIGH_DIST = 3; // Not used by GPU, and only read for conditional entropy
    public static final int RETURN_ARRAY_LENGTH_MI = 3;
    public static final int RETURN_ARRAY_LENGTH = 4;

    public MiKraskovThreadRunner(
        MutualInfoCalculatorMultiVariateKraskov miCalc,
        int myStartTimePoint, int numberOfTimePoints,
		double[][][] newObservations,
        boolean computeLocals) {
      this.miCalc = miCalc;
      this.myStartTimePoint = myStartTimePoint;
      this.numberOfTimePoints = numberOfTimePoints;
      this.computeLocals = computeLocals;
      this.newObservations = newObservations;
    }
    
    /**
     * Return the values from this part of the data,
     *  or throw any exception that was encountered by the 
     *  thread.
     * 
     * @return the relevant return values from this part of the data
     * @throws Exception an exception previously encountered by this thread.
     */
    public double[] getReturnValues() throws Exception {
      if (problem != null) {
        throw problem;
      }
      return returnValues;
    }
    
    /**
     * Start the thread for the given parameters
     */
    public void run() {
      try {
			if (newObservations == null) {
				// Computing on existing observations
				returnValues = miCalc.partialComputeFromObservations(
						myStartTimePoint, numberOfTimePoints, computeLocals);
			} else {
				// Computing on new observations
				returnValues = miCalc.partialComputeFromNewObservations(
						myStartTimePoint, numberOfTimePoints,
						newObservations[0], newObservations[1],
						computeLocals);
			}
      } catch (Exception e) {
        // Store the exception for later retrieval
        problem = e;
        return;
      }
    }
  }
  // end class MiKraskovThreadRunner  
  
  /**
   * Private class to handle multi-threading of the prediction from
   * k nearest neighbours.
   * Each instance calls partialComputePredictionErrorFromObservations()
   * to compute nearest neighbours for a part of the data.
   * 
   * 
   * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
   * <a href="http://lizier.me/joseph/">www</a>)
   * @author Ipek Özdemir
   */
  private class MiKraskovPredictionThreadRunner implements Runnable {
    protected MutualInfoCalculatorMultiVariateKraskov miCalc;
    protected int myStartTimePoint;
    protected int numberOfTimePoints;
    protected int kNNs;
    protected boolean predictFirstVariable;
    
    protected double[] returnValues = null;
    protected Exception problem = null;
    
    public MiKraskovPredictionThreadRunner(
        MutualInfoCalculatorMultiVariateKraskov miCalc,
        int myStartTimePoint, int numberOfTimePoints,
        int kNNs, boolean predictFirstVariable) {
      this.miCalc = miCalc;
      this.myStartTimePoint = myStartTimePoint;
      this.numberOfTimePoints = numberOfTimePoints;
      this.kNNs = kNNs;
      this.predictFirstVariable = predictFirstVariable;
    }
    
    /**
     * Return the sum of prediction errors from this part of the data,
     *  or throw any exception that was encountered by the 
     *  thread.
     * 
     * @return sum of prediction errors from this part of the data, for
     * each of the dimensions of the relevant variable
     * @throws Exception an exception previously encountered by this thread.
     */
    public double[] getReturnValues() throws Exception {
      if (problem != null) {
        throw problem;
      }
      return returnValues;
    }
    
    /**
     * Start the thread for the given parameters
     */
    public void run() {
      try {
        returnValues = miCalc.partialComputePredictionErrorFromObservations(
            myStartTimePoint, numberOfTimePoints, kNNs, predictFirstVariable);
      } catch (Exception e) {
        // Store the exception for later retrieval
        problem = e;
        return;
      }
    }
  }
  // end class MiKraskovPredictionThreadRunner

  /**
   * Protected method to be used internally for threaded implementations.
   * This method implements the guts of examining prediction errors in one variable
   *  from the k nearest neighbours of the other.
   *  It is intended to be called by one thread to work on that specific
   *  sub-set of the data.
   * 
   * @param startTimePoint start time for the partial set we examine
   * @param numTimePoints number of time points (including startTimePoint to examine)
   * @param kNNs number of nearest neighbours to use
   * @param predictFirstVariable whether to use the second variable to predict the first (true)
   *    or first variable to predict the second (false)
   * @return an array of the sum of square prediction errors for each dimension within the predicted
   *    variable.
   * @throws Exception
   */
  protected double[] partialComputePredictionErrorFromObservations(
      int startTimePoint, int numTimePoints, int kNNs, boolean predictFirstVariable) throws Exception {
    
    double startTime = Calendar.getInstance().getTimeInMillis();

    double[] totalErrors = new double[predictFirstVariable ? dimensionsSource : dimensionsDest];
    
    for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
      // Find the k nearest neighbours for the relevant predictor variable
      if (predictFirstVariable) {
        // First variable value to predict:
        double[] sourceValueToPredict = sourceObservations[t];
        // Find kNNs of second variable
        PriorityQueue<NeighbourNodeData> nnPQ = 
            nnSearcherDest.findKNearestNeighbours(kNNs, t, dynCorrExclTime);
        double[] predictedValue = new double[dimensionsSource];
        for (NeighbourNodeData kthNnData : nnPQ) {
          // Retrieve the source value this corresponds to
          double[] neighbourSourceValue = sourceObservations[kthNnData.sampleIndex];
          // And include it's contribution in the prediction
          for (int d = 0; d < dimensionsSource; d++) {
            predictedValue[d] += neighbourSourceValue[d];
          }
        }
        // Now add in the square prediction errors from the prediction: 
        for (int d = 0; d < dimensionsSource; d++) {
          predictedValue[d] /= (double) kNNs;
          totalErrors[d] += (sourceValueToPredict[d] - predictedValue[d]) * 
                    (sourceValueToPredict[d] - predictedValue[d]);
        }
      } else { // predict second variable
        // Second variable value to predict:
        double[] destValueToPredict = destObservations[t];
        // Find kNNs of first variable
        PriorityQueue<NeighbourNodeData> nnPQ = 
            nnSearcherSource.findKNearestNeighbours(kNNs, t, dynCorrExclTime);
        double[] predictedValue = new double[dimensionsDest];
        for (NeighbourNodeData kthNnData : nnPQ) {
          // Retrieve the dest value this corresponds to
          double[] neighbourDestValue = destObservations[kthNnData.sampleIndex];
          // And include it's contribution in the prediction
          for (int d = 0; d < dimensionsDest; d++) {
            predictedValue[d] += neighbourDestValue[d];
          }
        }
        // Now add in the square prediction errors from the prediction: 
        for (int d = 0; d < dimensionsDest; d++) {
          predictedValue[d] /= (double) kNNs;
          totalErrors[d] += (destValueToPredict[d] - predictedValue[d]) * 
                    (destValueToPredict[d] - predictedValue[d]);
        }
      }
    }
    
    if (debug) {
      Calendar rightNow2 = Calendar.getInstance();
      long endTime = rightNow2.getTimeInMillis();
      System.out.println("Subset " + startTimePoint + ":" +
          (startTimePoint + numTimePoints) + " Calculation time: " +
          ((endTime - startTime)/1000.0) + " sec" );
    }
    
    return totalErrors;
  } 
}
