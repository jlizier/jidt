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

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.KdTree;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.NearestNeighbourSearcher;
import infodynamics.utils.NeighbourNodeData;
import infodynamics.utils.UnivariateNearestNeighbourSearcher;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.NativeUtils;

/**
 * <p>Computes the differential conditional mutual information of two multivariate
 *  <code>double[][]</code> sets of observations, conditioned on another
 *  (implementing {@link ConditionalMutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see references below).
 *  The implementation is made using fast-neighbour searches with an
 *  underlying k-d tree algorithm.
 *  This is an abstract class, building on the common code base in
 *  {@link ConditionalMutualInfoMultiVariateCommon},
 *  to gather common functionality between the two
 *  algorithms defined by Kraskov et al.
 *  Two child classes {@link ConditionalMutualInfoCalculatorMultiVariateKraskov1} and
 *  {@link ConditionalMutualInfoCalculatorMultiVariateKraskov2} then
 *  actually implement the two KSG algorithms.</p>
 * 
 * <p>Crucially, the calculation is performed by examining
 * neighbours in the full joint space (as specified by Frenzel and Pompe,
 * and later by Vejmelka and Paluš, and Vlachos and Kugiumtzis)
 * rather than two MI calculators.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link ConditionalMutualInfoCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>For constructors see the child classes.</li>
 *  <li>Further properties are defined in {@link #setProperty(String, String)}.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 * 
 * <p>Finally, note that {@link Cloneable} is implemented allowing clone()
 *  to produce only an automatic shallow copy, which is fine
 *  for the statistical significance calculation it is intended for
 *  (none of the array 
 *  data will be changed there).</p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Frenzel and Pompe, <a href="http://dx.doi.org/10.1103/physrevlett.99.204101">
 * 	"Partial Mutual Information for Coupling Analysis of Multivariate Time Series"</a>,
 * 	Physical Review Letters, <b>99</b>, p. 204101+ (2007).</li>
 *  <li>Vejmelka and Paluš, <a href="http://dx.doi.org/10.1103/physreve.77.026214">
 *  "Inferring the directionality of coupling with conditional mutual information"</a>,
 *  Physical Review E, <b>77</b>, 026214, (2008)</li>
 *  <li>I. Vlachos and D. Kugiumtzis, <a href="http://dx.doi.org/10.1103/physreve.82.016207">
 *  "Nonuniform state-space reconstruction and coupling detection"</a>,
 *  Physical Review E, <b>82</b>, 016207, (2010)</li>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @author Ipek Özdemir
 */
public abstract class ConditionalMutualInfoCalculatorMultiVariateKraskov 
	extends ConditionalMutualInfoMultiVariateCommon
	implements Cloneable { // See comments on clonability above
	
	/**
	 * we compute distances to the kth neighbour in the joint space
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
	 * Property name for a dynamics exclusion time window 
	 * otherwise known as Theiler window (see Kantz and Schreiber).
	 * Default is 0 which means no dynamic exclusion window.
	 */
	public static final String PROP_DYN_CORR_EXCL_TIME = "DYN_CORR_EXCL";	
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
	 * Whether we use dynamic correlation exclusion
	 */
	protected boolean dynCorrExcl = false;
	/**
	 * Size of dynamic correlation exclusion window.
	 */
	protected int dynCorrExclTime = 0;
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
	 * protected k-d tree data structure (for fast nearest neighbour searches)
	 *  representing the (var1,conditional) space
	 */
	protected KdTree kdTreeVar1Conditional;
	/**
	 * protected univariate neighbour searcher data structure (for fast nearest neighbour searches)
	 *  representing the (var1) space; used only if var1 is univariate
	 */
	protected UnivariateNearestNeighbourSearcher uniNNSearcherVar1;
	/**
	 * protected k-d tree data structure (for fast nearest neighbour searches)
	 *  representing the (var2,conditional) space
	 */
	protected KdTree kdTreeVar2Conditional;
	/**
	 * protected univariate neighbour searcher data structure (for fast nearest neighbour searches)
	 *  representing the (var2) space; used only if var2 is univariate
	 */
	protected UnivariateNearestNeighbourSearcher uniNNSearcherVar2;
	/**
	 * protected data structure (for fast nearest neighbour searches)
	 *  representing the conditional space.
	 *  Could be implemented by either a k-d tree or sorted array, depending
	 *  on whether the conditional is multi-variate or not.
	 */
	protected NearestNeighbourSearcher nnSearcherConditional;
	
	/**
	 * Constant for digamma(k), with k the number of nearest neighbours selected
	 */
	protected double digammaK;
	
	/**
	 * Construct an instance of the KSG conditional MI calculator
	 */
	public ConditionalMutualInfoCalculatorMultiVariateKraskov() {
		super();
		addNoise = true;
		noiseLevel = 1e-8; // to match the noise order in MILCA toolkit.
	}

	@Override
	public void initialise(int dimensions1, int dimensions2, int dimensionsCond) {
		kdTreeJoint = null;
		kdTreeVar1Conditional = null;
		kdTreeVar2Conditional = null;
		nnSearcherConditional = null;
		uniNNSearcherVar1 = null;
		uniNNSearcherVar2 = null;
		super.initialise(dimensions1, dimensions2, dimensionsCond);
	}

	/**
	 * Sets properties for the KSG conditional MI calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 *  <li>{@link #PROP_K} -- number of k nearest neighbours to use in joint kernel space
	 *      in the KSG algorithm (default is 4).</li>
	 *  <li>{@link #PROP_NORM_TYPE}</li> -- normalization type to apply to 
	 *      working out the norms between the points in each marginal space.
	 *      Options are defined by {@link KdTree#setNormType(String)} -
	 *      default is {@link EuclideanUtils#NORM_MAX_NORM_STRING}.
	 *  <li>{@link #DYN_CORR_EXCL_TIME_NAME} -- a dynamics exclusion time window,
	 *      also known as Theiler window (see Kantz and Schreiber);
	 *      default is 0 which means no dynamic exclusion window.</li>
	 *  <li>{@link #PROP_NUM_THREADS} -- the integer number of parallel threads
	 *  	to use in the computation. Can be passed as a string "USE_ALL"
	 *      to use all available processors on the machine.
	 *      Default is "USE_ALL".
	 *  <li>any valid properties for {@link ConditionalMutualInfoMultiVariateCommon#setProperty(String, String)},
	 *     notably including {@link ConditionalMutualInfoMultiVariateCommon#PROP_NORMALISE}.</li>
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
			normType = KdTree.validateNormType(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_DYN_CORR_EXCL_TIME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
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
			// Assume this is a property for the common parent class
			super.setProperty(propertyName, propertyValue);
		}
	}

	@Override
	public String getProperty(String propertyName) {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			return Integer.toString(k);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			return KdTree.convertNormTypeToString(normType);
		} else if (propertyName.equalsIgnoreCase(PROP_DYN_CORR_EXCL_TIME)) {
			return Integer.toString(dynCorrExclTime);
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
	 * @see infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon#finaliseAddObservations()
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
	}

	/**
	 * {@inheritDoc} 
	 * 
	 * @return the average conditional MI in nats (not bits!)
	 */
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		// Compute the conditional MI
		double startTime = Calendar.getInstance().getTimeInMillis();
		lastAverage = computeFromObservations(false, null)[0];
		condMiComputed = true;
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
	 * If {@code reordering} is null, it is assumed there is no reordering of
	 *  the given variable.
	 *  
	 * @return the conditional MI under the new ordering, in nats (not bits!).
	 */
	@Override
	public double computeAverageLocalOfObservations(int variableToReorder, 
			int[] reordering) throws Exception {
		
		if (reordering == null) {
			return computeAverageLocalOfObservations();
		}
		double[][] originalData;
		KdTree originalKdTreeJoint = kdTreeJoint;
		kdTreeJoint = null; // So that it is rebuilt for the new ordering
		KdTree originalKdTreeVar1Conditional = kdTreeVar1Conditional;
		UnivariateNearestNeighbourSearcher originalUniNNSearcherVar1 = uniNNSearcherVar1;
		KdTree originalKdTreeVar2Conditional = kdTreeVar2Conditional;
		UnivariateNearestNeighbourSearcher originalUniNNSearcherVar2 = uniNNSearcherVar2;
		if (variableToReorder == 1) {
			originalData = var1Observations;
			kdTreeVar1Conditional = null; // So that it is rebuilt for the new ordering if required
			uniNNSearcherVar1 = null; // So that it is rebuilt for the new ordering if required
			var1Observations = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData, reordering);
		} else {
			originalData = var2Observations;
			kdTreeVar2Conditional = null; // So that it is rebuilt for the new ordering
			uniNNSearcherVar2 = null; // So that it is rebuilt for the new ordering if required
			var2Observations = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData, reordering);
		}
		// Compute the conditional MI
		double newCondMI = computeFromObservations(false, null)[0];
		// restore original data
		kdTreeJoint = originalKdTreeJoint;
		if (variableToReorder == 1) {
			var1Observations = originalData;
			kdTreeVar1Conditional = originalKdTreeVar1Conditional;
			uniNNSearcherVar1 = originalUniNNSearcherVar1;
		} else {
			var2Observations = originalData;
			kdTreeVar2Conditional = originalKdTreeVar2Conditional;
			uniNNSearcherVar2 = originalUniNNSearcherVar2;
		}
		return newCondMI;
	}

	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] localValues = computeFromObservations(true, null);
		lastAverage = MatrixUtils.mean(localValues);
		condMiComputed = true;
		return localValues;
	}

	/**
	 * @returns the series of local conditional MI values in nats, not bits.
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(double[][] states1,
			double[][] states2, double[][] condStates) throws Exception {
		// Do normalisation of the incoming data if required:
		double[][] states1ToUse, states2ToUse, condStatesToUse;
		if (normalise) {
			states1ToUse = MatrixUtils.normaliseIntoNewArray(states1, var1Means, var1Stds);
			states2ToUse = MatrixUtils.normaliseIntoNewArray(states2, var2Means, var2Stds);
			if (dimensionsCond != 0) {
				condStatesToUse = MatrixUtils.normaliseIntoNewArray(condStates, condMeans, condStds);
			} else {
				condStatesToUse = null;
			}
		} else {
			states1ToUse = states1;
			states2ToUse = states2;
			condStatesToUse = condStates;
		}
		// And call the algorithm:
		double[] localValues = computeFromObservations(true,
				new double[][][]{states1ToUse, states2ToUse, condStatesToUse});
		return localValues;
	}
	
	/**
	 * This protected method handles the multiple threads which
	 *  computes either the average or local conditional MI (over parts of the total
	 *  observations), computing the x, y and z 
	 *  distances between all tuples in time.
	 * 
	 * <p>The method returns:<ol>
	 *  <li>for (returnLocals == false), an array of size 1,
	 *      containing the average conditional MI </li>
	 *  <li>for local conditional MIs (returnLocals == true), the array of local conditional MI values</li>
	 *  </ol>
	 * 
	 * @param returnLocals whether to return an array or local values, or else
	 *  sums of these values
	 * @param newObservations set to null for computing for the observation set for the PDF, or pass in a new set
	 *  of observations to compute the average/locals for (using the existing observations to construct the PDF)
	 * @return either the average conditional MI, or array of local conditional MI value, in nats not bits
	 * @throws Exception
	 */
	protected double[] computeFromObservations(boolean returnLocals, double[][][] newObservations) throws Exception {
		int N = var1Observations.length; // number of observations for the PDFs
		
		double[] returnValues = null;
		
    // How many time points are we averaging over?
    int numTimePointsToComputeFor = (newObservations == null) ?
        N : newObservations[0].length;
    
    if (useGPU && (newObservations == null)) {
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
            newObservations[0], newObservations[1], newObservations[2], returnLocals);
      }
    } else {
      // We're going multithreaded:
      ensureKdTreesConstructed();

      if (returnLocals) {
        // We're computing local conditional MI
        returnValues = new double[numTimePointsToComputeFor];
      } else {
        // We're computing average conditional MI
        returnValues = new double[CondMiKraskovThreadRunner.RETURN_ARRAY_LENGTH];
      }
      
      // Distribute the observations to the threads for the parallel processing
      int lTimesteps = numTimePointsToComputeFor / numThreads; // each thread gets the same amount of data
      int res = numTimePointsToComputeFor % numThreads; // the first thread gets the residual data
      if (debug) {
        System.out.printf("Computing Kraskov conditional MI with %d threads (%d timesteps each, plus %d residual)\n",
            numThreads, lTimesteps, res);
      }
      Thread[] tCalculators = new Thread[numThreads];
      CondMiKraskovThreadRunner[] runners = new CondMiKraskovThreadRunner[numThreads];
      for (int t = 0; t < numThreads; t++) {
        int startTime = (t == 0) ? 0 : lTimesteps * t + res;
        int numTimesteps = (t == 0) ? lTimesteps + res : lTimesteps;
        if (debug) {
          System.out.println(t + ".Thread: from " + startTime +
              " to " + (startTime + numTimesteps)); // Trace Message
        }
        runners[t] = new CondMiKraskovThreadRunner(this, startTime, numTimesteps, newObservations, returnLocals);
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
          MatrixUtils.addInPlace(returnValues, runners[t].getReturnValues());
        }
      }
    }
		
		// Finalise the results:
		if (returnLocals) {
			return returnValues;
		} else {
			// Average out the components for the final equation(s) and for debugging:
			double averageDiGammas = returnValues[CondMiKraskovThreadRunner.INDEX_SUM_DIGAMMAS] / (double) numTimePointsToComputeFor;
			double avNxz = returnValues[CondMiKraskovThreadRunner.INDEX_SUM_NXZ] / (double) numTimePointsToComputeFor;
			double avNyz = returnValues[CondMiKraskovThreadRunner.INDEX_SUM_NYZ] / (double) numTimePointsToComputeFor;
			double avNz = returnValues[CondMiKraskovThreadRunner.INDEX_SUM_NZ] / (double) numTimePointsToComputeFor;
			if (debug) {
				System.out.printf("<n_xz>=%.3f, <n_yz>=%.3f, <n_z>=%.3f\n",
						avNxz, avNyz, avNz);
			}
			if (this.isAlgorithm1) {
				// Algorithm 1:
				if (debug) {
					System.out.printf("Av = digamma(k)=%.3f + <digammas>=%.3f = %.3f \n",
							MathsUtils.digamma(k), averageDiGammas, MathsUtils.digamma(k) + averageDiGammas);
				}
				double[] result = new double[1];
				result[0] = MathsUtils.digamma(k) + averageDiGammas;
				return result;
			} else {
				// Algorithm 2:
				// We also retrieve the sums of inverses for debugging purposes:
				double averageInverseCountInJointXZ =
						returnValues[CondMiKraskovThreadRunner.INDEX_SUM_INV_NXZ] / (double) numTimePointsToComputeFor;
				double averageInverseCountInJointYZ =
						returnValues[CondMiKraskovThreadRunner.INDEX_SUM_INV_NYZ] / (double) numTimePointsToComputeFor;
				double inverseKTerm;
				if (dimensionsCond > 0) {
					// We will add in the 2/K term as usual
					inverseKTerm = 2.0 / (double) k;
				} else {
					// We will only add in 1/K term since without a conditional there is
					//  technically one less variable in the full joint space
					inverseKTerm = 1.0 / (double) k;
				}
				double averageMeasure = MathsUtils.digamma(k) - inverseKTerm + averageDiGammas +
						averageInverseCountInJointXZ + averageInverseCountInJointYZ;
				if (debug) {
					System.out.printf("Av = digamma(k)=%.3f + <digammas>=%.3f +<inverses>=%.3f - $d/k=%.3f  = %.3f" +
							" (<1/n_yz>=%.3f, <1/n_xz>=%.3f)\n",
							MathsUtils.digamma(k), averageDiGammas,
							averageInverseCountInJointXZ + averageInverseCountInJointYZ,
							dimensionsCond > 0 ? 2 : 1,
							inverseKTerm, averageMeasure,
							averageInverseCountInJointYZ, averageInverseCountInJointXZ);
				}
				double[] result = new double[1];
				result[0] = averageMeasure;
				return result;
			}			
		}
	}
	
	/**
	 * Protected method to be used internally for threaded implementations.
	 * This method implements the guts of each Kraskov algorithm, computing the number of 
	 *  nearest neighbours in each dimension for a sub-set of the data points.
	 *  It is intended to be called by one thread to work on that specific
	 *  sub-set of the data.
	 * 
	 * <p>The method returns:<ol>
	 *  <li>for average conditional MIs (returnLocals == false), the relevant sums of digamma(n_{xz}), digamma(n_{yz})
	 *     and digamma(n_z)
	 *     for a partial set of the observations</li>
	 *  <li>for local conditional MIs (returnLocals == true), the array of local conditional MI values</li>
	 *  </ol>
	 * 
	 * @param startTimePoint start time for the partial set we examine
	 * @param numTimePoints number of time points (including startTimePoint to examine)
	 * @param returnLocals whether to return an array or local values, or else
	 *  sums of these values
	 * @return an array of the relevant sum of digamma(n_xz+1) and digamma(n_yz+1) and digamma(n_z), then
	 *  sum of n_xz, n_yz, n_z and for algorithm 2 only, sum of 1/n_xz and 1/n_yz
	 *  (these latter five are only for debugging purposes).
	 * @throws Exception
	 */
	protected abstract double[] partialComputeFromObservations(
			int startTimePoint, int numTimePoints, boolean returnLocals) throws Exception;

  /**
   * Protected method to be used internally for GPU implementations.
   * This method serves the same purpose as partialComputeFromObservations,
   *   but for GPU computation. Each algorithm must override this method
   *   and implement a GPU routine to calculate all values in a single call
   *   to the GPU code.
   */
  protected double[] gpuComputeFromObservations(int startTimePoint,
      int numTimePoints, boolean returnLocals, int nb_surrogates,
      int[][] newOrderings, int variableToReorder) throws Exception {

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
      res = CMIKraskov(totalObservations, var1Observations, dimensionsVar1,
          var2Observations, dimensionsVar2, condObservations, dimensionsCond,
          k, dynCorrExclTime, returnLocals, useMaxNorm, isAlgorithm1,
          nb_surrogates, null!=newOrderings, newOrderings, variableToReorder);
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
    return gpuComputeFromObservations(startTimePoint, numTimePoints, returnLocals, 0, null, -1);
  }

  /**
   * Native method to calculate CMI in GPU.
   */
  private native double[] CMIKraskov(
      int N, double[][] source, int dimx, double[][] dest, int dimy,
      double[][] cond, int dimz,
      int k, int theiler, boolean returnLocals, boolean useMaxNorm,
      boolean isAlgorithm1, int nb_surrogates, boolean orderingsGiven,
      int[][] newOrderings, int variableToReorder);

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
        String errmsg = "GPU library not found. To compile GPU code set the enablegpu flag to true in build.xml";
        if (gpuLibraryPath.length() > 0) {
          errmsg += "\nGPU library was not found in the path provided. Provide full path including library file name.";
          errmsg += "\nExample: /home/johndoe/myfolder/libKraskov.so";
        }
        throw new Exception(errmsg);
      }

      cudaLibraryLoaded = true;

    }
  }

  /**
   * Internal method to ensure that the Kd-tree data structures to represent the
   * observational data have been constructed (should be called prior to attempting
   * to use these data structures)
   */
  protected void ensureKdTreesConstructed() throws Exception {

    // We need to construct the k-d trees for use by the child
    //  classes. We check each tree for existence separately
    //  some can be used across original and surrogate data
    // TODO can parallelise these -- best done within the kdTree --
    //  though it's unclear if there's much point given that
    //  the tree construction itself afterwards can't really be well parallelised.
    if (kdTreeJoint == null) {
      kdTreeJoint = new KdTree(
          new int[] {dimensionsVar1, dimensionsVar2, dimensionsCond},
          new double[][][] {var1Observations, var2Observations, condObservations},
          observationSetIndices, observationTimePoints);
      kdTreeJoint.setNormType(normType);
    }
    if (dimensionsVar1 > 1) {
      if (kdTreeVar1Conditional == null) {
        kdTreeVar1Conditional = new KdTree(
            new int[] {dimensionsVar1, dimensionsCond},
            new double[][][] {var1Observations, condObservations},
            observationSetIndices, observationTimePoints);
        kdTreeVar1Conditional.setNormType(normType);
      } 
    } else { // Univariate variable 1, so we'll search its space alone as this is faster
      if (uniNNSearcherVar1 == null) {
        uniNNSearcherVar1 = new UnivariateNearestNeighbourSearcher(var1Observations,
                observationSetIndices, observationTimePoints);
      }
    }
    if (dimensionsVar2 > 1) {
      if (kdTreeVar2Conditional == null) {
        kdTreeVar2Conditional = new KdTree(
            new int[] {dimensionsVar2, dimensionsCond},
            new double[][][] {var2Observations, condObservations},
            observationSetIndices, observationTimePoints);
        kdTreeVar2Conditional.setNormType(normType);
      }
    } else { // Univariate variable 2, so we'll search its space alone as this is faster
      if (uniNNSearcherVar2 == null) {
        uniNNSearcherVar2 = new UnivariateNearestNeighbourSearcher(var2Observations,
                observationSetIndices, observationTimePoints);
      }
    }
    if ((nnSearcherConditional == null) && (dimensionsCond > 0)) {
      nnSearcherConditional = NearestNeighbourSearcher.create(condObservations,
              observationSetIndices, observationTimePoints);
      nnSearcherConditional.setNormType(normType);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public EmpiricalMeasurementDistribution computeSignificance(
      int variableToReorder, int numPermutationsToCheck) throws Exception {
    if (useGPU) {
      double[] res = gpuComputeFromObservations(0, totalObservations, false, numPermutationsToCheck, null, variableToReorder);
      return new EmpiricalMeasurementDistribution(
          MatrixUtils.select(res, 1, res.length - 1), res[0]);
    } else {
      return super.computeSignificance(variableToReorder, numPermutationsToCheck);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public EmpiricalMeasurementDistribution computeSignificance(
      int variableToReorder, int[][] newOrderings) throws Exception {
    if (useGPU) {
      double[] res = gpuComputeFromObservations(0, totalObservations, false, newOrderings.length, newOrderings, variableToReorder);
      return new EmpiricalMeasurementDistribution(
          MatrixUtils.select(res, 1, res.length - 1), res[0]);
    } else {
      return super.computeSignificance(variableToReorder, newOrderings);
    }
  }
  


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
	 *  <li>for average conditional MIs (returnLocals == false), the relevant sums of digamma(n_{xz}), digamma(n_{yz})
	 *     and digamma(n_z)
	 *     for a partial set of the observations</li>
	 *  <li>for local conditional MIs (returnLocals == true), the array of local conditional MI values</li>
	 *  </ol>
	 * 
	 * @param startTimePoint start time for the partial set we examine
	 * @param numTimePoints number of time points (including startTimePoint to examine)
	 * @param newVar1Observations new time series of observations for variable 1
	 * @param newVar2Observations new time series of observations for variable 1
	 * @param newCondObservations new time series of observations for variable 1
	 * @param returnLocals whether to return an array or local values, or else
	 *  sums of these values
	 * @return an array of the relevant sum of digamma(n_xz+1) and digamma(n_yz+1) and digamma(n_z), then
	 *  sum of n_xz, n_yz, n_z and for algorithm 2 only, sum of 1/n_xz and 1/n_yz
	 *  (these latter five are only for debugging purposes).
	 * @throws Exception
	 */
	protected abstract double[] partialComputeFromNewObservations(
			int startTimePoint, int numTimePoints,
			double[][] newVar1Observations, double[][] newVar2Observations, double[][] newCondObservations, 
			boolean returnLocals) throws Exception;

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
	private class CondMiKraskovThreadRunner implements Runnable {
		protected ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc;
		protected int myStartTimePoint;
		protected int numberOfTimePoints;
		protected double[][][] newObservations;
		protected boolean computeLocals;
		
		protected double[] returnValues = null;
		protected Exception problem = null;
		
		public static final int INDEX_SUM_DIGAMMAS = 0;
		public static final int INDEX_SUM_NXZ = 1;
		public static final int INDEX_SUM_NYZ = 2;
		public static final int INDEX_SUM_NZ = 3;
		public static final int INDEX_SUM_INV_NXZ = 4; // Only used for algorithm 2
		public static final int INDEX_SUM_INV_NYZ = 5; // Only used for algorithm 2
		public static final int RETURN_ARRAY_LENGTH = 6;
		
		public CondMiKraskovThreadRunner(
				ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc,
				int myStartTimePoint, int numberOfTimePoints,
				double[][][] newObservations,
				boolean computeLocals) {
			this.condMiCalc = condMiCalc;
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
		 * @return an exception previously encountered by this thread.
		 * @throws Exception
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
					returnValues = condMiCalc.partialComputeFromObservations(myStartTimePoint, numberOfTimePoints, computeLocals);
				} else {
					// Computing on new observations
					returnValues = condMiCalc.partialComputeFromNewObservations(
							myStartTimePoint, numberOfTimePoints,
							newObservations[0], newObservations[1],
							newObservations[2], computeLocals);
				}
			} catch (Exception e) {
				// Store the exception for later retrieval
				problem = e;
				return;
			}
		}
	}
	// end class MiKraskovThreadRunner

	// Note: no extra implementation of clone provided; we're simply
	//  allowing clone() to produce a shallow copy, which is fine
	//  for the statistical significance calculation (none of the array
	//  data will be changed there.
	//
	// public ConditionalMutualInfoCalculatorMultiVariateKraskov clone() {
	//	return this;
	// }
	
    /**
     * Debug method to return the k nearest neighbour distances that 
     *  would be utilised for each sample point here.
     * Note that this is specifically the max-norm across the three variables, which is used for each
     *  range search in algorithm 1 (although algorithm 2 would use the max distance
     *  for each variable within the kNNs in their separate range searches). 
     *  
     * @param startTimePoint
     * @param numTimePoints
     * @return
     * @throws Exception
     */
	public double[] kNNDistances(int startTimePoint, int numTimePoints) throws Exception {
		double[] kNNdistances = new double[numTimePoints];
		
		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps for this time step by
			//  finding the kth closest neighbour for point t:
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreeJoint.findKNearestNeighbours(k, t, dynCorrExclTime);
			// First element in the PQ is the kth NN,
			//  and epsilon = kthNnData.distance
			NeighbourNodeData kthNnData = nnPQ.poll();

			kNNdistances[t - startTimePoint] = kthNnData.distance;
		}
		return kNNdistances;
	}

	/**
	 * Debug method to return the k nearest neighbour distances that
	 *  would be utilised in {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}
	 *  for a cross conditional MI.
     * Note that this is specifically the max-norm across the three variables, which is used for each
     *  range search in algorithm 1 (although algorithm 2 would use the max distance
     *  for each variable within the kNNs in their separate range searches). 
	 *  
	 * @param startTimePoint
	 * @param numTimePoints
	 * @param newVar1Observations
	 * @param newVar2Observations
	 * @return
	 * @throws Exception
	 */
	public double[] kNNDistancesForNewSamples(int startTimePoint, int numTimePoints,
			double[][] newStates1, double[][] newStates2, double[][] newCondStates) throws Exception {

		double[] kNNdistances = new double[numTimePoints];
		
		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps for this time step by
			//  finding the kth closest neighbour for the new sample:
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreeJoint.findKNearestNeighbours(k,
							new double[][] {newStates1[t], newStates2[t], newCondStates[t]});
			// First element in the PQ is the kth NN,
			//  and epsilon = kthNnData.distance
			NeighbourNodeData kthNnData = nnPQ.poll();
			
			kNNdistances[t - startTimePoint] = kthNnData.distance;
		}
		return kNNdistances;
	}

	/**
	 * As per {@link #kNNDistancesForNewSamples(int, int, double[][], double[][], double[][])}
	 * but for univariates
	 *  
	 * @param startTimePoint
	 * @param numTimePoints
	 * @param newVar1Observations
	 * @param newVar2Observations
	 * @return
	 * @throws Exception
	 */
	public double[] kNNDistancesForNewSamples(int startTimePoint, int numTimePoints,
			double[] newStates1, double[] newStates2, double[][] newCondStates) throws Exception {

		if ((dimensionsVar1 != 1) || (dimensionsVar2 != 1)) {
			throw new Exception("The number of source and dest dimensions (having been initialised to " +
					dimensionsVar1 + " and " + dimensionsVar2 + ") can only be 1 when " +
					"the univariate kNNDistancesForNewSamples(int, int, double[],double[],double[][]) " + 
					"method is called");
		}
		return computeLocalUsingPreviousObservations(
				MatrixUtils.reshape(newStates1, newStates1.length, 1),
				MatrixUtils.reshape(newStates2, newStates2.length, 1),
				newCondStates);
	}

	/**
	 * As per {@link #kNNDistancesForNewSamples(int, int, double[][], double[][], double[][])}
	 * but for univariates
	 *  
	 * @param startTimePoint
	 * @param numTimePoints
	 * @param newVar1Observations
	 * @param newVar2Observations
	 * @return
	 * @throws Exception
	 */
	public double[] kNNDistancesForNewSamples(int startTimePoint, int numTimePoints,
			double[] newStates1, double[] newStates2, double[] newCondStates) throws Exception {

		if ((dimensionsVar1 != 1) || (dimensionsVar2 != 1) || (dimensionsCond != 1)) {
			throw new Exception("The number of source, dest and conditional dimensions (having been initialised to " +
					dimensionsVar1 + ", " + dimensionsVar2 + " and " + dimensionsCond + ") can only be 1 when " +
					"the univariate kNNDistancesForNewSamples(int, int, double[],double[],double[]) " + 
					"method is called");
		}
		return computeLocalUsingPreviousObservations(
				MatrixUtils.reshape(newStates1, newStates1.length, 1),
				MatrixUtils.reshape(newStates2, newStates2.length, 1),
				MatrixUtils.reshape(newCondStates, newCondStates.length, 1));
	}
}
