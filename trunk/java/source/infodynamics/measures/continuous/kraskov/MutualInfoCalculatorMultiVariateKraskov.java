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
import java.util.Random;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.MutualInfoMultiVariateCommon;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential mutual information of two given multivariate sets of
 *  observations (implementing {@link MutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below).
 *  This is an abstract class to gather common functionality between the two
 *  algorithms defined by Kraskov et al.
 *  Two child classes {@link MutualInfoCalculatorMultiVariateKraskov1} and
 *  {@link MutualInfoCalculatorMultiVariateKraskov2} then
 *  actually implement the two algorithms in the Kraskov et al. paper</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link MutualInfoCalculatorMultiVariate},
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
	 * Calculator for the norm between data points
	 */
	protected EuclideanUtils normCalculator;
	
	/**
	 * Property name for the number of K nearest neighbours used in
	 * the KSG algorithm in the full joint space (default 4).
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
	 * Property name for whether to normalise the incoming data to 
	 * mean 0, standard deviation 1 (default true)
	 */
	public static final String PROP_NORMALISE = "NORMALISE";
	/**
	 * Property name for an amount of random Gaussian noise to be
	 *  added to the data (default is 0).
	 */
	public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";
	/**
	 * Property name for the number of parallel threads to use in the
	 *  computation
	 */
	public static final String PROP_NUM_THREADS = "NUM_THREADS";
	/**
	 * Valid property value for {@link #PROP_NUM_THREADS} to indicate
	 *  that all available processors should be used. 
	 */
	public static final String USE_ALL_THREADS = "USE_ALL";
	
	/**
	 * Whether to normalise the incoming data 
	 */
	protected boolean normalise = true;
	/**
	 * Whether to add an amount of random noise to the incoming data
	 */
	protected boolean addNoise = false;
	/**
	 * Amount of random Gaussian noise to add to the incoming data
	 */
	protected double noiseLevel = 0.0;
	/**
	 * Number of parallel threads to use in the computation;
	 *  defaults to use all available.
	 */
	protected int numThreads = Runtime.getRuntime().availableProcessors();
	/**
	 * Private variable to record which algorithm this instance is implementing
	 */
	protected boolean isAlgorithm1 = false;
	
	/**
	 * Construct an instance of the KSG MI calculator
	 */
	public MutualInfoCalculatorMultiVariateKraskov() {
		super();
		normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
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
	 * 	<li>{@link #PROP_NORM_TYPE}</li> -- normalization type to apply to 
	 * 		working out the norms between the points in each marginal space.
	 * 		Options are defined by {@link EuclideanUtils#setNormToUse(String)} -
	 * 		default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 *  <li>{@link #PROP_NORMALISE} -- whether to normalise the incoming individual
	 *      variables to mean 0 and standard deviation 1 (true by default)</li>
	 *  <li>{@link #PROP_ADD_NOISE} -- a standard deviation for an amount of
	 *  	random Gaussian noise to add to
	 *      each variable, to avoid having neighbourhoods with artificially
	 *      large counts. The amount is added in after any normalisation,
	 *      so can be considered as a number of standard deviations of the data.
	 *      (Recommended by Kraskov. MILCA uses 1e-8; but adds in
	 *      a random amount of noise in [0,noiseLevel) ). Default 0.</li>
	 *  <li>{@link #PROP_NUM_THREADS} -- the integer number of parallel threads
	 *  	to use in the computation. Can be passed as a string "USE_ALL"
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
			normCalculator.setNormToUse(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
			addNoise = true;
			noiseLevel = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NUM_THREADS)) {
			if (propertyValue.equalsIgnoreCase(USE_ALL_THREADS)) {
				numThreads = Runtime.getRuntime().availableProcessors();
			} else { // otherwise the user has passed in an integer:
				numThreads = Integer.parseInt(propertyValue);
			}
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

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.MutualInfoMultiVariateCommon#finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		// Allow the parent to generate the data for us first
		super.finaliseAddObservations();
		
		// Normalise the data if required
		if (normalise) {
			// We can overwrite these since they're already
			//  a copy of the users' data.
			MatrixUtils.normalise(sourceObservations);
			MatrixUtils.normalise(destObservations);
		}
		
		if (addNoise) {
			Random random = new Random();
			// Add Gaussian noise of std dev noiseLevel to the data
			for (int r = 0; r < sourceObservations.length; r++) {
				for (int c = 0; c < dimensionsSource; c++) {
					sourceObservations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
				for (int c = 0; c < dimensionsDest; c++) {
					destObservations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
			}
		}
	}

	/**
	 * {@inheritDoc} 
	 * 
	 * @return the average MI in nats (not bits!)
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		// Compute the MI
		double startTime = Calendar.getInstance().getTimeInMillis();
		lastAverage = computeFromObservations(false)[0];
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
	 *  Returns NaN if any of the determinants are zero
	 *  (because this will make the denominator of the log 0).
	 */
	public double computeAverageLocalOfObservations(int[] reordering) throws Exception {
		double[][] originalData2 = destObservations;
		if (reordering != null) {
			// Generate a new re-ordered data2
			destObservations = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData2, reordering);
		}
		// Compute the MI
		double newMI = computeFromObservations(false)[0];
		// restore data2
		destObservations = originalData2;
		if (reordering == null) {
			// Only keep this average value if it was for
			//  the original data:
			miComputed = true;
			lastAverage = newMI;
		}
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
	 * @return the "time-series" of local MIs in bits
	 * @throws Exception
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] localValues = computeFromObservations(true);
		lastAverage = MatrixUtils.mean(localValues);
		miComputed = true;
		return localValues;
	}

	/**
	 * This method, specified in {@link MutualInfoCalculatorMultiVariate}
	 * is not implemented yet here.
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] states1, double[][] states2) throws Exception {
		// TODO If implemented, will need to incorporate any time difference here.
		// Will also need to handle normalisation of the incoming data
		//  appropriately
		throw new Exception("Local method not implemented yet");
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
	 * @return either the average MI, or array of local MI value, in nats not bits
	 * @throws Exception
	 */
	protected double[] computeFromObservations(boolean returnLocals) throws Exception {
		int N = sourceObservations.length; // number of observations
		
		double[] returnValues = null;
		
		if (numThreads == 1) {
			// Single-threaded implementation:
			returnValues = partialComputeFromObservations(0, N, returnLocals);
			
		} else {
			// We're going multithreaded:
			if (returnLocals) {
				// We're computing local MI
				returnValues = new double[N];
			} else {
				// We're computing average MI
				returnValues = new double[3];
			}
			
			// Distribute the observations to the threads for the parallel processing
			int lTimesteps = N / numThreads; // each thread gets the same amount of data
			int res = N % numThreads; // the first thread gets the residual data
			if (debug) {
				System.out.printf("Computing Kraskov MI alg1 with %d threads (%d timesteps each, plus %d residual)\n",
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
				runners[t] = new MiKraskovThreadRunner(this, startTime, numTimesteps, returnLocals);
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
			// Compute the average number of points within eps_x and eps_y
			double averageDiGammas = returnValues[MiKraskovThreadRunner.INDEX_SUM_DIGAMMAS] / (double) N;
			double avNx = returnValues[MiKraskovThreadRunner.INDEX_SUM_NX] / (double) N;
			double avNy = returnValues[MiKraskovThreadRunner.INDEX_SUM_NY] / (double) N;
			if (debug) {
				System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
			}
			
			// Finalise the average result, depending on which algorithm we are implementing:
			if (isAlgorithm1) {
				return new double[] { MathsUtils.digamma(k) - averageDiGammas + MathsUtils.digamma(N)};
			} else {
				return new double[] { MathsUtils.digamma(k) - (1.0 / (double)k) - averageDiGammas + MathsUtils.digamma(N)};
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
		protected boolean computeLocals;
		
		protected double[] returnValues = null;
		protected Exception problem = null;
		
		public static final int INDEX_SUM_DIGAMMAS = 0;
		public static final int INDEX_SUM_NX = 1;
		public static final int INDEX_SUM_NY = 2;
		
		public MiKraskovThreadRunner(
				MutualInfoCalculatorMultiVariateKraskov miCalc,
				int myStartTimePoint, int numberOfTimePoints,
				boolean computeLocals) {
			this.miCalc = miCalc;
			this.myStartTimePoint = myStartTimePoint;
			this.numberOfTimePoints = numberOfTimePoints;
			this.computeLocals = computeLocals;
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
				returnValues = miCalc.partialComputeFromObservations(myStartTimePoint, numberOfTimePoints, computeLocals);
			} catch (Exception e) {
				// Store the exception for later retrieval
				problem = e;
				return;
			}
		}
	}
	// end class MiKraskovThreadRunner
	
	/**
	 * Utility function used for debugging, printing digamma constants
	 * 
	 * @param N
	 * @return
	 * @throws Exception
	 */
	public abstract String printConstants(int N) throws Exception ;	
}
