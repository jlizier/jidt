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

package infodynamics.measures.continuous.kozachenko;

import java.util.Random;
import java.util.Vector;

import infodynamics.measures.continuous.EntropyCalculatorMultiVariate;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential entropy of a given set of observations
 *  (implementing {@link EntropyCalculatorMultiVariate}, using 
 *  the Kozachenko-Leonenko estimator.
 *  For details, see references below.
 *  This class computes it exactly as in the paper by Kraskov et al. below,
 *  i.e. using natural units and twice the minimum distance.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link EntropyCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #EntropyCalculatorMultiVariateKozachenko()}.</li>
 *  <li>An additional {@link #setObservations(double[][], double[][])} option;</li>
 * </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Kozachenko, L., Leonenko, N., "A statistical estimate for the entropy of a random vector",
 *   Problems of Information Transmission, 23 (1987) 9-16.</li>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 *  <li>George Mathews, Hugh Durrant-Whyte, and Mikhail Prokopenko,
 *     <a href="http://dx.doi.org/10.1007/11554028_81">"Measuring Global Behaviour of Multi-Agent Systems from
 *  	Pair-Wise Mutual Information"</a>, 
 *      Knowledge-Based Intelligent Information and Engineering Systems,
 *      Lecture Notes in Computer Science Volume 3684, 2005, pp 587-594.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EntropyCalculatorMultiVariateKozachenko  
	implements EntropyCalculatorMultiVariate {

	/**
	 * Total number of observations supplied.
	 */
	private int totalObservations = 0;
	/**
	 * Number of dimensions of our multivariate data set
	 */
	private int dimensions = 1;
	/**
	 * The set of observations, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 */
	protected double[][] rawData;
	/**
	 * Store the last computed average H
	 */
	private double lastAverage = 0.0;
	/**
	 * Store the last computed local H
	 */
	private double[] lastLocalEntropy;
	/**
	 * Track whether we've computed the average for the supplied
	 *  observations yet
	 */
	private boolean isComputed;
  /**
	 * Storage for observations supplied via {@link #addObservations(double[][])}
	 * type calls
	 */
	protected Vector<double[][]> vectorOfObservations;
	/**
	 * Whether to report debug messages or not
	 */
	protected boolean debug = false;

  /**
   * Property name for whether to normalise the incoming data to 
   * mean 0, standard deviation 1 (default true)
   */
  public static final String PROP_NORMALISE = "NORMALISE";
  /**
   * Property name for an amount of random Gaussian noise to be
   *  added to the data (default is 1e-8, matching the MILCA toolkit).
   */
  public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";

  /**
   * Whether to add an amount of random noise to the incoming data
   */
  protected boolean addNoise = true;
  /**
   * Amount of random Gaussian noise to add to the incoming data
   */
  protected double noiseLevel = (double) 1e-8;
  /**
   * Stored pre-computed value of the Euler-Mascheroni constant
   */
	public static final double EULER_MASCHERONI_CONSTANT = 0.5772156;
	
	/**
	 * Construct an instance
	 */
	public EntropyCalculatorMultiVariateKozachenko() {
		totalObservations = 0;
		isComputed = false;
		lastLocalEntropy = null;
	}

	@Override
	public void initialise() throws Exception {
		initialise(dimensions);
	}

	public void initialise(int dimensions) {
		this.dimensions = dimensions;
		rawData = null;
		totalObservations = 0;
		isComputed = false;
		lastLocalEntropy = null;
	}

	@Override
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(NUM_DIMENSIONS_PROP_NAME)) {
			dimensions = Integer.parseInt(propertyValue);
    } else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
      if (propertyValue.equals("0") ||
          propertyValue.equalsIgnoreCase("false")) {
        addNoise = false;
        noiseLevel = 0;
      } else {
        addNoise = true;
        noiseLevel = Double.parseDouble(propertyValue);
      }
		} else {
			// No property was set, and no superclass to call.
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(NUM_DIMENSIONS_PROP_NAME)) {
			return Integer.toString(dimensions);
    } else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
      return Double.toString(noiseLevel);
		} else {
			// No property was set, and no superclass to call.
			return null;
		}
	}

	public void startAddObservations() {
		isComputed = false;
    totalObservations = 0;
		lastLocalEntropy = null;
    rawData = null;
		vectorOfObservations = new Vector<double[][]>();
	}

	public void finaliseAddObservations() {

		rawData = new double[totalObservations][dimensions];
		
		// Construct the joint vectors from the given observations
		//  (removing redundant data which is outside any timeDiff)
		int startObservation = 0;
		for (double[][] obs : vectorOfObservations) {
			// Copy the data from these given observations into our master array
			MatrixUtils.arrayCopy(obs, 0, 0,
					rawData, startObservation, 0,
					obs.length, dimensions);
			startObservation += obs.length;
		}

		// We don't need to keep the vector of observation sets anymore:
		vectorOfObservations = null;

    if (addNoise) {
      Random random = new Random();
      // Add Gaussian noise of std dev noiseLevel to the data
      for (int r = 0; r < totalObservations; r++) {
        for (int c = 0; c < dimensions; c++) {
          rawData[r][c] += random.nextGaussian()*noiseLevel;
        }
      }
    }
	}

	@Override
  public void setObservations(double[][] observations) {
    startAddObservations();
    addObservations(observations);
    finaliseAddObservations();
  }

  public void setObservations(double[][] observations1, double[][] observations2)
      throws Exception {
    startAddObservations();
    addObservations(observations1, observations2);
    finaliseAddObservations();
  }

	/*
	 * (non-Javadoc)
	 * @see infodynamics.measures.continuous.EntropyCalculator#setObservations(double[])
	 * 
	 * This method here to ensure we make compatibility with the 
	 *  EntropyCalculator interface.
	 */
	@Override
  public void setObservations(double[] observations) {
    startAddObservations();
    addObservations(observations);
    finaliseAddObservations();
  }

  public void addObservations(double[][] observations) {
		if (vectorOfObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
    vectorOfObservations.add(observations);
    totalObservations += observations.length;
  }

	public void addObservations(double[] observations) {
		rawData = MatrixUtils.reshape(observations, observations.length, 1);
		addObservations(rawData);
	}
	
	/**
	 * Each row of the data is an observation; each column of
	 *  the row is a new variable in the multivariate observation.
	 * This method signature allows the user to call setObservations for
	 *  joint time series without combining them into a single joint time
	 *  series (we do the combining for them).
	 * 
	 * @param data1 first few variables in the joint data
	 * @param data2 the other variables in the joint data
	 * @throws Exception When the length of the two arrays of observations do not match.
	 * @see #addObservations(double[][])
	 */
	public void addObservations(double[][] data1,
			double[][] data2) throws Exception {
		int timeSteps = data1.length;
		if ((data1 == null) || (data2 == null)) {
			throw new Exception("Cannot have null data arguments");
		}
		if (data1.length != data2.length) {
			throw new Exception("Length of data1 (" + data1.length + ") is not equal to the length of data2 (" +
					data2.length + ")");
		}
		int data1Variables = data1[0].length;
		int data2Variables = data2[0].length;
		double[][] data = new double[timeSteps][data1Variables + data2Variables];
		for (int t = 0; t < timeSteps; t++) {
			System.arraycopy(data1[t], 0, data[t], 0, data1Variables);
			System.arraycopy(data2[t], 0, data[t], data1Variables, data2Variables);
		}
		// Now defer to the normal setObservations method
		addObservations(data);
	}

	/**
	 * @return entropy in natural units
	 */
	@Override
	public double computeAverageLocalOfObservations() {
		if (isComputed) {
			return lastAverage;
		}
		double sdTermHere = sdTerm(totalObservations, dimensions);
		double emConstHere = eulerMascheroniTerm(totalObservations);
		double[] minDistance = EuclideanUtils.computeMinEuclideanDistances(rawData);
		double entropy = 0.0;
		if (debug) {
			System.out.println("t,\tminDist,\tlogMinDist,\tsum");
		}
		for (int t = 0; t < rawData.length; t++) {
			entropy += Math.log(2.0 * minDistance[t]);
			if (debug) {
				System.out.println(t + ",\t" + 
						minDistance[t] + ",\t" +
						Math.log(minDistance[t]) + ",\t" +
						entropy);
			}
		}
		// Using natural units
		// entropy /= Math.log(2);
		entropy *= (double) dimensions / (double) totalObservations;
		if (debug) {
			System.out.println("Sum part:   " + entropy);
			System.out.println("Euler part: " + emConstHere);
			System.out.println("Sd term:    " + sdTermHere);
		}
		entropy += emConstHere;
		entropy += sdTermHere;
		lastAverage = entropy;
		isComputed = true;
		return entropy;
	}

	/**
	 * @return local entropies in natural units
	 */
	@Override
	public double[] computeLocalOfPreviousObservations() {
		if (lastLocalEntropy != null) {
			return lastLocalEntropy;
		}

		double sdTermHere = sdTerm(totalObservations, dimensions);
		double emConstHere = eulerMascheroniTerm(totalObservations);
		double constantToAddIn = sdTermHere + emConstHere;
		
		double[] minDistance = EuclideanUtils.computeMinEuclideanDistances(rawData);
		double entropy = 0.0;
		double[] localEntropy = new double[rawData.length];
		if (debug) {
			System.out.println("t,\tminDist,\tlogMinDist,\tlocal,\tsum");
		}
		for (int t = 0; t < rawData.length; t++) {
			localEntropy[t] = Math.log(2.0 * minDistance[t]) * (double) dimensions;
			// using natural units
			// localEntropy[t] /= Math.log(2);
			localEntropy[t] += constantToAddIn;
			entropy += localEntropy[t];
			if (debug) {
				System.out.println(t + ",\t" + 
						minDistance[t] + ",\t" +
						Math.log(minDistance[t]) + ",\t" +
						localEntropy[t] + ",\t" +
						entropy);
			}
		}
		entropy /= (double) totalObservations;
		lastAverage = entropy;
		lastLocalEntropy = localEntropy;
		return localEntropy;
	}

	/**
	 * Not implemented yet
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(double[][] states) throws Exception {
		throw new Exception("Local method for other data not implemented");
	}

	/**
	 * Not implemented yet
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] states1, double[][] states2) throws Exception {
		throw new Exception("Local method for other data not implemented");
	}

	/**
	 * Returns the value of the Euler-Mascheroni term.
	 * Public for debugging purposes
	 * 
	 * @return
	 */
	public double eulerMascheroniTerm(int N) {
		// Using natural units
		// return EULER_MASCHERONI_CONSTANT / Math.log(2);
		try {
			return -MathsUtils.digamma(1) + MathsUtils.digamma(N);
		} catch (Exception e) {
			// Exception will only be thrown if N < 0
			return 0;
		}
	}
	
	/**
	 * Returns the value of the Sd term
	 * Public for debugging purposes
	 * 
	 * @param numObservations
	 * @param dimensions
	 * @return
	 */
	public double sdTerm(int numObservations, int dimensions) {
		// To compute directly:
		// double unLoggedSdTerm = 
		//	Math.pow(Math.PI/4.0, ((double) dimensions) / 2.0) /	// Brought 2^d term from denominator into Pi term
		//	MathsUtils.gammaOfArgOn2Plus1(dimensions);
		
		// But we need to compute it carefully, to allow the maximum range of dimensions
		// double unLoggedSdTerm = 
		//	1.0 / MathsUtils.gammaOfArgOn2Plus1IncludeDivisor(dimensions, 
		//			Math.pow(Math.PI, ((double) dimensions) / 2.0));
		// Don't include the 2^d in the above divisor, since that makes the divisor < 1, which 
		//  doesn't help at all.
		// unLoggedSdTerm /= Math.pow(2, dimensions);
		// return Math.log(unLoggedSdTerm) / Math.log(2);
		// Using natural units
		// return Math.log(unLoggedSdTerm);
		
		// But even that method falls over by about d = 340.
		// Break down the log into the log of a factorial term and the log of a constant term
		double constantTerm = Math.pow(Math.PI / 4.0, (double) dimensions / 2.0);
		double result = 0.0;
		if (dimensions % 2 == 0) {
			// d even
			// Now take log (1/(d/2)!) = -log (1/(d/2)!) = -sum(d/2 --) {log d/2}
			for (int d = dimensions/2; d > 1; d--) {
				result -= Math.log(d);
			}
		} else {
			// d odd
			constantTerm *= Math.pow(2.0, (double) (dimensions + 1) / 2.0);
			constantTerm /= Math.sqrt(Math.PI);
			// Now take log (1/d!!) = - log (d!!) = - sum(d -= 2) {log d}
			for (int d = dimensions; d > 1; d -= 2) {
				result -= Math.log(d);
			}
		}
		result += Math.log(constantTerm);
		return result;
	}
	
	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	@Override
	public double getLastAverage() {
		return lastAverage;
	}

	@Override
	public int getNumObservations() {
		return totalObservations;
	}
}
