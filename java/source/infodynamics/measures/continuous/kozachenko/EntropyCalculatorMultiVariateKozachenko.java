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

import infodynamics.measures.continuous.EntropyCalculatorMultiVariate;
import infodynamics.measures.continuous.EntropyCalculatorMultiVariateCommon;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MathsUtils;

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
	extends EntropyCalculatorMultiVariateCommon 
	implements EntropyCalculatorMultiVariate {


  /**
   * Stored pre-computed value of the Euler-Mascheroni constant
   */
	public static final double EULER_MASCHERONI_CONSTANT = 0.5772156;
	
	/**
	 * Construct an instance
	 */
	public EntropyCalculatorMultiVariateKozachenko() {
		super();
		noiseLevel = (double) 1e-8; // Default to align with KSG estimators
	    addNoise = true;
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
		double[] minDistance = EuclideanUtils.computeMinEuclideanDistances(observations);
		double entropy = 0.0;
		if (debug) {
			System.out.println("t,\tminDist,\tlogMinDist,\tsum");
		}
		for (int t = 0; t < observations.length; t++) {
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

		double sdTermHere = sdTerm(totalObservations, dimensions);
		double emConstHere = eulerMascheroniTerm(totalObservations);
		double constantToAddIn = sdTermHere + emConstHere;
		
		double[] minDistance = EuclideanUtils.computeMinEuclideanDistances(observations);
		double entropy = 0.0;
		double[] localEntropy = new double[observations.length];
		if (debug) {
			System.out.println("t,\tminDist,\tlogMinDist,\tlocal,\tsum");
		}
		for (int t = 0; t < observations.length; t++) {
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
		isComputed = true;
		return localEntropy;
	}

	/**
	 * Not implemented yet
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(double[][] states) throws Exception {
		double sdTermHere = sdTerm(totalObservations, dimensions);
		double emConstHere = eulerMascheroniTerm(totalObservations);
		double constantToAddIn = sdTermHere + emConstHere;
		
		double[] minDistance = EuclideanUtils.computeMinEuclideanDistances(observations);
		double entropy = 0.0;
		double[] localEntropy = new double[observations.length];
		if (debug) {
			System.out.println("t,\tminDist,\tlogMinDist,\tlocal,\tsum");
		}
		for (int t = 0; t < observations.length; t++) {
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
		return localEntropy;
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
	
}
