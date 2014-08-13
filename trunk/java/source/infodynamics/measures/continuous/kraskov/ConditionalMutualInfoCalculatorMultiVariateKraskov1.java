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

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential conditional mutual information of two multivariate
 *  <code>double[][]</code> sets of observations, conditioned on another
 *  (implementing {@link ConditionalMutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see references below)
 *  <b>algorithm 1</b>.
 *  Most of the functionality is defined by the parent class 
 *  {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}.</p>
 *
 * <p>Crucially, the calculation is performed by examining
 * neighbours in the full joint space (as specified by Frenzel and Pompe)
 * rather than two MI calculators.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link ConditionalMutualInfoCalculatorMultiVariate},
 * and expanded on in {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}.
 * </p>
 *  
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Frenzel and Pompe, <a href="http://dx.doi.org/10.1103/physrevlett.99.204101">
 * 	"Partial Mutual Information for Coupling Analysis of Multivariate Time Series"</a>,
 * 	Physical Review Letters, <b>99</b>, p. 204101+ (2007).</li>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class ConditionalMutualInfoCalculatorMultiVariateKraskov1
	extends ConditionalMutualInfoCalculatorMultiVariateKraskov {
	
	/**
	 * Multiplier used as hueristic for determining whether to use a linear search
	 *  for kth nearest neighbour or a binary search.
	 */
	protected static final double CUTOFF_MULTIPLIER = 1.5;

	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		return computeAverageLocalOfObservations(1, null);
	}

	/**
	 * See description in {@link ConditionalMutualInfoCalculatorMultiVariate#computeAverageLocalOfObservations(int, int[])} 
	 *   
	 * If {@code reordering} is null, it is assumed there is no reordering of
	 *  the given variable.
	 */
	@Override
	public double computeAverageLocalOfObservations(int variableToReorder, 
			int[] reordering) throws Exception {
		if (!tryKeepAllPairsNorms || (var1Observations.length > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			double[][] originalData;
			if (variableToReorder == 1) {
				originalData = var1Observations;
			} else {
				originalData = var2Observations;
			}
			if (reordering != null) {
				// Generate a new re-ordered data array
				if (variableToReorder == 1) {
					var1Observations = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData, reordering);
				} else {
					var2Observations = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData, reordering);
				}
			}
			// Compute the MI
			double newMI = computeAverageLocalOfObservationsWhileComputingDistances();
			// restore original data
			if (variableToReorder == 1) {
				var1Observations = originalData;
			} else {
				var2Observations = originalData;
			}
			return newMI;
		}
		
		// Else we'll use the arrays of marginal distances
		
		if (xNorms == null) {
			computeNorms();
		}
		int N = var1Observations.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		
		// Count the average number of points within eps_xz and eps_yz and eps_z
		double averageDiGammas = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  using x, y and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			
			int tForReorderedVar = (reordering == null) ? t : reordering[t];

			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				int t2ForReorderedVar = (reordering == null) ? t2 : reordering[t2];
				// Joint norm is the max of all three marginals
				if (variableToReorder == 1) {
					jointNorm[t2] = Math.max(xNorms[tForReorderedVar][t2ForReorderedVar],
							Math.max(yNorms[t][t2],
									zNorms[t][t2]));
				} else {
					jointNorm[t2] = Math.max(xNorms[t][t2],
							Math.max(yNorms[tForReorderedVar][t2ForReorderedVar],
									zNorms[t][t2]));
				}
			}
			// Then find the kth closest neighbour, using a heuristic to 
			// select whether to keep the k mins only or to do a sort.
			double epsilon = 0.0;
			if (k <= cutoffForKthMinLinear) {
				// just do a linear search for the minimum
				epsilon = MatrixUtils.kthMin(jointNorm, k);
			} else {
				// Sort the array of joint norms first
				java.util.Arrays.sort(jointNorm);
				// And find the distance to it's kth closest neighbour
				// (we subtract one since the array is indexed from zero)
				epsilon = jointNorm[k-1];
			}
			
			// Count the number of points whose (x,z) distance is less
			//  than eps, and whose (y,z) distance is less than eps, and
			//  whose z distance is less than eps.
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (zNorms[t][t2] < epsilon) {
					n_z++;
					int t2ForReorderedVar = (reordering == null) ? t2 : reordering[t2];
					if (variableToReorder == 1) {
						if (xNorms[tForReorderedVar][t2ForReorderedVar] < epsilon) {
							n_xz++;
						}
						if (yNorms[t][t2] < epsilon) {
							n_yz++;
						}
					} else {
						if (xNorms[t][t2] < epsilon) {
							n_xz++;
						}
						if (yNorms[tForReorderedVar][t2ForReorderedVar] < epsilon) {
							n_yz++;
						}
					}
				}
			}
			// And take the digamma before adding into the 
			//  average:
			// Note: we're using digamma function which has opposite sign to the harmonic
			//  number used by Frenzel and Pompe, and is also offset by a constant (though
			//  this cancels out)
			averageDiGammas += MathsUtils.digamma(n_z+1) - MathsUtils.digamma(n_xz+1)
							- MathsUtils.digamma(n_yz+1);
			if (debug) {
				avNxz += n_xz;
				avNyz += n_yz;
				avNz += n_z;
				double localCondMi = MathsUtils.digamma(k) +
						MathsUtils.digamma(n_z+1) - MathsUtils.digamma(n_xz+1)
						- MathsUtils.digamma(n_yz+1);
				System.out.printf("t=%d, n_xz=%d, n_yz=%d, n_z=%d, local=%.4f\n",
						t, n_xz, n_yz, n_z, localCondMi);
			}
		}
		averageDiGammas /= (double) N;
		lastAverage = MathsUtils.digamma(k) + averageDiGammas;
		condMiComputed = true;
		
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.printf("<n_xz>=%.3f, <n_yz>=%.3f, <n_z>=%.3f\n",
					avNxz, avNyz, avNz);
			System.out.printf("Av = digamma(k)=%.3f + <digammas>=%.3f = %.3f \n",
					MathsUtils.digamma(k), averageDiGammas, lastAverage);
		}
		
		return lastAverage;
	}
	
	/**
	 * This method correctly computes the average conditional MI, but recomputes the
	 *  distances between all tuples in time.
	 * Kept here for cases where we have too many observations
	 *  to keep the norm between all pairs, and for testing purposes.
	 * 
	 * @see #computeAverageLocalOfObservations()
	 * @return average conditional MI value in nats not bits
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservationsWhileComputingDistances() throws Exception {
		int N = var1Observations.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		
		// Count the average number of points within eps_xz, eps_yz and eps_z
		double averageDiGammas = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  First get x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xyzNorms = normCalculator.computeNorms(var1Observations, var2Observations, condObservations, t);
			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2] = Math.max(xyzNorms[t2][0], Math.max(xyzNorms[t2][1], xyzNorms[t2][2]));
			}
			// Then find the kth closest neighbour, using a heuristic to 
			// select whether to keep the k mins only or to do a sort.
			double epsilon = 0.0;
			if (k <= cutoffForKthMinLinear) {
				// just do a linear search for the minimum
				epsilon = MatrixUtils.kthMin(jointNorm, k);
			} else {
				// Sort the array of joint norms first
				java.util.Arrays.sort(jointNorm);
				// And find the distance to it's kth closest neighbour
				// (we subtract one since the array is indexed from zero)
				epsilon = jointNorm[k-1];
			}
			
			// Count the number of points whose (x,z) distance is less
			//  than eps, whose (y,z) distance is less than eps, and whose
			//  z distance is less than eps
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyzNorms[t2][2] < epsilon) {
					n_z++;
					if (xyzNorms[t2][0] < epsilon) {
						n_xz++;
					}
					if (xyzNorms[t2][1] < epsilon) {
						n_yz++;
					}
				}
			}
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_z+1) - MathsUtils.digamma(n_xz+1)
						- MathsUtils.digamma(n_yz+1);
			if (debug) {
				avNxz += n_xz;
				avNyz += n_yz;
				avNz += n_z;
				double localCondMi = MathsUtils.digamma(k) +
						MathsUtils.digamma(n_z+1) - MathsUtils.digamma(n_xz+1)
						- MathsUtils.digamma(n_yz+1);
				System.out.printf("t=%d, n_xz=%d, n_yz=%d, n_z=%d, local=%.4f\n",
						t, n_xz, n_yz, n_z, localCondMi);
			}
		}
		averageDiGammas /= (double) N;
		lastAverage = MathsUtils.digamma(k) + averageDiGammas;
		condMiComputed = true;

		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.printf("<n_xz>=%.3f, <n_yz>=%.3f, <n_z>=%.3f\n",
					avNxz, avNyz, avNz);
			System.out.printf("Av = digamma(k)=%.3f + <digammas>=%.3f = %.3f \n",
					MathsUtils.digamma(k), averageDiGammas, lastAverage);
		}
		
		return lastAverage;
	}

	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		int N = var1Observations.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		double[] localCondMi = new double[N];
		
		// Constants:
		double digammaK = MathsUtils.digamma(k);
		
		// Count the average number of points within eps_xz and eps_yz and eps_z
		double averageDiGammas = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps for this time step:
			//  First get x and y and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity.
			double[][] xyzNorms = normCalculator.computeNorms(var1Observations, var2Observations, condObservations, t);
			double[] jointNorm = new double[N];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2] = Math.max(xyzNorms[t2][0], Math.max(xyzNorms[t2][1], xyzNorms[t2][2]));
			}
			// Then find the kth closest neighbour, using a heuristic to 
			// select whether to keep the k mins only or to do a sort.
			double epsilon = 0.0;
			if (k <= cutoffForKthMinLinear) {
				// just do a linear search for the minimum
				epsilon = MatrixUtils.kthMin(jointNorm, k);
			} else {
				// Sort the array of joint norms first
				java.util.Arrays.sort(jointNorm);
				// And find the distance to it's kth closest neighbour
				// (we subtract one since the array is indexed from zero)
				epsilon = jointNorm[k-1];
			}
			
			// Count the number of points whose x distance is less
			//  than eps, and whose y distance is less than eps
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyzNorms[t2][2] < epsilon) {
					n_z++;
					if (xyzNorms[t2][0] < epsilon) {
						n_xz++;
					}
					if (xyzNorms[t2][1] < epsilon) {
						n_yz++;
					}
				}
			}
			// And take the digamma:
			double digammaNxzPlusOne = MathsUtils.digamma(n_xz+1);
			double digammaNyzPlusOne = MathsUtils.digamma(n_yz+1);
			double digammaNzPlusOne = MathsUtils.digamma(n_z+1);
			
			localCondMi[t] = digammaK - digammaNxzPlusOne - digammaNyzPlusOne + digammaNzPlusOne;
			
			// And keep track of the average
			averageDiGammas += digammaNzPlusOne - digammaNxzPlusOne - digammaNyzPlusOne;
			if (debug) {
				avNxz += n_xz;
				avNyz += n_yz;
				avNz += n_z;
				System.out.printf("t=%d, n_xz=%d, n_yz=%d, n_z=%d, local=%.4f\n",
						t, n_xz, n_yz, n_z, localCondMi[t]);
			}
		}
		averageDiGammas /= (double) N;
		lastAverage = digammaK + averageDiGammas;
		condMiComputed = true;
		
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.printf("<n_xz>=%.3f, <n_yz>=%.3f, <n_z>=%.3f\n",
					avNxz, avNyz, avNz);
			System.out.printf("Av = digamma(k)=%.3f + <digammas>=%.3f = %.3f \n",
					digammaK, averageDiGammas, lastAverage);
		}
		
		return localCondMi;
	}

	@Override
	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e",
				k, MathsUtils.digamma(k));
		return constants;
	}
}
