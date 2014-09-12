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

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential mutual information of two given multivariate sets of
 *  observations (implementing {@link MutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below),
 *  <b>algorithm 2</b>.
 *  Most of the functionality is defined by the parent class 
 *  {@link MutualInfoCalculatorMultiVariateKraskov}.</p>
 *
 * <p>Usage is as per the paradigm outlined for {@link MutualInfoCalculatorMultiVariate},
 * and expanded on in {@link MutualInfoCalculatorMultiVariateKraskov}.
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
public class MutualInfoCalculatorMultiVariateKraskov2
	extends MutualInfoCalculatorMultiVariateKraskov {

	protected static final int JOINT_NORM_VAL_COLUMN = 0;
	protected static final int JOINT_NORM_TIMESTEP_COLUMN = 1;

	/**
	 * Multiplier used as hueristic for determining whether to use a linear search
	 *  for kth nearest neighbour or a binary search.
	 */
	protected static final double CUTOFF_MULTIPLIER = 1.5;
	
	public MutualInfoCalculatorMultiVariateKraskov2() {
		super();
		isAlgorithm1 = false;
	}

	protected double[] partialComputeFromObservations(
			int startTimePoint, int numTimePoints, boolean returnLocals) throws Exception {
		
		double startTime = Calendar.getInstance().getTimeInMillis();

		int N = sourceObservations.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));
		
		double[] localMi = null;
		if (returnLocals) {
			localMi = new double[numTimePoints];
		}
		
		// Constants:
		double digammaK = MathsUtils.digamma(k);
		double invK = 1.0 / (double)k;
		double digammaN = MathsUtils.digamma(N);
		
		// Count the average number of points within eps_x and eps_y of each point
		double sumDiGammas = 0;
		double sumNx = 0;
		double sumNy = 0;
				
		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps_x and eps_y for this time step:
			//  First get x and y norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xyNorms = normCalculator.computeNorms(sourceObservations, destObservations, t);
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = Math.max(xyNorms[t2][0], xyNorms[t2][1]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours, using a heuristic to 
			// select whether to keep the k mins only or to do a sort.
			double eps_x = 0.0;
			double eps_y = 0.0;
			int[] timeStepsOfKthMins = null;
			if (k <= cutoffForKthMinLinear) {
				// just do a linear search for the minimum epsilon value
				timeStepsOfKthMins = MatrixUtils.kMinIndices(jointNorm, JOINT_NORM_VAL_COLUMN, k);
			} else {
				// Sort the array of joint norms
				java.util.Arrays.sort(jointNorm, FirstIndexComparatorDouble.getInstance());
				// and now we have the closest k points.
				timeStepsOfKthMins = new int[k];
				for (int j = 0; j < k; j++) {
					timeStepsOfKthMins[j] = (int) jointNorm[j][JOINT_NORM_TIMESTEP_COLUMN];
				}
			}
			// and now we have the closest k points.
			// Find eps_{x,y} as the maximum x and y norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xyNorms[timeStepOfJthPoint][0] > eps_x) {
					eps_x = xyNorms[timeStepOfJthPoint][0];
				}
				if (xyNorms[timeStepOfJthPoint][1] > eps_y) {
					eps_y = xyNorms[timeStepOfJthPoint][1];
				}
			}

			// Count the number of points whose x distance is less
			//  than or equal to eps_x, and whose y distance is less
			//  than or equal to eps_y
			int n_x = 0;
			int n_y = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyNorms[t2][0] <= eps_x) {
					n_x++;
				}
				if (xyNorms[t2][1] <= eps_y) {
					n_y++;
				}
			}
			sumNx += n_x;
			sumNy += n_y;
			// And take the digammas:
			double digammaNx = MathsUtils.digamma(n_x);
			double digammaNy = MathsUtils.digamma(n_y);
			sumDiGammas += digammaNx + digammaNy;

			if (returnLocals) {
				localMi[t-startTimePoint] = digammaK - invK - digammaNx - digammaNy + digammaN;
			}
		}

		if (debug) {
			Calendar rightNow2 = Calendar.getInstance();
			long endTime = rightNow2.getTimeInMillis();
			System.out.println("Subset " + startTimePoint + ":" +
					(startTimePoint + numTimePoints) + " Calculation time: " +
					((endTime - startTime)/1000.0) + " sec" );
		}
		
		// Select what to return:
		if (returnLocals) {
			return localMi;
		} else {
			return new double[] {sumDiGammas, sumNx, sumNy};
		}
	}

	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e - 1/k=%.3e + digamma(N=%d)=%.3e => %.3e",
				k, MathsUtils.digamma(k), 1.0/(double)k, N, MathsUtils.digamma(N),
				(MathsUtils.digamma(k) - 1.0/(double)k + MathsUtils.digamma(N)));
		return constants;
	}
}
