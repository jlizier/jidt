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

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential conditional mutual information of two multivariate
 *  <code>double[][]</code> sets of observations, conditioned on another
 *  (implementing {@link ConditionalMutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see references below)
 *  <b>algorithm 2</b>.
 *  Most of the functionality is defined by the parent class 
 *  {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}.</p>
 *
 * <p>Crucially, the calculation is performed by examining
 * neighbours in the full joint space
 * rather than two MI calculators.
 * This is roughly as specified by Frenzel and Pompe (who only
 * specified this for algorithm 1), but mathematically
 * adapted to algorithm 2 by Wibral et al. (see below).</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link ConditionalMutualInfoCalculatorMultiVariate},
 * and expanded on in {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}.
 * </p>
 *  
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>M. Wibral, R. Vicente, and M. Lindner, <a href="http://dx.doi.org/10.1007/978-3-642-54474-3_1">
 *  "Transfer Entropy in Neuroscience"</a>,
 *  in "Directed Information Measures in Neuroscience",
 *  Understanding Complex Systems series, edited by M. Wibral, R. Vicente, and J. T. Lizier
 *  (Springer, Berlin/Heidelberg, 2014) pp. 3--36.</li>
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
 * @author Ipek Özdemir
 */
public class ConditionalMutualInfoCalculatorMultiVariateKraskov2
	extends ConditionalMutualInfoCalculatorMultiVariateKraskov {

	protected static final int JOINT_NORM_VAL_COLUMN = 0;
	protected static final int JOINT_NORM_TIMESTEP_COLUMN = 1;

	// Multiplier used in hueristic for determining whether to use a linear search
	//  for min kth element or a binary search.
	protected static final double CUTOFF_MULTIPLIER = 1.5;
	
	@Override
	protected double[] partialComputeFromObservations(int startTimePoint,
			int numTimePoints, boolean returnLocals) throws Exception {
		
		double startTime = Calendar.getInstance().getTimeInMillis();

		int N = var1Observations.length; // number of observations
		int cutoffForKthMinLinear = (int) (CUTOFF_MULTIPLIER * Math.log(N) / Math.log(2.0));

		double[] localCondMi = null;
		if (returnLocals) {
			localCondMi = new double[numTimePoints];
		}
		
		// Constants:
		double digammaK = MathsUtils.digamma(k);
		double twoOnK = 2.0 / (double) k;
		
		// Count the average number of points within eps_x and eps_y
		double sumDiGammas = 0;
		double sumNxz = 0;
		double sumNyz = 0;
		double sumNz = 0;
		double sumInverseCountInJointYZ = 0;
		double sumInverseCountInJointXZ = 0;

		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {

			// Compute eps_x and eps_y and eps_z for this time step:
			//  First get x and y and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xyzNorms = normCalculator.computeNorms(var1Observations, var2Observations, condObservations, t);
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][JOINT_NORM_VAL_COLUMN] = Math.max(xyzNorms[t2][0],
							Math.max(xyzNorms[t2][1], xyzNorms[t2][2]));
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][JOINT_NORM_TIMESTEP_COLUMN] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_y = 0.0;
			double eps_z = 0.0;
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
			// Find eps_{x,y,z} as the maximum x and y norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xyzNorms[timeStepOfJthPoint][0] > eps_x) {
					eps_x = xyzNorms[timeStepOfJthPoint][0];
				}
				if (xyzNorms[timeStepOfJthPoint][1] > eps_y) {
					eps_y = xyzNorms[timeStepOfJthPoint][1];
				}
				if (xyzNorms[timeStepOfJthPoint][2] > eps_z) {
					eps_z = xyzNorms[timeStepOfJthPoint][2];
				}
			}

			// Count the number of points whose z distance is less
			//  than or equal to eps_z, and whose x and z distances are less
			//  than or equal to eps_z and eps_x, and whose y and z distance are less
			//  than or equal to eps_z and eps_y:
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xyzNorms[t2][2] <= eps_z) {
					n_z++;
					if (xyzNorms[t2][0] <= eps_x) {
						n_xz++;
					}
					if (xyzNorms[t2][1] <= eps_y) {
						n_yz++;
					}
				}
			}
			sumNxz += n_xz;
			sumNyz += n_yz;
			sumNz += n_z;
			// And take the digammas:
			double digammaNxz = MathsUtils.digamma(n_xz);
			double digammaNyz = MathsUtils.digamma(n_yz);
			double digammaNz = MathsUtils.digamma(n_z);
			double invN_xz = 1.0/(double) n_xz;
			double invN_yz = 1.0/(double) n_yz;
			sumInverseCountInJointXZ += invN_xz;
			sumInverseCountInJointYZ += invN_yz;
			double contributionDigammas = digammaNz - digammaNxz - digammaNyz;
			sumDiGammas += contributionDigammas;

			if (returnLocals) {
				localCondMi[t-startTimePoint] = digammaK - twoOnK + contributionDigammas + invN_xz + invN_yz;
				if (debug) {
					// Only tracking this for debugging purposes: 
					System.out.printf("t=%d, n_xz=%d, n_yz=%d, n_z=%d, 1/n_yz=%.3f, 1/n_xz=%.3f, local=%.4f\n",
						t, n_xz, n_yz, n_z, invN_yz, invN_xz, localCondMi[t-startTimePoint]);
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
		
		// Select what to return:
		if (returnLocals) {
			return localCondMi;
		} else {
			return new double[] {sumDiGammas, sumNxz, sumNyz, sumNz,
					sumInverseCountInJointXZ, sumInverseCountInJointYZ};
		}		
	}

	@Override
	public String printConstants(int N) throws Exception {
		String constants = String.format("digamma(k=%d)=%.3e",
				k, MathsUtils.digamma(k));
		return constants;
	}
}
