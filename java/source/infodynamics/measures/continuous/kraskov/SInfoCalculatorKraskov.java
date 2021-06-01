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

import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.NeighbourNodeData;
import infodynamics.utils.KdTree;
import infodynamics.utils.UnivariateNearestNeighbourSearcher;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

import java.util.PriorityQueue;
import java.util.Calendar;
import java.util.Random;
import java.util.Arrays;

/**
 * <p>Computes the differential S-information of a given multivariate
 * set of observations, using Kraskov-Stoegbauer-Grassberger (KSG) estimation
 * (see Kraskov et al., below).</p>
 *
 * <p>Usage is as per the paradigm outlined for
 * {@link MultiVariateInfoMeasureCalculatorCommon}.</p>
 *
 * <p>Finally, note that {@link Cloneable} is implemented allowing clone()
 *  to produce only an automatic shallow copy, which is fine
 *  for the statistical significance calculation it is intended for
 *  (none of the array 
 *  data will be changed there).
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>Rosas, F., Mediano, P., Gastpar, M., Jensen, H.,
 *   "Quantifying high-order effects via multivariate extensions of the
 *   mutual information".</li>
 *
 *  <li>Kraskov, A., Stoegbauer, H., Grassberger, P.,
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class SInfoCalculatorKraskov
  extends MultiVariateInfoMeasureCalculatorKraskov
  implements Cloneable { // See comments on clonability above


  protected double[] partialComputeFromObservations(
      int startTimePoint, int numTimePoints, boolean returnLocals) throws Exception {

    double startTime = Calendar.getInstance().getTimeInMillis();

    double[] localMi = null;
    if (returnLocals) {
      localMi = new double[numTimePoints];
    }
    
    // Constants:
    double dimensionsMinus1TimesDiGammaN = (double) (dimensions - 1) * digammaN;

    // Count the average number of points within eps_x for each marginal x of each point
    double totalSumF = 0.0;
        
    for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
      // Compute eps for this time step by
      //  finding the kth closest neighbour for point t:
      PriorityQueue<NeighbourNodeData> nnPQ =
          kdTreeJoint.findKNearestNeighbours(k, t, dynCorrExclTime);
      // First element in the PQ is the kth NN,
      //  and epsilon = kthNnData.distance
      NeighbourNodeData kthNnData = nnPQ.poll();

      // Distance to kth neighbour in joint space
      double eps = kthNnData.distance;

      double sumF = 0.0;

      sumF += (digammaK - digammaN);
      
      for (int d = 0; d < dimensions; d++) {
        int n_small = rangeSearchersInSmallMarginals[d].countPointsStrictlyWithinR(
                    t, eps, dynCorrExclTime);

        int n_big = rangeSearchersInBigMarginals[d].countPointsStrictlyWithinR(
                    t, eps, dynCorrExclTime);


        sumF -= (MathsUtils.digamma(n_big + 1) - digammaN)/dimensions;
        sumF -= (MathsUtils.digamma(n_small + 1) - digammaN)/dimensions;

        if (debug) {
          // Only tracking this for debugging purposes: 
          System.out.printf("t=%d, d=%d, n_small=%d, n_big=%d, sumF=%.3f%n",
            t, d, n_small, n_big, sumF);
        }

      }

      sumF *= dimensions;

      totalSumF += sumF;

      if (returnLocals) {
        localMi[t-startTimePoint] = sumF;
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
      double[] returnArray = new double[] {totalSumF/((double) totalObservations)};
      return returnArray;
    }

  }

}

