/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2017, Joseph T. Lizier, Ipek Oezdemir and Pedro Mediano
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

package infodynamics.measures.continuous.gaussian;

import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential O-information of a given multivariate
 *  <code>double[][]</code> set of
 *  observations (extending {@link MultiVariateInfoMeasureCalculatorGaussian}),
 *  assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution.</p>
 *
 * <p>Usage is as per the paradigm outlined for {@link MultiVariateInfoMeasureCalculatorCommon}.
 * </p>
 *
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Rosas, F., Mediano, P., Gastpar, M, Jensen, H.,
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.100.032305">"Quantifying high-order
 *   interdependencies via multivariate extensions of the mutual information"</a>,
 *   Physical Review E 100, (2019) 032305.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class OInfoCalculatorGaussian
  extends MultiVariateInfoMeasureCalculatorGaussian {

  /**
   * Constructor.
   */
  public OInfoCalculatorGaussian() {
    // Nothing to do
  }

  /**
   * {@inheritDoc}
   *
   * @return the average O-info in nats (not bits!)
   * @throws Exception if not sufficient data have been provided, or if the
   *   supplied covariance matrix is invalid.
   */
  public double computeAverageLocalOfObservations() throws Exception {

    if (covariance == null) {
      throw new Exception("Cannot calculate O-Info without having " +
          "a covariance either supplied or computed via setObservations()");
    }

    if (!isComputed) {
      double oinfo = (dimensions - 2)*Math.log(MatrixUtils.determinantSymmPosDefMatrix(covariance));
      for (int i = 0; i < dimensions; i++) {
        int[] idx = allExcept(i, dimensions);
        double[][] marginal_cov = MatrixUtils.selectRowsAndColumns(covariance, idx, idx);
        oinfo += Math.log(covariance[i][i]) - Math.log(MatrixUtils.determinantSymmPosDefMatrix(marginal_cov));
      }
      // This "0.5" comes from the entropy formula for Gaussians: h = 0.5*logdet(2*pi*e*Sigma)
      lastAverage = 0.5*oinfo;;
      isComputed = true;
    }

    return lastAverage;
  }

  /**
   * {@inheritDoc}
   *
   * @return the "time-series" of local O-info values in nats (not bits!)
   *   for the supplied states.
   * @throws Exception if not sufficient data have been provided, or if the
   *   supplied covariance matrix is invalid.
   */
  public double[] computeLocalUsingPreviousObservations(double[][] states) throws Exception {

		if ((means == null) || (covariance == null)) {
			throw new Exception("Cannot compute local values without having means " +
          "and covariance either supplied or computed via setObservations()");
		}

    EntropyCalculatorMultiVariateGaussian hCalc = new EntropyCalculatorMultiVariateGaussian();
    hCalc.initialise(dimensions);
    hCalc.setCovarianceAndMeans(covariance, means);
    double[] localValues = MatrixUtils.multiply(hCalc.computeLocalUsingPreviousObservations(states), dimensions - 2);

    for (int i = 0; i < dimensions; i++) {
      int[] idx = allExcept(i, dimensions);

      // Local entropy of this variable (i) only
      double[][] this_cov   = MatrixUtils.selectRowsAndColumns(covariance, i, 1, i, 1);
      double[]   this_means = MatrixUtils.select(means, i, 1);
      double[][] this_state = MatrixUtils.selectColumns(states, i, 1);

      hCalc.initialise(1);
      hCalc.setCovarianceAndMeans(this_cov, this_means);
      double[] thisLocals = hCalc.computeLocalUsingPreviousObservations(this_state);

      // Local entropy of the rest of the variables (0, ... i-1, i+1, ... D)
      double[][] rest_cov   = MatrixUtils.selectRowsAndColumns(covariance, idx, idx);
      double[]   rest_means = MatrixUtils.select(means, idx);
      double[][] rest_state = MatrixUtils.selectColumns(states, idx);

      hCalc.initialise(dimensions - 1);
      hCalc.setCovarianceAndMeans(rest_cov, rest_means);
      double[] restLocals = hCalc.computeLocalUsingPreviousObservations(rest_state);

      localValues = MatrixUtils.add(localValues, MatrixUtils.subtract(thisLocals, restLocals));

    }

    return localValues;

  }

}


