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

import infodynamics.measures.continuous.MultiVariateInfoMeasureCalculatorCommon;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Base class with common functionality for child class implementations of
 * multivariate information measures on a given multivariate
 *  <code>double[][]</code> set of
 *  observations (extending {@link MultiVariateInfoMeasureCalculatorCommon}),
 *  assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution.</p>
 *
 * <p>Usage is as per the paradigm outlined for {@link MultiVariateInfoMeasureCalculatorCommon},
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
 * 	<li>Rosas, F., Mediano, P., Gastpar, M, Jensen, H.,
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.100.032305">"Quantifying high-order
 *   interdependencies via multivariate extensions of the mutual information"</a>,
 *   Physical Review E 100, (2019) 032305.</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public abstract class MultiVariateInfoMeasureCalculatorGaussian
  extends MultiVariateInfoMeasureCalculatorCommon {

  /**
   * Covariance of the system. Can be calculated from supplied observations
   * or supplied directly by the user.
   */
  double[][] covariance = null;

  /**
   * Means of the system. Can be calculated from supplied observations
   * or supplied directly by the user.
   */
  double[] means = null;

  /**
   * Whether the current covariance matrix has been determined from data or
   * supplied directly. This changes the approach to local measures and
   * significance testing.
   */
  boolean covFromObservations;

  @Override
  public void finaliseAddObservations() throws Exception {
    super.finaliseAddObservations();

    this.means = MatrixUtils.means(observations);
    setCovariance(MatrixUtils.covarianceMatrix(observations), true);

    return;
  }

  /**
   * <p>Set the covariance of the distribution for which we will compute the
   *  measure directly, without supplying observations.</p>
   *
   *  <p>See {@link #setCovariance(double[][], boolean)}.
   *
   * @param covariance covariance matrix of the system
   * @throws Exception for covariance matrix not matching the expected dimensions,
   *  being non-square, asymmetric or non-positive definite
   */
  public void setCovariance(double[][] covariance) throws Exception {
    setCovariance(covariance, false);
  }

  /**
   * <p>Set the covariance of the distribution for which we will compute the
   *  measure.</p>
   *
   * <p>This is an alternative to sequences of calls to {@link #setObservations(double[][])} or
   * {@link #addObservations(double[][])} etc.
   * Note that without setting any observations, you cannot later
   *  call {@link #computeLocalOfPreviousObservations()}.</p>
   *
   * @param covariance covariance matrix of the system
   * @param means mean of the system
   */
  public void setCovarianceAndMeans(double[][] covariance, double[] means)
      throws Exception {
    this.means = means;
    setCovariance(covariance, false);
  }

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  measure.</p>
	 *
	 * <p>This is an alternative to sequences of calls to {@link #setObservations(double[][])} or
	 * {@link #addObservations(double[][])} etc.
	 * Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][])}.</p>
	 *
	 * @param covariance covariance matrix of the system
	 * @param covFromObservations whether the covariance matrix
	 *  was determined internally from observations or not
	 * @throws Exception for covariance matrix not matching the expected dimensions,
	 *  being non-square, asymmetric or non-positive definite
	 */
  public void setCovariance(double[][] cov, boolean covFromObservations) throws Exception {

    if (!covFromObservations) {
      // Make sure we're not keeping any observations
      observations = null;
    }
    if (cov.length != dimensions) {
      throw new Exception("Supplied covariance matrix does not match initialised number of dimensions");
    }
    if (cov.length != cov[0].length) {
      throw new Exception("Covariance matrices must be square");
    }

    this.covFromObservations = covFromObservations;
    this.covariance = cov;

  }

}

