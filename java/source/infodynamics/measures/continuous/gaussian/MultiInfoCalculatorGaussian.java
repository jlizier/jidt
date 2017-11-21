package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.MultiInfoCalculatorCommon;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential multi-information of a given multivariate
 *  <code>double[][]</code> set of
 *  observations (implementing {@link MultiInfoCalculator}),
 *  assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution.</p>
 *
 * <p>Usage is as per the paradigm outlined for {@link MultiInfoCalculator},
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
 *  <li>T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991).</li>
 * </ul>
 *
 * @author Pedro A.M. Mediano (<a href="pmediano at imperial.ac.uk">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class MultiInfoCalculatorGaussian
  extends MultiInfoCalculatorCommon {

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


  /**
   * Constructor.
   */
  public MultiInfoCalculatorGaussian() {
    // Nothing to do
  }

  @Override
  public void finaliseAddObservations() throws Exception {
    super.finaliseAddObservations();

    setCovariance(MatrixUtils.covarianceMatrix(observations), true);

    return;
  }

  /**
   * <p>Set the covariance of the distribution for which we will compute the
   *  multi-information directly, without supplying observations.</p>
   *  
   *  <p>See {@link #setCovariance(double[][], boolean)}.
   * 
   * @param covariance covariance matrix of the system
   * @throws Exception for covariance matrix not matching the expected dimensions,
   *  being non-square, asymmetric or non-positive definite
   */
  public void setCovariance(double[][] cov) throws Exception {
    setCovariance(cov, false);
  }

  /**
   * <p>Set the covariance of the distribution for which we will compute the
   *  multi-information.</p>
   * 
   * <p>This is an alternative to sequences of calls to {@link #setObservations(double[][])} or
   * {@link #addObservations(double[][])} etc.
   * Note that without setting any observations, you cannot later
   *  call {@link #computeLocalOfPreviousObservations()}.</p>
   * 
   * @param covariance covariance matrix of the system
   *  variables, considered together (variable indices start with the source
   *  and continue into the destination).
   * @param means mean of the system
   * @param numObservations the number of observations that the mean and covariance
   *  were determined from. This is used for later significance calculations
   */
  public void setCovarianceAndMeans(double[][] covariance, double[] means)
      throws Exception {
    this.means = means;
    setCovariance(covariance, false);
  }

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  multi-information.</p>
	 *  
	 * <p>This is an alternative to sequences of calls to {@link #setObservations(double[][])} or
	 * {@link #addObservations(double[][])} etc.
	 * Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][])}.</p>
	 * 
	 * @param covariance covariance matrix of the system
	 * @param determinedFromObservations whether the covariance matrix
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

  /**
   * {@inheritDoc} 
   * 
   * @return the average multi-info in nats (not bits!)
   */
  public double computeAverageLocalOfObservations() throws Exception {

    if (covariance == null) {
      throw new Exception("Cannot calculate multi-information without having " +
          "a covariance either supplied or computed via setObservations()");
    }

    if (!miComputed) {
      double mi = - Math.log(MatrixUtils.determinantSymmPosDefMatrix(covariance));
      for (int i = 0; i < dimensions; i++) {
        mi += Math.log(covariance[i][i]);
      }
      lastAverage = 0.5*mi;;
    }

    return lastAverage;
  }

  /**
   * {@inheritDoc}
   * 
   * @return the "time-series" of local multi-info values in nats (not bits!)
   * @throws Exception
   */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		// Cannot do if destObservations haven't been set
		if (observations == null) {
			throw new Exception("Cannot compute local values of previous observations " +
					"if they have not been set!");
		}
		
		return computeLocalUsingPreviousObservations(observations);
	}

  /**
   * {@inheritDoc}
   * 
   * @return the "time-series" of local multi-info values in nats (not bits!)
   *   for the supplied states.
   * @throws Exception
   */
  public double[] computeLocalUsingPreviousObservations(double[][] states) throws Exception {

		if ((means == null) || (covariance == null)) {
			throw new Exception("Cannot compute local values without having means " +
          "and covariance either supplied or computed via setObservations()");
		}

    double[][] L = MatrixUtils.CholeskyDecomposition(covariance);
    double[][] invCovariance = MatrixUtils.solveViaCholeskyResult(L, MatrixUtils.identityMatrix(L.length));
    double detCovariance = MatrixUtils.determinantViaCholeskyResult(L);
		double[] localValues = new double[states.length];

    for (int t = 0; t < states.length; t++) {

      double[] deviationsFromMean = MatrixUtils.subtract(states[t], means);
			double jointExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(deviationsFromMean,
							invCovariance),
					deviationsFromMean);
			double adjustedPJoint = Math.exp(-0.5 * jointExpArg) /
					Math.sqrt(detCovariance);

      double localValue = Math.log(adjustedPJoint);

      for (int i = 0; i < dimensions; i++) {

        double adjustedPVar = Math.exp(-0.5*(states[t][i] - means[i])*
            (states[t][i] - means[i])/covariance[i][i])/Math.sqrt(covariance[i][i]);
        localValue -= Math.log(adjustedPVar);
      }
			// Returning results in nats:
			localValues[t] = localValue;			

    }

    return localValues;
  }

}

