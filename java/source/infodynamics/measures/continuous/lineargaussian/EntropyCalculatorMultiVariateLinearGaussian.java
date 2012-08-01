package infodynamics.measures.continuous.lineargaussian;

import infodynamics.measures.continuous.EntropyCalculatorMultiVariate;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential entropy of a given multivariate set of observations,
 *  assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution.</p>
 *  
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct</li>
 *		<li>initialise()</li>
 * 		<li>setObservations(), or setCovariance().</li> 
 * 		<li>computeAverageLocalOfObservations() to return the average differential
 *          entropy based on either the set variance or the variance of
 *          the supplied observations, or computeLocalUsingPrevious</li>
 * 	</ol>
 * </p>
 * 
 * @see Differential entropy for Gaussian random variables defined at 
 *      {@link http://mathworld.wolfram.com/DifferentialEntropy.html}
 * @author Joseph Lizier joseph.lizier_at_gmail.com
 *
 */
public class EntropyCalculatorMultiVariateLinearGaussian implements EntropyCalculatorMultiVariate {

	/**
	 * Covariance matrix of the most recently supplied observations
	 */
	protected double[][] covariance;
	
	/**
	 * The set of observations, retained in case the user wants to retrieve the local
	 *  entropy values of these
	 */
	protected double[][] observations;
	
	/**
	 * Number of dimenions for our multivariate data
	 */
	protected int dimensions;
	
	protected double lastAverage;
	
	protected boolean debug;
	
	/**
	 * Constructor
	 */
	public EntropyCalculatorMultiVariateLinearGaussian() {
		// Nothing to do
	}
	
	/**
	 * Initialise the calculator ready for reuse
	 */
	public void initialise(int dimensions) {
		covariance = null;
		observations = null;
		this.dimensions = dimensions;
	}

	/**
	 * Provide the multivariate observations from which to compute the entropy
	 * 
	 * @param observations the observations to compute the entropy from.
	 *        First index is time, second index is variable number.
	 */
	public void setObservations(double[][] observations) {
		this.observations = observations;
		covariance = MatrixUtils.covarianceMatrix(observations);
		// Check that the observations was of the correct number of dimensions:
		//  (done afterwards since the covariance matrix computation checks that
		//  all rows had the right number of columns
		if (covariance.length != dimensions) {
			throw new RuntimeException("Supplied observations does not match initialised number of dimensions");
		}
	}

	/**
	 * Set the covariance of the distribution for which we will compute the
	 *  entropy.
	 * 
	 * @param covariance covariance matrix
	 */
	public void setCovariance(double[][] covariance) throws Exception {
		observations = null;
		// Make sure the supplied covariance matrix is square:
		int rows = covariance.length;
		if (rows != dimensions) {
			throw new Exception("Supplied covariance matrix does not match initialised number of dimensions");
		}
		for (int r = 0; r < rows; r++) {
			if (covariance[r].length != rows) {
				throw new Exception("Cannot compute the determinant of a non-square matrix");
			}
		}

		this.covariance = covariance;
	}
	
	/**
	 * <p>The joint entropy for a multivariate Gaussian-distribution of dimension n
	 *  with covariance matrix C is 0.5*\log_e{(2*pi*e)^n*|det(C)|},
	 *  where det() is the matrix determinant of C.</p>
	 * 
	 * <p>Here we compute the joint entropy assuming that the recorded estimation of the
	 *  covariance is correct (i.e. we will not make a bias correction for limited
	 *  observations here).</p>
	 * 
	 * @return the joint entropy of the previously provided observations or from the
	 *  supplied covariance matrix.
	 */
	public double computeAverageLocalOfObservations() {
		try {
			lastAverage =  0.5 * (dimensions* (1 + Math.log(2.0*Math.PI)) +
				Math.log(Math.abs(MatrixUtils.determinant(covariance))));
			return lastAverage;
		} catch (Exception e) {
			// Should not happen, since we check the validity of the supplied
			//  matrix beforehand; so we'll throw an Error in this case
			throw new Error(e);
		}
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	/**
	 * <p>Set the given property to the given value.</p>
	 * 
	 * <p>There are currently no properties to set for this calculator</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		// No properties to set here
	}

	/**
	 * @return the lastAverage
	 */
	public double getLastAverage() {
		return lastAverage;
	}

	public double[] computeLocalUsingPreviousObservations(double[][] states)
			throws Exception {
		// TODO Implement me
		throw new RuntimeException("Not implemented yet");
	}

	public double[] computeLocalOfPreviousObservations() {
		
		// TODO Implement this function
		if (true)
			throw new RuntimeException("Not implemented yet");
		
		if (observations == null) {
			throw new RuntimeException("Cannot compute local values since no observations were supplied");
		}
		double[] localEntropy = new double[observations.length];
		for (int t=0; t < observations.length; t++) {
			// Compute the probability for the given observation, based on 
			//  the assumption of a multivariate Gaussian PDF:
				
		}
		return null;
	}
}
