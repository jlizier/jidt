package infodynamics.measures.continuous.gaussian;

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
public class EntropyCalculatorMultiVariateGaussian
	implements EntropyCalculatorMultiVariate, Cloneable {

	/**
	 * Covariance matrix of the most recently supplied observations
	 */
	protected double[][] covariance;
	
	/**
	 * Means of the most recently supplied observations (source variables
	 *  listed first, destination variables second).
	 */
	protected double[] means;
	
	/**
	 * The set of observations, retained in case the user wants to retrieve the local
	 *  entropy values of these
	 */
	protected double[][] observations;
	
	/**
	 * Number of dimenions for our multivariate data
	 */
	protected int dimensions;
	
	/**
	 * Determinant of the covariance matrix; stored to save computation time
	 */
	protected double detCovariance;

	protected double lastAverage;
	
	protected boolean debug;
	
	/**
	 * Constructor
	 */
	public EntropyCalculatorMultiVariateGaussian() {
		// Nothing to do
	}
	
	/**
	 * Initialise the calculator ready for reuse
	 */
	public void initialise(int dimensions) {
		covariance = null;
		means = null;
		observations = null;
		this.dimensions = dimensions;
		detCovariance = 0;
	}

	/**
	 * Provide the multivariate observations from which to compute the entropy
	 * 
	 * @param observations the observations to compute the entropy from.
	 *        First index is time, second index is variable number.
	 */
	public void setObservations(double[][] observations) {
		this.observations = observations;
		means = MatrixUtils.means(observations);
		covariance = MatrixUtils.covarianceMatrix(observations, means);
		detCovariance = 0;
		// Check that the observations was of the correct number of dimensions:
		//  (done afterwards since the covariance matrix computation checks that
		//  all rows had the right number of columns
		if (covariance.length != dimensions) {
			throw new RuntimeException("Supplied observations does not match initialised number of dimensions");
		}
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  entropy.</p>
	 * 
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}.</p>
	 *  
	 * @param covariance covariance matrix
	 */
	public void setCovariance(double[][] covariance) throws Exception {
		detCovariance = 0;
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
	 * <p>Set the covariance and mean of the distribution for which we will compute the
	 *  entropy.</p>
	 * 
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}.</p>
	 * 
	 * @param covariance covariance matrix of the variables
	 * @param means mean of the variables
	 * @throws Exception where the dimensions of the covariance or means are not correct
	 */
	public void setCovarianceAndMeans(double[][] covariance, double[] means) throws Exception {
		if (means.length != dimensions) {
			throw new Exception("Supplied mean matrix does not match initialised number of dimensions");
		}
		this.means = means;
		setCovariance(covariance);
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
			detCovariance = MatrixUtils.determinant(covariance);
			lastAverage =  0.5 * (dimensions* (1 + Math.log(2.0*Math.PI)) +
				Math.log(Math.abs(detCovariance)));
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

	/**
	 * Compute the local entropy for each of the given supplied
	 *  observations
	 * 
	 * @param states joint vectors of observations; first index
	 *  is observation number, second is variable number.
	 * @return an array of local values
	 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Multivariate normal distribution on Wikipedia</a>
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] states)
			throws Exception {
		if (means == null) {
			throw new Exception("Cannot compute local values without having means either supplied or computed via setObservations()");
		}
		
		// Check that the covariance matrix was positive definite:
		//  it is known (i think because it is symmetric) that the covariance
		//  matrix will be positive definite unless one variable is an exact
		//  linear combination of the others (ref: wikipedia, above).
		// We can check for linear dependence by ensuring that the determinant
		//  of the covariance matrix is non-zero:
		if (detCovariance == 0) {
			// We need to check:
			detCovariance = MatrixUtils.determinant(covariance); 
			if (detCovariance == 0) {
				throw new Exception("Covariance matrix is not positive definite");
			}
		}
		// Now we are clear to take the matrix inverse (via Cholesky decomposition,
		//  since we have a symmetric positive definite matrix):
		double[][] invCovariance = MatrixUtils.invertSymmPosDefMatrix(covariance);
		
		// If we have a time delay, slide the local values
		double[] localValues = new double[states.length];
		for (int t = 0; t < states.length; t++) {
			double[] deviationsFromMean =
					MatrixUtils.subtract(states[t], means);
			// Computing PDF
			// (see the PDF defined at the wikipedia page referenced in the method header)
			double jointExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(deviationsFromMean,
							invCovariance),
					deviationsFromMean);
			double pJoint = Math.pow(2.0 * Math.PI, -(double) dimensions / 2.0) *
					Math.exp(-0.5 * jointExpArg) /
					Math.sqrt(detCovariance);
			localValues[t] = - Math.log(pJoint);
		}
		
		// Don't set average if this was the previously supplied observations,
		//  since it won't be the same as what would have been computed 
		//  analytically.
		
		return localValues;
	}

	public double[] computeLocalOfPreviousObservations() throws Exception {
		if (observations == null) {
			throw new Exception("Cannot compute local values since no observations were supplied");
		}
		return computeLocalUsingPreviousObservations(observations);
	}

	/**
	 * Provide an implementation of the clone() method.
	 * This does not deeply copy all of the underlying data, just providing
	 *  a copy of the references to it all.
	 * This is enough to protect the integrity of the calculator
	 *  however if the clone is supplied different data (though the 
	 *  clone should not alter the data).
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}
}
