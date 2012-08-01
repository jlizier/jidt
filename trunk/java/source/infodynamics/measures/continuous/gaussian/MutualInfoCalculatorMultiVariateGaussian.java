package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * <p>Computes the differential mutual information of two given multivariate sets of
 *  observations,
 *  assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution.</p>
 *  
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct {@link #MutualInfoCalculatorMultiVariateLinearGaussian()}</li>
 *		<li>{@link #initialise(int, int)}</li>
 * 		<li>Provide the observations to the calculator using:
 * 			{@link #setObservations(double[][], double[][])}, or
 * 			{@link #setCovariance(double[][])}, or
 * 			a sequence of:
 * 			{@link #startAddObservations()},
 *          multiple calls to either {@link #addObservations(double[][], double[][])}
 *          or @{link {@link #addObservations(double[][], double[][], int, int)}, and then
 *          @{link {@link #finaliseAddObservations()}.</li>
 * 		<li>Compute the required information-theoretic results, primarily:
 * 			@{link #computeAverageLocalOfObservations()} to return the average differential
 *          entropy based on either the set variance or the variance of
 *          the supplied observations; or other calls to compute
 *          local values or statistical significance.</li>
 * 	</ol>
 * </p>
 * 
 * @see Differential entropy for Gaussian random variables defined at 
 *      {@link http://mathworld.wolfram.com/DifferentialEntropy.html}
 * @author Joseph Lizier joseph.lizier_at_gmail.com
 *
 */
public class MutualInfoCalculatorMultiVariateGaussian implements
		MutualInfoCalculatorMultiVariate {

	/**
	 * Covariance matrix of the most recently supplied observations
	 */
	protected double[][] covariance;

	/**
	 * The set of source observations, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link addObservations(double[][], double[][])} functions.
	 */
	protected double[][] sourceObservations;

	/**
	 * The set of destination observations, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link addObservations(double[][], double[][])} functions.
	 */
	protected double[][] destObservations;

	/**
	 * Number of dimenions for each of our multivariate data sets
	 */
	protected int dimensionsDest;
	protected int dimensionsSource;
	
	protected double lastAverage;
	
	protected boolean debug;

	public MutualInfoCalculatorMultiVariateGaussian() {
		// Nothing to do
	}
	
	/**
	 * Clear any previously supplied probability distributions and prepare
	 * the calculator to be used again.
	 * 
	 * @param sourceDimensions number of joing variables in the source
	 * @param destDimensions number of joing variables in the destination
	 */
	public void initialise(int sourceDimensions, int destDimensions) {
		covariance = null;
		sourceObservations = null;
		destObservations = null;
		dimensionsSource = sourceDimensions;
		dimensionsDest = destDimensions;
	}

	/**
	 * Provide the complete set of observations to use to compute the
	 *  mutual information.
	 * One cannot use the {@link addObservations(double[][], double[][])}
	 *  style methods after this without calling initialise again first.
	 * 
	 */
	public void setObservations(double[][] source, double[][] destination)
			throws Exception {
		sourceObservations = source;
		destObservations = destination;
		covariance = MatrixUtils.covarianceMatrix(source, destination);
		// Check that the observations was of the correct number of dimensions:
		//  (done afterwards since the covariance matrix computation checks that
		//  all rows had the right number of columns
		if (covariance.length != dimensionsSource + dimensionsDest) {
			throw new RuntimeException("Supplied observations do not match initialised number of dimensions");
		}
	}

	public void addObservations(double[][] source, double[][] destination)
			throws Exception {
		// TODO implement these addObservations style functions.
		// This will not be hard to implement - see the implementation
		//  for TE in TransferEntropyCommon. It might be useful
		//  to have a MutualInfoCommon which pulls the same functionality
		//  together for the MI calculators anyway.
		throw new RuntimeException("Not implemented yet");
	}

	public void addObservations(double[][] source, double[][] destination,
			int startTime, int numTimeSteps) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void startAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void finaliseAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	/**
	 * Set the covariance of the distribution for which we will compute the
	 *  mutual information.
	 * 
	 * @param covariance covariance matrix of the source and destination
	 *  variables, considered together (variable indices start with the source
	 *  and continue into the destination).
	 */
	public void setCovariance(double[][] covariance) throws Exception {
		sourceObservations = null;
		destObservations = null;
		// Make sure the supplied covariance matrix is square:
		int rows = covariance.length;
		if (rows != dimensionsSource + dimensionsDest) {
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
	 * <p>Here we compute the mutual information from the joint entropies
	 *  of the source variables (H_s), destination variables (H_d), and all variables
	 *  taken together (H_sd), giving MI = H_s + H_d - H_sd.
	 *  We assume that the recorded estimation of the
	 *  covariance is correct (i.e. we will not make a bias correction for limited
	 *  observations here).</p>
	 * 
	 * @return the mutual information of the previously provided observations or from the
	 *  supplied covariance matrix.
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		try {
			int[] sourceIndicesInCovariance = MatrixUtils.range(0, dimensionsSource - 1);
			int[] destIndicesInCovariance = MatrixUtils.range(dimensionsSource,
						dimensionsSource + dimensionsDest - 1);
			double[][] sourceCovariance =
				MatrixUtils.selectRowsAndColumns(covariance, 
						sourceIndicesInCovariance, sourceIndicesInCovariance);
			double[][] destCovariance =
				MatrixUtils.selectRowsAndColumns(covariance, 
					destIndicesInCovariance, destIndicesInCovariance);
			double sourceEntropy = 0.5 *
				Math.log(Math.abs(MatrixUtils.determinant(sourceCovariance)));
			double destEntropy = 0.5 *
				Math.log(Math.abs(MatrixUtils.determinant(destCovariance)));
			double jointEntropy = 0.5 *
				Math.log(Math.abs(MatrixUtils.determinant(covariance)));
			lastAverage = sourceEntropy + destEntropy - jointEntropy;
			return lastAverage;
		} catch (Exception e) {
			// Should not happen, since we check the validity of the supplied
			//  matrix beforehand; so we'll throw an Error in this case
			throw new Error(e);
		}
	}

	public double[] computeLocalOfPreviousObservations() throws Exception {
		// TODO Implement me
		throw new RuntimeException("Not implemented yet");
	}

	/**
	 * <p>Compute the statistical significance of the mutual information 
	 *  result analytically, without creating a distribution
	 *  under the null hypothesis by bootstrapping.</p>
	 *
	 * <p>Brillinger (see reference below) shows that under the null hypothesis
	 *  of no source-destination relationship, the MI for two
	 *  Gaussian distributions follows a chi-square distribution with
	 *  degrees of freedom equal to the product of the number of variables
	 *  in each joint variable.</p>
	 *
	 * @return MeasurementDistribution object with only the 
	 *  {@link EmpiricalMeasurementDistribution#actualValue} and
	 *  {@link EmpiricalMeasurementDistribution#pValue} fields filled out.
	 *  This object contains the proportion of MI scores from the distribution
	 *  which have higher or equal MIs to ours.
	 *  
	 * @see Brillinger, "Some data analyses using mutual information",
	 * {@link http://www.stat.berkeley.edu/~brill/Papers/MIBJPS.pdf}
	 * @see Cheng et al., "Data Information in Contingency Tables: A
	 *  Fallacy of Hierarchical Loglinear Models",
	 *  {@link http://www.jds-online.com/file_download/112/JDS-369.pdf}
	 * @see Barnett and Bossomaier, "Transfer Entropy as a Log-likelihood Ratio" 
	 *  {@link http://arxiv.org/abs/1205.6339}
	 */
	public ChiSquareMeasurementDistribution computeSignificance() {
		// TODO Check that the null distribution actually follows chi with
		//  these degrees of freedom
		return new ChiSquareMeasurementDistribution(lastAverage,
				dimensionsSource * dimensionsDest);
	}
	
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		// TODO Implement me
		throw new RuntimeException("Not implemented yet");
	}

	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings)
			throws Exception {
		// TODO Implement me
		throw new RuntimeException("Not implemented yet");
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

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	/**
	 * @return the previously computed average mutual information
	 */
	public double getLastAverage() {
		return lastAverage;
	}

	/**
	 * @return the number of previously supplied observations for which
	 *  the mutual information will be / was computed.
	 */
	public int getNumObservations() {
		return destObservations.length;
	}

	public double computeAverageLocalOfObservations(int[] newOrdering)
			throws Exception {
		// TODO Implement me
		throw new RuntimeException("Not implemented yet");
	}

	public double[] computeLocalUsingPreviousObservations(double[][] states1,
			double[][] states2) throws Exception {
		// TODO Implement me
		throw new RuntimeException("Not implemented yet");
	}

}
