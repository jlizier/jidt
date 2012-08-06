package infodynamics.measures.continuous.gaussian;

import java.util.Iterator;

import infodynamics.measures.continuous.AnalyticNullDistributionComputer;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.MutualInfoMultiVariateCommon;
import infodynamics.utils.ChiSquareMeasurementDistribution;
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
public class MutualInfoCalculatorMultiVariateGaussian 
		extends MutualInfoMultiVariateCommon
		implements MutualInfoCalculatorMultiVariate, AnalyticNullDistributionComputer {

	/**
	 * Covariance matrix of the most recently supplied observations.
	 * Is a matrix [C_ss, C_sd; C_ds, C_dd], where C_ss is the covariance
	 * matrix of the source observations, C_dd is the covariance matrix
	 * of the destination observations, and C_sd and C_ds are the covariances
	 * of source to destination and destination to source observations.
	 */
	protected double[][] covariance;
	
	/**
	 * Means of the most recently supplied observations (source variables
	 *  listed first, destination variables second).
	 */
	protected double[] means;

	protected double detCovariance;
	
	public MutualInfoCalculatorMultiVariateGaussian() {
		// Nothing to do
	}
	
	/**
	 * Clear any previously supplied probability distributions and prepare
	 * the calculator to be used again.
	 * 
	 * @param sourceDimensions number of joint variables in the source
	 * @param destDimensions number of joint variables in the destination
	 */
	public void initialise(int sourceDimensions, int destDimensions) {
		super.initialise(sourceDimensions, destDimensions);
		covariance = null;
		means = null;
		sourceObservations = null;
		destObservations = null;
		detCovariance = 0;
	}

	public void setObservations(double[][] source, double[][] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		// TODO implement me
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		// TODO implement me
		throw new RuntimeException("Not implemented yet");
	}

	/**
	 * Finalise the addition of multiple observation sets.
	 * 
	 */
	public void finaliseAddObservations() {

		// Get the observations properly stored in the sourceObservations[][] and
		//  destObservations[][] arrays.
		super.finaliseAddObservations();

		// Store the means of each variable (useful for local values later)
		means = new double[dimensionsSource + dimensionsDest];
		double[] sourceMeans = MatrixUtils.means(sourceObservations);
		double[] destMeans = MatrixUtils.means(destObservations);
		System.arraycopy(sourceMeans, 0, means, 0, dimensionsSource);
		System.arraycopy(destMeans, 0, means, dimensionsSource, dimensionsDest);
		
		// Store the covariances of the variables
		covariance = MatrixUtils.covarianceMatrix(sourceObservations,
				destObservations);
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  mutual information.</p>
	 *  
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][])}.</p>
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
				throw new Exception("Covariance matrix must be square");
			}
			// Check that is is symmetric
			for (int c = 0; c < r; c++) {
				if (covariance[r][c] != covariance[c][r]) {
					throw new Exception("Covariance matrix is not symmetric!");
				}
			}
		}

		this.covariance = covariance;
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  mutual information.</p>
	 * 
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}.</p>
	 * 
	 * @param covariance covariance matrix of the source and destination
	 *  variables, considered together (variable indices start with the source
	 *  and continue into the destination).
	 */
	public void setCovarianceAndMeans(double[][] covariance, double[] means) throws Exception {
		this.means = means;
		setCovariance(covariance);
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
	 *  supplied covariance matrix, in nats (not bits!).
	 *  Returns NaN if any of the determinants are zero
	 *  (because this will make the denominator of the log zero)
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		int[] sourceIndicesInCovariance = MatrixUtils.range(0, dimensionsSource - 1);
		int[] destIndicesInCovariance = MatrixUtils.range(dimensionsSource,
					dimensionsSource + dimensionsDest - 1);
		double[][] sourceCovariance =
			MatrixUtils.selectRowsAndColumns(covariance, 
					sourceIndicesInCovariance, sourceIndicesInCovariance);
		double[][] destCovariance =
			MatrixUtils.selectRowsAndColumns(covariance, 
				destIndicesInCovariance, destIndicesInCovariance);
		detCovariance = MatrixUtils.determinant(covariance);
		lastAverage = 0.5 * Math.log(Math.abs(
					MatrixUtils.determinant(sourceCovariance) *
						MatrixUtils.determinant(destCovariance) /
						detCovariance));
		return lastAverage;
	}

	/**
	 * <p>Compute the local or pointwise mutual information for each of the previously
	 * supplied observations</p>
	 * 
	 * @return array of the local values in nats (not bits!)
	 * @see {@link http://en.wikipedia.org/wiki/Positive-definite_matrix}
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		// Cannot do if destObservations haven't been set
		if (destObservations == null) {
			throw new Exception("Cannot compute local values of previous observations " +
					"if they have not been set!");
		}
		
		return computeLocalUsingPreviousObservations(sourceObservations,
				destObservations, true);
		
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
	 * @return ChiSquareMeasurementDistribution object 
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
	 * @return the number of previously supplied observations for which
	 *  the mutual information will be / was computed.
	 */
	public int getNumObservations() throws Exception {
		if (destObservations == null) {
			throw new Exception("Cannot return number of observations because either " +
					"this calculator has not had observations supplied or " +
					"the user supplied the covariance matrix instead of observations");
		}
		return destObservations.length;
	}

	/**
	 * Compute the mutual information if the second variable were ordered as per the ordering
	 *  specified in newOrdering
	 * 
	 * @param newOrdering
	 * @return
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[] newOrdering)
			throws Exception {
		// TODO Implement me
		// Cannot do if means haven't been set
		if (means == null) {
			throw new Exception("Cannot compute local values of previous observations " +
					"without the means having been supplied (either directly or by " +
					"supplying observations)");
		}
		throw new RuntimeException("Not implemented yet");
	}

	/**
	 * Compute the local mutual information for a new series of 
	 *  observations, based on variances computed with the previously
	 *  supplied observations.
	 * 
	 * @param newSourceObs provided source observations
	 * @param newDestObs provided destination observations
	 * @return the local values in nats (not bits).
	 *  If the {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}
	 *  property was set to say k, then the local values align with the
	 *  destination value (i.e. after the given delay k). As such, the
	 *  first k values of the array will be zeros.  
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] newSourceObs,
			double[][] newDestObs) throws Exception {
		return computeLocalUsingPreviousObservations(newSourceObs, newDestObs, false);
	}

	/**
	 * Compute the local mutual information for a new series of 
	 *  observations, based on variances computed with the previously
	 *  supplied observations.
	 * 
	 * @param newSourceObs provided source observations
	 * @param newDestObs provided destination observations
	 * @param isPreviousObservations whether these are our previous
	 *  observations - this determines whether to add zeros for the first
	 *  timeDiff local values, and also 
	 *  whether to set the internal lastAverage field,
	 *  which is returned by later calls to {@link #getLastAverage()}
	 * @return the local values in nats (not bits).
	 *  If the {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}
	 *  property was set to say k, then the local values align with the
	 *  destination value (i.e. after the given delay k). As such, the
	 *  first k values of the array will be zeros.  
	 * @throws Exception
	 */
	protected double[] computeLocalUsingPreviousObservations(double[][] newSourceObs,
			double[][] newDestObs, boolean isPreviousObservations) throws Exception {
		
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
		
		// Pull out the source and dest components of the covariance,
		//  and take their inverses:
		int[] sourceIndicesInCovariance = MatrixUtils.range(0, dimensionsSource - 1);
		double[][] sourceCovariance = MatrixUtils.selectRowsAndColumns(covariance,
				sourceIndicesInCovariance, sourceIndicesInCovariance);
		double detSourceCovariance = MatrixUtils.determinant(sourceCovariance);
		double[][] invSourceCovariance = MatrixUtils.invertSymmPosDefMatrix(sourceCovariance);
		int[] destIndicesInCovariance = MatrixUtils.range(dimensionsSource,
			dimensionsSource + dimensionsDest - 1);
		double[][] destCovariance = MatrixUtils.selectRowsAndColumns(covariance,
				destIndicesInCovariance, destIndicesInCovariance);
		double detDestCovariance = MatrixUtils.determinant(destCovariance);
		double[][] invDestCovariance = MatrixUtils.invertSymmPosDefMatrix(destCovariance);
		
		double[] sourceMeans = MatrixUtils.select(means, 0, dimensionsSource);
		double[] destMeans = MatrixUtils.select(means, dimensionsSource, dimensionsDest);
		
		int lengthOfReturnArray, offset;
		if (isPreviousObservations && addedMoreThanOneObservationSet) {
			// We're returning the local values for a set of disjoint
			//  observations. So we don't add timeDiff zeros to the start,
			//  and note that the required timeDiff is already 
			//  built into the supplied observations.
			lengthOfReturnArray = newDestObs.length;
			offset = 0;
		} else {
			lengthOfReturnArray = newDestObs.length + timeDiff;
			offset = timeDiff;
		}
		// If we have a time delay, slide the local values
		double[] localValues = new double[lengthOfReturnArray];
		double newAverage = 0.0;
		for (int t = offset; t < newDestObs.length; t++) {
			// Computing local values for:
			//  a. sourceObservations[t - offset]
			//  b. destObservations[t]
			
			double[] sourceDeviationsFromMean =
					MatrixUtils.subtract(newSourceObs[t - offset],
							sourceMeans);
			double[] destDeviationsFromMean =
					MatrixUtils.subtract(newDestObs[t], destMeans);
			double[] deviationsFromMean =
					MatrixUtils.append(sourceDeviationsFromMean,
							destDeviationsFromMean);
			
			// Computing PDFs WITHOUT (2*pi)^dim factor, since these will cancel:
			double sourceExpArg = MatrixUtils.dotProduct(
						MatrixUtils.matrixProduct(sourceDeviationsFromMean,
								invSourceCovariance),
						sourceDeviationsFromMean);
			double adjustedPSource = Math.exp(-0.5 * sourceExpArg) /
						Math.sqrt(detSourceCovariance);
			double destExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(destDeviationsFromMean,
							invDestCovariance),
					destDeviationsFromMean);
			double adjustedPDest = Math.exp(-0.5 * destExpArg) /
					Math.sqrt(detDestCovariance);
			double jointExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(deviationsFromMean,
							invCovariance),
					deviationsFromMean);
			double adjustedPJoint = Math.exp(-0.5 * jointExpArg) /
					Math.sqrt(detCovariance);
			
			// Returning results in nats:
			double localValue = Math.log(adjustedPJoint /
					(adjustedPSource * adjustedPDest));
			localValues[t] = localValue;
			
			if (isPreviousObservations) {
				newAverage += localValue;
			}
		}
		
		if (isPreviousObservations) {
			lastAverage = newAverage / (newDestObs.length - timeDiff); 
		}
		
		return localValues;
	}

}
