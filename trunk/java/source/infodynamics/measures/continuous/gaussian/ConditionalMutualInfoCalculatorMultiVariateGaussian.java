package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.NonPositiveDefiniteMatrixException;

/**
 * <p>Computes the differential conditional mutual information of two given multivariate sets of
 *  observations,
 *  assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution.</p>
 *  
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct {@link #ConditionalMutualInfoCalculatorMultiVariateLinearGaussian()}</li>
 *		<li>{@link #initialise(int, int)}</li>
 * 		<li>Set properties using {@link #setProperty(String, String)}</li>
 * 		<li>Provide the observations to the calculator using:
 * 			{@link #setObservations(double[][], double[][])}, or
 * 			{@link #setCovariance(double[][])}, or
 * 			a sequence of:
 * 			{@link #startAddObservations()},
 *          multiple calls to either {@link #addObservations(double[][], double[][])}
 *          or {@link #addObservations(double[][], double[][], int, int)}, and then
 *          {@link #finaliseAddObservations()}.</li>
 * 		<li>Compute the required information-theoretic results, primarily:
 * 			{@link #computeAverageLocalOfObservations()} to return the average differential
 *          entropy based on either the set variance or the variance of
 *          the supplied observations; or other calls to compute
 *          local values or statistical significance.</li>
 * 	</ol>
 * </p>
 * 
 * <p>
 * Alters behaviour slightly from parent class {@link ConditionalMutualInfoMultiVariateCommon}
 * in that property {@link ConditionalMutualInfoMultiVariateCommon#PROP_NORMALISE}
 * is set to false by default here (since this makes more sense for
 * linear-Gaussian analysis). 
 * </p>
 * 
 * @see <a href="http://mathworld.wolfram.com/DifferentialEntropy.html">Differential entropy for Gaussian random variables at Mathworld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Differential_entropy">Differential entropy for Gaussian random variables at Wikipedia</a>
 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Multivariate normal distribution on Wikipedia</a>
 * @author Joseph Lizier joseph.lizier_at_gmail.com
 *
 */
public class ConditionalMutualInfoCalculatorMultiVariateGaussian 
		extends ConditionalMutualInfoMultiVariateCommon
		implements ConditionalMutualInfoCalculatorMultiVariate,
			AnalyticNullDistributionComputer, Cloneable {

	/**
	 * Cached Cholesky decomposition of the covariance matrix
	 * of the most recently supplied observations.
	 * Is a matrix [C_11, C_12, C_1c; C_21, C_22, C_2c; C_c1, C_c2, C_cc],
	 * where C_xy represents the covariance matrix of variable x to variable y
	 * where x,y are either variable 1, 2 or the conditional.
	 * The covariance matrix is symmetric, and should be positive definite
	 * (otherwise we have linealy dependent variables).
	 */
	protected double[][] L;
	/**
	 * Cached Cholesky decomposition of the (var1, conditional) covariance matrix
	 */
	protected double[][] L_1c;
	/**
	 * Cached Cholesky decomposition of the (var2, conditional) covariance matrix
	 */
	protected double[][] L_2c;
	/**
	 * Cached Cholesky decomposition of the conditional covariance matrix
	 */
	protected double[][] L_cc;
	
	/**
	 * Means of the most recently supplied observations (source variables
	 *  listed first, destination variables second).
	 */
	protected double[] means;

	/**
	 * Cached determinants of the covariance matrices
	 */
	protected double detCovariance;
	protected double det1cCovariance;
	protected double det2cCovariance;
	protected double detccCovariance;
	
	/**
	 * Cache which sub-variables for each variable are a linearly-independent set
	 *  and are used in the covariances
	 */
	protected int[] condIndicesInCovariance;
	protected int[] var1IndicesInCovariance;
	protected int[] var2IndicesInCovariance;
	
	public ConditionalMutualInfoCalculatorMultiVariateGaussian() {
		// Normalising data makes less sense for linear-Gaussian estimation,
		//  so we turn this off by default.
		normalise = false;
	}
	
	public void initialise(int var1Dimensions, int var2Dimensions, int condDimensions) {
		super.initialise(var1Dimensions, var2Dimensions, condDimensions);
		L = null;
		L_1c = null;
		L_2c = null;
		L_cc = null;
		means = null;
		detCovariance = 0;
		det1cCovariance = 0;
		det2cCovariance = 0;
		detccCovariance = 0;
		condIndicesInCovariance = null;
		var1IndicesInCovariance = null;
		var2IndicesInCovariance = null;
	}

	/**
	 * Finalise the addition of multiple observation sets.
	 * 
	 * @throws Exception if the observation variables are not linearly independent
	 *  (leading to a non-positive definite covariance matrix).
	 */
	public void finaliseAddObservations() throws Exception {

		// Get the observations properly stored in the sourceObservations[][] and
		//  destObservations[][] arrays.
		super.finaliseAddObservations();

		// Store the means of each variable (useful for local values later)
		means = new double[dimensionsVar1 + dimensionsVar2 + dimensionsCond];
		double[] var1Means = MatrixUtils.means(var1Observations);
		double[] var2Means = MatrixUtils.means(var2Observations);
		double[] condMeans = MatrixUtils.means(condObservations);
		System.arraycopy(var1Means, 0, means, 0, dimensionsVar1);
		System.arraycopy(var2Means, 0, means, dimensionsVar1, dimensionsVar2);
		System.arraycopy(condMeans, 0, means, dimensionsVar1 + dimensionsVar2,
				dimensionsCond);
		
		// Store the covariances of the variables
		// Generally, this should not throw an exception, since we checked
		//  the observations had the correct number of variables
		//  on receiving them, and in constructing the covariance matrix
		//   ourselves we know it should be symmetric.
		// It could occur however if the covariance matrix was not
		//  positive definite, which would occur if one variable
		//  is linearly redundant.
		setCovariance(
				MatrixUtils.covarianceMatrix(var1Observations, var2Observations, condObservations),
				true);
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  conditional mutual information.</p>
	 *  
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][])}.</p>
	 * 
	 * @param covariance covariance matrix of var1, var2, conditional
	 *  variables, considered together.
	 * @param numObservations the number of observations that the covariance
	 *  was determined from. This is used for later significance calculations
	 * @throws Exception for covariance matrix not matching the expected dimensions,
	 *  being non-square, asymmetric or non-positive definite
	 */
	public void setCovariance(double[][] covariance, int numObservations) throws Exception {
		setCovariance(covariance, false);
		totalObservations = numObservations;
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  conditional mutual information.</p>
	 *  
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][])}.</p>
	 * 
	 * @param covariance covariance matrix of var1, var2 and the conditional
	 *  variables, considered together.
	 * @param determinedFromObservations whether the covariance matrix
	 *  was determined internally from observations or not
	 * @throws Exception for covariance matrix not matching the expected dimensions,
	 *  being non-square, asymmetric or non-positive definite
	 */
	protected void setCovariance(double[][] covariance, boolean determinedFromObservations)
			throws Exception {
		if (!determinedFromObservations) {
			// Make sure we're not keeping any observations
			var1Observations = null;
			var2Observations = null;
			condObservations = null;
		}
		// Make sure the supplied covariance matrix matches the required dimenions:
		int rows = covariance.length;
		if (rows != dimensionsVar1 + dimensionsVar2 + dimensionsCond) {
			throw new Exception("Supplied covariance matrix does not match initialised number of dimensions");
		}
		
		// Now store the Cholesky decompositions for computing the cond MI later.
		// Start with the conditional variable:
		//  In case there are linear redundancies amongst the conditional variable,
		//  remove some of its components until the linear redundancy is gone.
		//  This allows a conditional MI computation to still take place.
		boolean redundanciesRemoved = false;
		condIndicesInCovariance = null;
		for (int v = dimensionsCond; v >= 1; v--) {
			// Select the first v variables in the conditional joint variable:
			condIndicesInCovariance = MatrixUtils.range(dimensionsVar1 + dimensionsVar2,
					dimensionsVar1 + dimensionsVar2 + v - 1);
			double[][] condCovariance =
				MatrixUtils.selectRowsAndColumns(covariance, 
						condIndicesInCovariance, condIndicesInCovariance);
			try {
				L_cc = MatrixUtils.CholeskyDecomposition(condCovariance);
			} catch (NonPositiveDefiniteMatrixException e) {
				// There is a linear redundancy between the variables - so allow
				//  one to be removed in the next loop iteration
				continue;
			}
			// Allow exceptions indicating asymmetric and non-square to be propagated
			// Otherwise, we've found a linearly independent subset of the conditioned variable
			redundanciesRemoved = true;
			break;
		}
		if (!redundanciesRemoved) {
			// This should not happen, because a univariate conditional should
			//  be positive definite
			throw new Exception("Conditional variables are linearly redundant, even " +
					"apparently when using only one of them (this should not happen)");
		}
		
		// And next store the Cholesky decompositions for var1 with
		//  the conditional variable:
		//  In case there are linear redundancies amongst the conditional variable
		//  with variable 1, remove some of the components of variable 1
		//  until the linear redundancy is gone.
		//  This allows a conditional MI computation to still take place.
		var1IndicesInCovariance = null;
		redundanciesRemoved = false;
		for (int v = dimensionsVar1; v >= 1; v--) {
			var1IndicesInCovariance = MatrixUtils.range(0, v - 1);
			int[] var1AndCondIndicesInCovariance = MatrixUtils.append(var1IndicesInCovariance, condIndicesInCovariance);
			double[][] var1AndCondCovariance =
					MatrixUtils.selectRowsAndColumns(covariance, 
							var1AndCondIndicesInCovariance, var1AndCondIndicesInCovariance);
			try {
				L_1c = MatrixUtils.CholeskyDecomposition(var1AndCondCovariance);
			} catch (NonPositiveDefiniteMatrixException e) {
				// There is a linear redundancy between the variables - so allow
				//  one to be removed in the next loop iteration
				continue;
			}
			// Allow exceptions indicating asymmetric and non-square to be propagated
			// Otherwise, we've found a linearly independent subset of the conditioned variable
			redundanciesRemoved = true;
			break;
		}
		if (!redundanciesRemoved) {
			// This means that variable 1 is *fully* linearly dependent on the conditioned
			//  variable (i.e. not even cutting it down to just one variable helped).
			// Flag this by setting:
			L_1c = null;
		}
		
		// And next store the Cholesky decompositions for var2 with
		//  the conditional variable:
		//  In case there are linear redundancies amongst the conditional variable
		//  with variable 2, remove some of the components of variable 2
		//  until the linear redundancy is gone.
		//  This allows a conditional MI computation to still take place.
		var2IndicesInCovariance = null;
		int[] var2AndCondIndicesInCovariance = null;
		redundanciesRemoved = false;
		for (int v = dimensionsVar2; v >= 1; v--) {
			var2IndicesInCovariance = MatrixUtils.range(dimensionsVar1, dimensionsVar1 + v - 1);
			var2AndCondIndicesInCovariance = MatrixUtils.append(var2IndicesInCovariance, condIndicesInCovariance);
			double[][] var2AndCondCovariance =
				MatrixUtils.selectRowsAndColumns(covariance, 
						var2AndCondIndicesInCovariance, var2AndCondIndicesInCovariance);
			try {
				L_2c = MatrixUtils.CholeskyDecomposition(var2AndCondCovariance);
			} catch (NonPositiveDefiniteMatrixException e) {
				// There is a linear redundancy between the variables - so allow
				//  one to be removed in the next loop iteration
				continue;
			}
			// Allow exceptions indicating asymmetric and non-square to be propagated
			// Otherwise, we've found a linearly independent subset of the conditioned variable
			redundanciesRemoved = true;
			break;
		}
		if (!redundanciesRemoved) {
			// This means that variable 2 is *fully* linearly dependent on the conditioned
			//  variable (i.e. not even cutting it down to just one variable helped).
			// Flag this by setting:
			L_2c = null;
		}

		// Finally, store the Cholesky decomposition for the whole covariance matrix:
		//  first prune the covariance in line with the variable removed
		//  from variable 1, 2 and the conditional:
		int[] prunedIndicesInCovariance = MatrixUtils.append(var1IndicesInCovariance, var2AndCondIndicesInCovariance);
		double[][] prunedCovariance =
				MatrixUtils.selectRowsAndColumns(covariance, 
						prunedIndicesInCovariance, prunedIndicesInCovariance);
		try {
			L = MatrixUtils.CholeskyDecomposition(prunedCovariance);
		} catch (NonPositiveDefiniteMatrixException e) {
			// There is a linear redundancy between the variables - 
			// Flag this by setting:
			L = null;
		}
		// Allow exceptions indicating asymmetric and non-square to be propagated

		// Postcondition: L's contain Cholesky decompositions of covariance
		//  matrices with linearly dependent variables removed (except for 
		//  the whole covariance matrix), using null to flag where this was not possible
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  mutual information.</p>
	 * 
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}.</p>
	 * 
	 * @param covariance covariance matrix of var1, var2 and conditional
	 *  variables, considered together.
	 * @param means mean of var1, var2 and conditional variables (as per
	 *  covariance)
	 * @param numObservations the number of observations that the mean and covariance
	 *  were determined from. This is used for later significance calculations
	 */
	public void setCovarianceAndMeans(double[][] covariance, double[] means,
			int numObservations) throws Exception {
		
		this.means = means;
		setCovariance(covariance, numObservations);
	}

	/**
	 * <p>The joint differential entropy for a multivariate Gaussian-distribution of dimension n
	 *  with covariance matrix C is -0.5*\log_e{(2*pi*e)^n*|det(C)|},
	 *  where det() is the matrix determinant of C.</p>
	 * 
	 * <p>Here we compute the conditional mutual information from the joint entropies
	 *  of all variables (H_12c), variable 1 and conditional (H_1c),
	 *  variable 2 and conditional (H_2c) and conditional (H_c),
	 *  giving MI = H_1c + H_2c - H_c - H_12c.
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
		// Simple way:
		// detCovariance = MatrixUtils.determinantSymmPosDefMatrix(covariance);
		// Using cached Cholesky decomposition:
		// And also with extended checks for linear redundancies:
		
		// Should always have a valid L_cc (since we can reduce it down to one variable:
		detccCovariance = MatrixUtils.determinantViaCholeskyResult(L_cc);
		if (L_1c == null) {
			// Variable 1 is fully linearly redundant with conditional, so
			//  we will have zero conditional MI:
			lastAverage = 0;
		} else {
			det1cCovariance = MatrixUtils.determinantViaCholeskyResult(L_1c);
			if (L_2c == null) {
				// Variable 2 is fully linearly redundant with conditional, so
				//  we will have zero conditional MI:
				lastAverage = 0;
			} else {
				det2cCovariance = MatrixUtils.determinantViaCholeskyResult(L_2c);
				if (L == null) {
					// There is a linear dependence amongst variables 1 and 2 given the
					//  conditional which did not exist for either with the conditional alone,
					//  so conditional MI diverges:
					lastAverage = Double.POSITIVE_INFINITY;
				} else {
					detCovariance = MatrixUtils.determinantViaCholeskyResult(L);

					// Else all the covariance matrices were ok, so compute as normal:
					lastAverage = 0.5 * Math.log(Math.abs(
							det1cCovariance * det2cCovariance /
									(detCovariance * detccCovariance)));
				}
			}
		}
		
		condMiComputed = true;
		return lastAverage;
	}

	/**
	 * <p>Compute the local or pointwise mutual information for each of the previously
	 * supplied observations</p>
	 * 
	 * @return array of the local values in nats (not bits!)
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		// Cannot do if destObservations haven't been set
		if (var2Observations == null) {
			throw new Exception("Cannot compute local values of previous observations " +
					"if they have not been set!");
		}
		
		return computeLocalUsingPreviousObservations(var1Observations,
				var2Observations, condObservations, true);
	}

	/**
	 * <p>Compute the statistical significance of the conditional mutual information 
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
	public ChiSquareMeasurementDistribution computeSignificance() throws Exception {
		if (!condMiComputed) {
			computeAverageLocalOfObservations();
		}
		// Number of extra parameters in the model incorporating the
		//  extra variable is independent of the number of variables
		//  in the conditional:
		// (Assuming that all variables went into the calculation:)
		// return new ChiSquareMeasurementDistribution(2.0*((double)totalObservations)*lastAverage,
		//		dimensionsVar1 * dimensionsVar2);
		// Taking the subsets into account:
		return new ChiSquareMeasurementDistribution(2.0*((double)totalObservations)*lastAverage,
				var1IndicesInCovariance.length * var2IndicesInCovariance.length);
	}
	
	
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon#computeSignificance(int, int)
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int numPermutationsToCheck) throws Exception {
		if (var2Observations == null) {
			throw new Exception("Cannot compute empirical statistical significance " +
					"if user passed in covariance matrix rather than observations.");
		}

		return super.computeSignificance(variableToReorder, numPermutationsToCheck);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon#computeSignificance(int, int[][])
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int[][] newOrderings) throws Exception {
		if (var2Observations == null) {
			throw new Exception("Cannot compute empirical statistical significance " +
					"if user passed in covariance matrix rather than observations.");
		}

		return super.computeSignificance(variableToReorder, newOrderings);
	}

	/**
	 * @return the number of previously supplied observations for which
	 *  the conditional mutual information will be / was computed.
	 */
	public int getNumObservations() throws Exception {
		if (var2Observations == null) {
			throw new Exception("Cannot return number of observations because either " +
					"this calculator has not had observations supplied or " +
					"the user supplied the covariance matrix instead of observations");
		}
		return super.getNumObservations();
	}

	/**
	 * Compute the conditional mutual information if the given variable was
	 *  ordered as per the ordering specified in newOrdering
	 * 
	 * @param newOrdering array of time indices with which to reorder the data
	 * @return a surrogate conditional MI evaluated for the given ordering of the source variable
	 * @throws Exception if the user previously supplied covariance directly rather
	 *  than by setting observations (this means we have no observations
	 *  to reorder).
	 */
	public double computeAverageLocalOfObservations(int variableToReorder, 
			int[] newOrdering) throws Exception {
		// Cannot do if observations haven't been set (i.e. the variances
		//  were directly supplied)
		if (var1Observations == null) {
			throw new Exception("Cannot compute local values of previous observations " +
					"without supplying observations");
		}
		return super.computeAverageLocalOfObservations(variableToReorder, newOrdering);
	}

	/**
	 * Compute the local conditional mutual information for a new series of 
	 *  observations, based on variances computed with the previously
	 *  supplied observations.
	 * 
	 * @param newVar1Obs provided variable 1 observations
	 * @param newVar2Obs provided variable 2 observations
	 * @param newCondObs provided conditional observations
	 * @return the local values in nats (not bits).
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] newVar1Obs,
			double[][] newVar2Obs, double[][] newCondObs) throws Exception {
		return computeLocalUsingPreviousObservations(
				newVar1Obs, newVar2Obs, newCondObs, false);
	}

	/**
	 * Compute the local conditional mutual information for a new series of 
	 *  observations, based on variances computed with the previously
	 *  supplied observations.
	 * 
	 * @param newVar1Obs provided variable 1 observations
	 * @param newVar2Obs provided variable 2 observations
	 * @param newCondObs provided conditional observations
	 * @param isPreviousObservations whether these are our previous
	 *  observations - this determines whether to 
	 *  set the internal lastAverage field,
	 *  which is returned by later calls to {@link #getLastAverage()}
	 * @return the local values in nats (not bits).
	 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Multivariate normal distribution on Wikipedia</a>
	 * @see <a href="http://en.wikipedia.org/wiki/Positive-definite_matrix>"Positive definite matrix in Wikipedia"</a>
	 * @throws Exception if means were not defined by {@link #setObservations(double[][], double[][])} etc
	 *  or {@link #setCovarianceAndMeans(double[][], double[])}
	 * TODO Need to rewrite this to handle the subset of linearly independent variables only
	 */
	protected double[] computeLocalUsingPreviousObservations(double[][] newVar1Obs,
			double[][] newVar2Obs, double[][] newCondObs, boolean isPreviousObservations) throws Exception {
		
		if (means == null) {
			throw new Exception("Cannot compute local values without having means either supplied or computed via setObservations()");
		}

		// Check that the covariance matrix was positive definite:
		// (this was done earlier in computing the Cholesky decomposition,
		//  we may still need to compute the determinant)
		if (detCovariance == 0) {
			// The determinant has not been computed yet
			// Simple way:
			// detCovariance = MatrixUtils.determinantSymmPosDefMatrix(covariance);
			// Using cached Cholesky decomposition:

			// Should always have a valid L_cc (since we can reduce it down to one variable:
			detccCovariance = MatrixUtils.determinantViaCholeskyResult(L_cc);
			if (L_1c == null) {
				// Variable 1 is fully linearly redundant with conditional, so
				//  we will have zero conditional MI:
				return MatrixUtils.constantArray(newVar2Obs.length, 0);
			} else {
				det1cCovariance = MatrixUtils.determinantViaCholeskyResult(L_1c);
				if (L_2c == null) {
					// Variable 2 is fully linearly redundant with conditional, so
					//  we will have zero conditional MI:
					return MatrixUtils.constantArray(newVar2Obs.length, 0);
				} else {
					det2cCovariance = MatrixUtils.determinantViaCholeskyResult(L_2c);
					if (L == null) {
						// There is a linear dependence amongst variables 1 and 2 given the
						//  conditional which did not exist for either with the conditional alone,
						//  so conditional MI diverges:
						return MatrixUtils.constantArray(newVar2Obs.length, Double.POSITIVE_INFINITY);
					} else {
						detCovariance = MatrixUtils.determinantViaCholeskyResult(L);
					}
				}
			}
		}

		// Now we are clear to take the matrix inverse (via Cholesky decomposition,
		//  since we have a symmetric positive definite matrix):
		double[][] invCovariance = MatrixUtils.solveViaCholeskyResult(L,
				MatrixUtils.identityMatrix(L.length));
		double[][] invVar1CondCovariance = MatrixUtils.solveViaCholeskyResult(L_1c,
				MatrixUtils.identityMatrix(L_1c.length));
		double[][] invVar2CondCovariance = MatrixUtils.solveViaCholeskyResult(L_2c,
				MatrixUtils.identityMatrix(L_2c.length));
		double[][] invCondCovariance = MatrixUtils.solveViaCholeskyResult(L_cc,
				MatrixUtils.identityMatrix(L_cc.length));
			
		// Now, only use the means from the subsets of linearly independent variables:
		// double[] var1Means = MatrixUtils.select(means, 0, dimensionsVar1);
		double[] var1Means = MatrixUtils.select(means, var1IndicesInCovariance);
		// double[] var2Means = MatrixUtils.select(means, dimensionsVar1, dimensionsVar2);
		double[] var2Means = MatrixUtils.select(means, var2IndicesInCovariance);
		// double[] condMeans = MatrixUtils.select(means, dimensionsVar1 + dimensionsVar2, dimensionsCond);
		double[] condMeans = MatrixUtils.select(means, condIndicesInCovariance);
		
		int lengthOfReturnArray;
		lengthOfReturnArray = newVar2Obs.length;

		double[] localValues = new double[lengthOfReturnArray];
		int[] var2IndicesSelected = MatrixUtils.subtract(var2IndicesInCovariance, dimensionsVar1);
		int[] condIndicesSelected = MatrixUtils.subtract(condIndicesInCovariance, dimensionsVar1 + dimensionsVar2);
		for (int t = 0; t < newVar2Obs.length; t++) {
			
			double[] var1DeviationsFromMean =
					MatrixUtils.subtract(
							MatrixUtils.select(newVar1Obs[t], var1IndicesInCovariance),
							var1Means);
			double[] var2DeviationsFromMean =
					MatrixUtils.subtract(
							MatrixUtils.select(newVar2Obs[t], var2IndicesSelected),
							var2Means);
			double[] condDeviationsFromMean =
					MatrixUtils.subtract(
							MatrixUtils.select(newCondObs[t], condIndicesSelected),
							condMeans);
			double[] var1CondDeviationsFromMean =
					MatrixUtils.append(var1DeviationsFromMean,
							condDeviationsFromMean);
			double[] var2CondDeviationsFromMean =
					MatrixUtils.append(var2DeviationsFromMean,
							condDeviationsFromMean);
			double[] tempDeviationsFromMean =
					MatrixUtils.append(var1DeviationsFromMean,
							var2DeviationsFromMean);
			double[] deviationsFromMean =
					MatrixUtils.append(tempDeviationsFromMean,
							condDeviationsFromMean);
			
			// Computing PDFs WITHOUT (2*pi)^dim factor, since these will cancel:
			// (see the PDFs defined at the wikipedia page referenced in the method header)
			double var1CondExpArg = MatrixUtils.dotProduct(
						MatrixUtils.matrixProduct(var1CondDeviationsFromMean,
								invVar1CondCovariance),
						var1CondDeviationsFromMean);
			double adjustedPVar1Cond = Math.exp(-0.5 * var1CondExpArg) /
						Math.sqrt(det1cCovariance);
			double var2CondExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(var2CondDeviationsFromMean,
							invVar2CondCovariance),
					var2CondDeviationsFromMean);
			double adjustedPVar2Cond = Math.exp(-0.5 * var2CondExpArg) /
					Math.sqrt(det2cCovariance);
			double condExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(condDeviationsFromMean,
							invCondCovariance),
					condDeviationsFromMean);
			double adjustedPCond = Math.exp(-0.5 * condExpArg) /
					Math.sqrt(detccCovariance);
			double jointExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(deviationsFromMean,
							invCovariance),
					deviationsFromMean);
			double adjustedPJoint = Math.exp(-0.5 * jointExpArg) /
					Math.sqrt(detCovariance);
			
			// Returning results in nats:
			double localValue = Math.log(adjustedPJoint * adjustedPCond /
					(adjustedPVar1Cond * adjustedPVar2Cond));
			localValues[t] = localValue;
		}
		
		// if (isPreviousObservations) {
		// Don't store the average value here, since it won't be exactly 
		//  the same as what would have been computed under the analytic expression
		// }
		
		return localValues;
	}
}
