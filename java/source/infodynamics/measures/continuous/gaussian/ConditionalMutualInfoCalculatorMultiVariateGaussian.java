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

package infodynamics.measures.continuous.gaussian;

import java.util.ArrayList;

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.NonPositiveDefiniteMatrixException;

/**
 * <p>Computes the differential conditional mutual information of two multivariate
 *  <code>double[][]</code> sets of observations, conditioned on another
 *  (implementing {@link ConditionalMutualInfoCalculatorMultiVariate}),
 *  assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution.
 *  This is achieved by extending the common code base in {@link ConditionalMutualInfoMultiVariateCommon}.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link ConditionalMutualInfoCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #ConditionalMutualInfoCalculatorMultiVariateGaussian()}.</li>
 *  <li>The property {@link ConditionalMutualInfoMultiVariateCommon#PROP_NORMALISE}
 *     is set to false by default here (since this makes more sense for
 *     linear-Gaussian analysis), which is different to the parent class.</li>
 * 	<li>The user can call {@link #setCovariance(double[][], boolean)} or
 *     {@link #setCovariance(double[][], int)} or {@link #setCovarianceAndMeans(double[][], double[], int)}
 *     instead of supplying observations via {@link #setObservations(double[][], double[][], double[][])} or
 *     {@link #addObservations(double[][], double[][], double[][])} etc.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  <li>Additional method {@link #computeSignificance()} to compute null distribution analytically.</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991).</li>
 * </ul>
 * 
 * @see <a href="http://mathworld.wolfram.com/DifferentialEntropy.html">Differential entropy for Gaussian random variables at Mathworld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Differential_entropy">Differential entropy for Gaussian random variables at Wikipedia</a>
 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Multivariate normal distribution on Wikipedia</a>
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class ConditionalMutualInfoCalculatorMultiVariateGaussian 
		extends ConditionalMutualInfoMultiVariateCommon
		implements ConditionalMutualInfoCalculatorMultiVariate,
			AnalyticNullDistributionComputer, Cloneable {

	/**
	 * Property name for whether analytically-determined bias
	 *  is to be corrected out of estimated provided by the calculator.
	 */
	public static final String PROP_BIAS_CORRECTION = "BIAS_CORRECTION";
	/**
	 * Whether to analytically bias correct the returned values 
	 */
	protected boolean biasCorrection = false;
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
	 * Cached Cholesky decomposition of the (conditional, var1) covariance matrix
	 */
	protected double[][] L_c1;
	/**
	 * Cached Cholesky decomposition of the (conditional, var2) covariance matrix
	 */
	protected double[][] L_c2;
	/**
	 * Cached Cholesky decomposition of the conditional covariance matrix
	 */
	protected double[][] L_cc;
	
	/**
	 * Cached determinants of the joint covariance matrix
	 */
	protected double detCovariance;
	/**
	 * Cached determinants of the covariance matrix for
	 *  variable 1 and conditional
	 */
	protected double detc1Covariance;
	/**
	 * Cached determinants of the covariance matrix for
	 *  variable 2 and the conditional
	 */
	protected double detc2Covariance;
	/**
	 * Cached determinants of the covariance matrix
	 *  for the conditional
	 */
	protected double detccCovariance;
	
	/**
	 * Cache the sub-variables which are a linearly-independent set
	 *  (and so are used in the covariances) for variable 1 (given conditionals)
	 */
	protected int[] var1IndicesInCovariance;
	/**
	 * Cache the sub-variables which are a linearly-independent set
	 *  (and so are used in the covariances) for variable 2 (given conditionals)
	 */
	protected int[] var2IndicesInCovariance;
	/**
	 * Cache the sub-variables which are a linearly-independent set
	 *  (and so are used in the covariances) for the conditional
	 */
	protected int[] condIndicesInCovariance;

	/**
	 * Construct an instance
	 */
	public ConditionalMutualInfoCalculatorMultiVariateGaussian() {
		// Normalising data makes less sense for linear-Gaussian estimation,
		//  so we turn this off by default.
		normalise = false;
	}
	
	@Override
	public void initialise(int var1Dimensions, int var2Dimensions, int condDimensions) {
		super.initialise(var1Dimensions, var2Dimensions, condDimensions);
		L = null;
		L_c1 = null;
		L_c2 = null;
		L_cc = null;
		detCovariance = 0;
		detc1Covariance = 0;
		detc2Covariance = 0;
		detccCovariance = 0;
		condIndicesInCovariance = null;
		var1IndicesInCovariance = null;
		var2IndicesInCovariance = null;
	}

	/**
	 * Sets properties for the Gaussian CMI calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 *  <li>{@link #PROP_BIAS_CORRECTION} -- if set to "true", then the analytically determined bias
	 *      (as the mean of the surrogate distribution) will be subtracted from all
	 *      calculated values.
	 *      Default is "false".
	 *  <li>any valid properties for {@link ConditionalMutualInfoMultiVariateCommon#setProperty(String, String)}</li>
	 * </ul>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_BIAS_CORRECTION)) {
			biasCorrection = Boolean.parseBoolean(propertyValue);
		} else {
			// No property was set here
			propertySet = false;
			// try the superclass:
			super.setProperty(propertyName, propertyValue);
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/**
	 * Get property values for the calculator.
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, are the same as those for
	 * {@link #setProperty(String, String)}</p>
	 * 
	 * <p>Unknown property values are responded to with a null return value.</p>
	 * 
	 * @param propertyName name of the property
	 * @return current value of the property
	 * @throws Exception for invalid property values
	 */
	public String getProperty(String propertyName) {
	
		if (propertyName.equalsIgnoreCase(PROP_BIAS_CORRECTION)) {
			return Boolean.toString(biasCorrection);
		} else {
			// try the superclass:
			return super.getProperty(propertyName);
		}
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

	/**
	 * @throws Exception if the observation variables are not linearly independent
	 *  (leading to a non-positive definite covariance matrix).
	 */
	@Override
	public void finaliseAddObservations() throws Exception {

		// Get the observations properly stored in the sourceObservations[][] and
		//  destObservations[][] arrays.
		super.finaliseAddObservations();

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
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}.</p>
	 * 
	 * @param covariance covariance matrix of var1, var2, conditional
	 *  variables, considered jointly together.
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
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}.</p>
	 * 
	 * @param covariance covariance matrix of var1, var2 and the conditional
	 *  variables, considered jointly together.
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
		// Make sure the supplied covariance matrix matches the required dimensions:
		if (covariance.length != dimensionsVar1 + dimensionsVar2 + dimensionsCond) {
			throw new Exception("Supplied covariance matrix does not match initialised number of dimensions");
		}
		// Check also for non-square matrix: the later calls for Cholesky decompositions
		//  will only check the sub-matrices supplied to them
		for (int r = 0; r < covariance.length; r++) {
			if (covariance[r].length != dimensionsVar1 + dimensionsVar2 + dimensionsCond) {
				throw new Exception("Number of columns for row " + r +
						" of supplied covariance matrix does not match initialised number of dimensions");
			}
		}

		// Now store the Cholesky decompositions for computing the cond MI later.
		// Start with the conditional variable:
		ArrayList<Integer> condIndicesSet = MatrixUtils.createArrayList(
				MatrixUtils.range(dimensionsVar1 + dimensionsVar2,
						dimensionsVar1 + dimensionsVar2 + dimensionsCond - 1));
		L_cc = MatrixUtils.makeCholeskyOfIndependentComponents(covariance, condIndicesSet, null);
		// Postcondition: condIndicesSet stores a set of independent components of the conditionals
		condIndicesInCovariance = MatrixUtils.toArray(condIndicesSet);

		// Next var1, given conditionals:
		ArrayList<Integer> var1IndicesSet = MatrixUtils.createArrayList(
				MatrixUtils.range(0, dimensionsVar1 - 1));
		L_c1 = MatrixUtils.makeCholeskyOfIndependentComponents(covariance, var1IndicesSet, condIndicesSet);
		// Postcondition: var1IndicesSet stores a set of independent components of var1, given cond.
		// L_c1 is computed for the variables in order condIndicesSet then var1IndicesSet
		var1IndicesInCovariance = MatrixUtils.toArray(var1IndicesSet);
		
		// Next var2, given conditionals:
		ArrayList<Integer> var2IndicesSet = MatrixUtils.createArrayList(
				MatrixUtils.range(dimensionsVar1, dimensionsVar1 + dimensionsVar2 - 1));
		L_c2 = MatrixUtils.makeCholeskyOfIndependentComponents(covariance, var2IndicesSet, condIndicesSet);
		// Postcondition: var2IndicesSet stores a set of independent components of var2, given cond.
		// L_c2 is computed for the variables in order condIndicesSet then var2IndicesSet
		var2IndicesInCovariance = MatrixUtils.toArray(var2IndicesSet);
		
		// Finally, store the Cholesky decomposition for the whole covariance matrix.
		// The order var1-var2-conditional is important for the local evaluations later.
		ArrayList<Integer> varsForCovarianceSet = new ArrayList<Integer>(var1IndicesSet);
		varsForCovarianceSet.addAll(var2IndicesSet);
		varsForCovarianceSet.addAll(condIndicesSet);
		double[][] prunedCovariance =
				MatrixUtils.selectRowsAndColumns(covariance, 
						varsForCovarianceSet, varsForCovarianceSet);
		try {
			L = MatrixUtils.CholeskyDecomposition(prunedCovariance);
		} catch (NonPositiveDefiniteMatrixException e) {
			// There is a linear redundancy between the var1 and var2 given the conditional
			//  (it's not possible that it was for var1 or var 2 with conditional, given 
			//  that we've already pruned these above, so we need not try to remove the redundancy.
			// Flag this by setting:
			L = null;
		}
		// Allow exceptions indicating asymmetric and non-square to be propagated
		
		// Postcondition: L's contain Cholesky decompositions of covariance
		//  matrices with linearly dependent variables removed,
		//  using null to flag where this was not possible
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  mutual information, as well as the means for each variable.</p>
	 * 
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}
	 *  (but you can call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}).</p>
	 * 
	 * @param covariance covariance matrix of var1, var2 and conditional
	 *  variables, considered jointly together.
	 * @param means mean of var1, var2 and conditional variables (as per
	 *  <code>covariance</code>)
	 * @param numObservations the number of observations that the mean and covariance
	 *  were determined from. This is used for later significance calculations
	 */
	public void setCovarianceAndMeans(double[][] covariance, double[] means,
			int numObservations) throws Exception {
		
		var1Means = MatrixUtils.select(means, 0, dimensionsVar1);
		var2Means = MatrixUtils.select(means, dimensionsVar1, dimensionsVar2);
		if (dimensionsCond > 0 ) {
			condMeans = MatrixUtils.select(means, dimensionsVar1 + dimensionsVar2, dimensionsCond);
		}
		
		setCovariance(covariance, numObservations);

		// set the std deviations from the covariance:
		var1Stds = new double[dimensionsVar1];
		var2Stds = new double[dimensionsVar2];
		condStds = new double[dimensionsCond];
		for (int i = 0; i < covariance.length; i++) {
			if (i < dimensionsVar1) {
				var1Stds[i] = Math.sqrt(covariance[i][i]);
			} else if (i < dimensionsVar1 + dimensionsVar2) {
				var2Stds[i - dimensionsVar1] = Math.sqrt(covariance[i][i]);
			} else {
				condStds[i - dimensionsVar1 - dimensionsVar2] = Math.sqrt(covariance[i][i]);
			}
		
		}
	}

	/**
	 * <p>Computes the local values of the conditional mutual information,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>The joint differential entropy for a multivariate Gaussian-distribution of dimension n
	 *  with covariance matrix C is -0.5*\log_e{(2*pi*e)^n*|det(C)|},
	 *  where det() is the matrix determinant of C.</p>
	 * 
	 * <p>Here we compute the conditional mutual information from the joint entropies
	 *  of all variables (H_12c), conditional and variable 1 (H_c1),
	 *  conditional and variable 2 (H_c2) and conditional (H_c),
	 *  giving MI = H_c1 + H_c2 - H_c - H_12c.
	 *  We assume that the recorded estimation of the
	 *  covariance is correct (i.e. we will not make a bias correction for limited
	 *  observations here).</p>
	 * 
	 * @return the conditional mutual information of the previously provided observations or from the
	 *  supplied covariance matrix, in <b>nats</b> (not bits!).
	 *  Returns NaN if any of the determinants are zero
	 *  (because this will make the denominator of the log zero)
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		// Simple way:
		// detCovariance = MatrixUtils.determinantSymmPosDefMatrix(covariance);
		// Using cached Cholesky decomposition:
		// And also with extended checks for linear redundancies:
		
		if (L_c1 == null) {
			// Variable 1 is fully linearly redundant with conditional, so
			//  we will have zero conditional MI:
			lastAverage = 0;
		} else {
			detc1Covariance = MatrixUtils.determinantViaCholeskyResult(L_c1);
			if (L_c2 == null) {
				// Variable 2 is fully linearly redundant with conditional, so
				//  we will have zero conditional MI:
				lastAverage = 0;
			} else {
				detc2Covariance = MatrixUtils.determinantViaCholeskyResult(L_c2);
				if (L == null) {
					// There is a linear dependence amongst variables 1 and 2 given the
					//  conditional which did not exist for either with the conditional alone,
					//  so conditional MI diverges:
					lastAverage = Double.POSITIVE_INFINITY;
				} else {
					detCovariance = MatrixUtils.determinantViaCholeskyResult(L);

					// Else all the covariance matrices were ok, except perhaps the
					//  conditional covariance
					if (L_cc == null) {
						// The conditional variables had no covariance, so
						//  just return an MI:
						lastAverage = 0.5 * Math.log(Math.abs(
								detc1Covariance * detc2Covariance /
										detCovariance));
					} else {
						detccCovariance = MatrixUtils.determinantViaCholeskyResult(L_cc);
						// So compute as normal:
						lastAverage = 0.5 * Math.log(Math.abs(
									detc1Covariance * detc2Covariance /
											(detCovariance * detccCovariance)));
					}
				}
			}
		}
		
		if (biasCorrection) {
			// Need to play a slight trick here so that computeSignificance()
			//  thinks the average has already been computed, otherwise
			//  we will get an infinite loop where it calls this method again, and so on:
			
			ChiSquareMeasurementDistribution analyticMeasDist = computeSignificance(true);
			lastAverage -= analyticMeasDist.getMeanOfUncorrectedDistribution();
		}
		condMiComputed = true;
		return lastAverage;
	}

	/**
	 * @return array of the local values in nats (not bits!)
	 * @throws Exception if the user had previously supplied covariances
	 *  directly (ie had not supplied observations).
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
	 * Generate an <b>analytic</b> distribution of what the conditional MI would look like,
	 * under a null hypothesis that our variables had no relation
	 * (in the context of the conditional value).
	 * This is performed without bootstrapping (which is done in
	 * {@link #computeSignificance(int, int)} and {@link #computeSignificance(int, int[][])}).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below, and the other papers referenced in
	 * {@link AnalyticNullDistributionComputer#computeSignificance()}
	 * (in particular Geweke),
	 * for a description of how this is done for conditional MI.
	 * Basically, the null distribution is a chi-square distribution 
	 * with degrees of freedom equal to the product of the number of variables
	 * in each joint variable 1 and 2.
	 * </p>
	 * 
	 * @return ChiSquareMeasurementDistribution object which describes
	 * the proportion of conditional MI scores from the null distribution
	 *  which have higher or equal conditional MIs to our actual value.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception
	 */
	@Override
	public ChiSquareMeasurementDistribution computeSignificance() throws Exception {
		return computeSignificance(false);
	}

	/**
	 * As per {@link #computeSignificance()} except allows the caller
	 *  to request that the average is not first computed (if we don't have
	 *  it already). This is required internally to avoid infinite looping
	 *  between computeAverage and computeSignificance, when we just want 
	 *  the null distribution and don't need to have the pValue
	 * 
	 * @param skipComputingThisAverage
	 * @return
	 * @throws Exception
	 */
	protected ChiSquareMeasurementDistribution computeSignificance(boolean skipComputingThisAverage) throws Exception {
		double averageToUse = lastAverage;
		if (!condMiComputed) {
			if (skipComputingThisAverage) {
				averageToUse = 0; // Caller only wants the distribution
			} else {
				averageToUse = computeAverageLocalOfObservations();
			}
		}
		// else use 0 for now in the distribution
		
		// Number of extra parameters in the model incorporating the
		//  extra variable is independent of the number of variables
		//  in the conditional:
		// (Assuming that all variables went into the calculation:)
		// return new ChiSquareMeasurementDistribution(2.0*((double)totalObservations)*lastAverage,
		//		dimensionsVar1 * dimensionsVar2);
		// Taking the subsets into account:
		return new ChiSquareMeasurementDistribution(averageToUse,
				totalObservations,
				var1IndicesInCovariance.length * var2IndicesInCovariance.length,
				biasCorrection);
	}

	/**
	 * @throws Exception if user passed in covariance matrix rather than observations
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

	/**
	 * @throws Exception if user passed in covariance matrix rather than observations
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

	public int getNumObservations() throws Exception {
		if (var2Observations == null) {
			throw new Exception("Cannot return number of observations because either " +
					"this calculator has not had observations supplied or " +
					"the user supplied the covariance matrix instead of observations");
		}
		return super.getNumObservations();
	}

	/**
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
	 * @return the series of local conditional MI values in <b>nats</b> (not bits).
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] states1,
			double[][] states2, double[][] condStates) throws Exception {
		return computeLocalUsingPreviousObservations(
				states1, states2, condStates, false);
	}

	/**
	 * Utility function to implement either {@link #computeLocalOfPreviousObservations()}
	 * or {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}
	 * depending on value of <code>isPreviousObservations</code>
	 * 
	 * @param newVar1Obs provided variable 1 observations
	 * @param newVar2Obs provided variable 2 observations
	 * @param newCondObs provided conditional observations
	 * @param isPreviousObservations whether these are our previous
	 *  observations - this determines whether we are implementing
	 *  {@link #computeLocalOfPreviousObservations()} if true 
	 *  or {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])} 
	 *  if false. Also indicates whether to 
	 *  set the internal lastAverage field,
	 *  which is returned by later calls to {@link #getLastAverage()} 
	 *  (in case of true).
	 * @return the local conditional MI values in nats (not bits).
	 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Multivariate normal distribution on Wikipedia</a>
	 * @see <a href="http://en.wikipedia.org/wiki/Positive-definite_matrix>"Positive definite matrix in Wikipedia"</a>
	 * @throws Exception if means were not defined by {@link #setObservations(double[][], double[][])} etc
	 *  or {@link #setCovarianceAndMeans(double[][], double[])}
	 */
	protected double[] computeLocalUsingPreviousObservations(double[][] newVar1Obs,
			double[][] newVar2Obs, double[][] newCondObs, boolean isPreviousObservations) throws Exception {
		
		if (var1Means == null) {
			throw new Exception("Cannot compute local values without having means either supplied or computed via setObservations()");
		}


		if ((!isPreviousObservations) && normalise) {
			// Need to normalise new observations
			newVar1Obs = MatrixUtils.normaliseIntoNewArray(newVar1Obs, var1Means, var1Stds);
			newVar2Obs = MatrixUtils.normaliseIntoNewArray(newVar2Obs, var2Means, var2Stds);
			if (dimensionsCond > 0) {
				newCondObs = MatrixUtils.normaliseIntoNewArray(newCondObs, condMeans, condStds);
			}
		}
		
		// And in case we need this for bias correction:
		ChiSquareMeasurementDistribution analyticMeasDist = computeSignificance(true);

		// Check that the covariance matrix was positive definite:
		// (this was done earlier in computing the Cholesky decomposition,
		//  we may still need to compute the determinant)
		if (detCovariance == 0) {
			// The determinant has not been computed yet
			// Simple way:
			// detCovariance = MatrixUtils.determinantSymmPosDefMatrix(covariance);
			// Using cached Cholesky decomposition:

			if (L_cc == null) {
				// We will compute local MIs
				detccCovariance = 0;
			} else {
				detccCovariance = MatrixUtils.determinantViaCholeskyResult(L_cc);
			}
			if (L_c1 == null) {
				// Variable 1 is fully linearly redundant with conditional, so
				//  we will have zero conditional MI:
				return MatrixUtils.constantArray(newVar2Obs.length,
						biasCorrection ? -analyticMeasDist.getMeanOfUncorrectedDistribution() : 0);
			} else {
				detc1Covariance = MatrixUtils.determinantViaCholeskyResult(L_c1);
				if (L_c2 == null) {
					// Variable 2 is fully linearly redundant with conditional, so
					//  we will have zero conditional MI:
					return MatrixUtils.constantArray(newVar2Obs.length, 
							biasCorrection ? -analyticMeasDist.getMeanOfUncorrectedDistribution() : 0);
				} else {
					detc2Covariance = MatrixUtils.determinantViaCholeskyResult(L_c2);
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
		double[][] invCondVar1Covariance = MatrixUtils.solveViaCholeskyResult(L_c1,
				MatrixUtils.identityMatrix(L_c1.length));
		double[][] invCondVar2Covariance = MatrixUtils.solveViaCholeskyResult(L_c2,
				MatrixUtils.identityMatrix(L_c2.length));
		double[][] invCondCovariance = null;
		if (L_cc != null) {
			invCondCovariance = MatrixUtils.solveViaCholeskyResult(L_cc,
					MatrixUtils.identityMatrix(L_cc.length));
		}
		
		// Use the following arrays to index directly into dimensions of the var2 and cond sample vectors.
		// We don't need a var1IndicesSelected because there is no offset
		//  from zero for var1IndicesInCovariance (unlike for the others)
		int[] var2IndicesSelected = MatrixUtils.subtract(var2IndicesInCovariance, dimensionsVar1);
		int[] condIndicesSelected = MatrixUtils.subtract(condIndicesInCovariance, dimensionsVar1 + dimensionsVar2);
		
		// Now, only use the means from the subsets of linearly independent variables:
		double[] var1SelectedMeans = MatrixUtils.select(var1Means, var1IndicesInCovariance);
		double[] var2SelectedMeans = MatrixUtils.select(var2Means, var2IndicesSelected);
		double[] condSelectedMeans = MatrixUtils.select(condMeans, condIndicesSelected);
		
		int lengthOfReturnArray;
		lengthOfReturnArray = newVar2Obs.length;

		double[] localValues = new double[lengthOfReturnArray];

		if (isPreviousObservations) {
			lastAverage = 0;
		}
		
		for (int t = 0; t < newVar2Obs.length; t++) {
			
			double[] var1DeviationsFromMean =
					MatrixUtils.subtract(
							MatrixUtils.select(newVar1Obs[t], var1IndicesInCovariance),
							var1SelectedMeans);
			double[] var2DeviationsFromMean =
					MatrixUtils.subtract(
							MatrixUtils.select(newVar2Obs[t], var2IndicesSelected),
							var2SelectedMeans);
			double[] condDeviationsFromMean =
					MatrixUtils.subtract(
							MatrixUtils.select(newCondObs[t], condIndicesSelected),
							condSelectedMeans);
			double[] condVar1DeviationsFromMean =
					MatrixUtils.append(condDeviationsFromMean,
							var1DeviationsFromMean);
			double[] condVar2DeviationsFromMean =
					MatrixUtils.append(condDeviationsFromMean,
							var2DeviationsFromMean);
			double[] tempDeviationsFromMean =
					MatrixUtils.append(var1DeviationsFromMean,
							var2DeviationsFromMean);
			double[] deviationsFromMean =
					MatrixUtils.append(tempDeviationsFromMean,
							condDeviationsFromMean);
			
			// Computing PDFs WITHOUT (2*pi)^dim factor, since these will cancel:
			// (see the PDFs defined at the wikipedia page referenced in the method header)
			double condVar1ExpArg = MatrixUtils.dotProduct(
						MatrixUtils.matrixProduct(condVar1DeviationsFromMean,
								invCondVar1Covariance),
						condVar1DeviationsFromMean);
			double adjustedPCondVar1 = Math.exp(-0.5 * condVar1ExpArg) /
						Math.sqrt(detc1Covariance);
			double condVar2ExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(condVar2DeviationsFromMean,
							invCondVar2Covariance),
					condVar2DeviationsFromMean);
			double adjustedPCondVar2 = Math.exp(-0.5 * condVar2ExpArg) /
					Math.sqrt(detc2Covariance);
			double condExpArg = 0;
			double adjustedPCond = 0;
			if (L_cc != null) {
				condExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(condDeviationsFromMean,
							invCondCovariance),
					condDeviationsFromMean);
				adjustedPCond = Math.exp(-0.5 * condExpArg) /
					Math.sqrt(detccCovariance);
			}
			double jointExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(deviationsFromMean,
							invCovariance),
					deviationsFromMean);
			double adjustedPJoint = Math.exp(-0.5 * jointExpArg) /
					Math.sqrt(detCovariance);
			
			if (L_cc != null) {
				// Returning results in nats:
				localValues[t] = Math.log(adjustedPJoint * adjustedPCond /
						(adjustedPCondVar1 * adjustedPCondVar2));
			} else {
				// Return an MI (no linearly independent, non-zero conditional vars):
				localValues[t] = Math.log(adjustedPJoint /
						(adjustedPCondVar1 * adjustedPCondVar2));
			}

			if (biasCorrection) {
				// Remove the average bias from every local estimate.
				// Note that we do the same thing even if these are new observations,
				//  because the variances have been computed from the same
				//  number of samples.
				localValues[t] -= analyticMeasDist.getMeanOfUncorrectedDistribution();
			}

			if (isPreviousObservations) {
				lastAverage += localValues[t];
			}
		}
		
		if (isPreviousObservations) {
			// Originally we didn't store the average value here, worrying that it wouldn't be exactly 
			//  the same as what would have been computed under the analytic expression; 
			//  however it should be, any difference is only numerical
			lastAverage /= (double) newVar2Obs.length;
			condMiComputed = true;
		}
		
		return localValues;
	}
}
