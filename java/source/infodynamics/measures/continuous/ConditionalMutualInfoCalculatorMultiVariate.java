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

package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.EmpiricalNullDistributionComputer;

/**
 * <p>Interface for implementations of the <b>conditional mutual information</b>,
 * which may be applied to either multivariate or merely univariate
 * continuous data.
 * That is, it is applied to <code>double[][]</code> data, where the first index
 * is observation number or time, and the second is variable number.</p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)}
 * 			which may now include the property {@link #PROP_TIME_DIFF}
 * 			to specify a source-destination time difference (default 0);</li>
 *		<li>Initialise the calculator using
 *			{@link #initialise(int, int, int)};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[][], double[][], double[][])} or
 * 					{@link #setObservations(double[][], double[][], double[][], boolean[], boolean[], boolean[])} or
 * 					{@link #setObservations(double[][], double[][], double[][], boolean[][], boolean[][], boolean[][])}
 * 					for calculations based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to
 * 							{@link #addObservations(double[][], double[][], double[][])} or
 * 							{@link #addObservations(double[][], double[][], double[][], int, int)}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average MI: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local MI values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local MI values for a specific set of samples:
 * 				{@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])} </li>
 * 				<li>the distribution of MI values under the null hypothesis
 * 					of no relationship between source and
 * 					destination values: {@link #computeSignificance(int, int)} or
 * 					{@link #computeSignificance(int, int[][])}.</li>
 * 			</ul>
 * 		</li>
 * 		<li>
 * 		Return to step 2 or 3 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see "T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991)."
 */
public interface ConditionalMutualInfoCalculatorMultiVariate
	extends InfoMeasureCalculatorContinuous, EmpiricalNullDistributionComputer {

	/**
	 * Property name for 
	 *  specifying whether the data is normalised or not (to mean 0,
	 *  variance 1, for each of the multiple variables)
	 *  before the calculation is made.
	 */
	public static final String PROP_NORMALISE = "NORMALISE";
	/**
	 * Property name for an amount of random Gaussian noise to be
	 *  added to the data (default 0, except for Kraskov/KSG estimator
	 *  where it is 1e-8, matching the MILCA toolkit)
	 */
	public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";
	/**
	 * Property name for the seed for the random number generator for noise to be
	 *  added to the data (default is no seed)
	 */
	public static final String PROP_NOISE_SEED = "NOISE_SEED";
	/**
	 * Property value to indicate no seed for the random number generator for noise to be
	 *  added to the data
	 */
	public static final String NOISE_NO_SEED_VALUE = "NONE";

	/**
	 * Initialise the calculator for (re-)use, clearing PDFs,
	 * with the existing
	 * (or default) values of parameters (except the numbers
	 * of joint variables) as specified below
	 *
	 * @param var1Dimensions the number of joint variables in variable 1
	 * @param var2Dimensions the number of joint variables in variable 2
	 * @param condDimensions the number of joint variables in the conditional
	 */
	public void initialise(int var1Dimensions, int var2Dimensions, int condDimentions) throws Exception;

	/**
	 * Sets a single series from which to compute the PDF.
	 * Cannot be called in conjunction with
	 * {@link #startAddObservations()} / {@link #addObservations(double[][], double[][], double[][])} or
	 * {@link #addObservations(double[][], double[][], double[][], int, int)} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * <p>The supplied series may be (multivariate) time-series or 
	 * simply a set of separate observations without a time interpretation.
	 * 
	 * @param var1 multivariate observations for variable 1
	 *  (first index is time or observation index, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond multivariate observations for the conditional
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond) throws Exception;

	/**
	 * As per {@link #setObservations(double[][], double[][], double[][])}
	 * but where dimensions of all variables have been set to 1 
	 * 
	 * @param var1 univariate observations for variable 1
	 * @param var2 univariate observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond univariate observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void setObservations(double[] var1, double[] var2,
			double[] cond) throws Exception;

	/**
	 * As per {@link #setObservations(double[][], double[][], double[][])}
	 * but where dimensions of var2 and cond are 1 
	 * 
	 * @param var1 observations for variable 1
	 * @param var2 univariate observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond univariate observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void setObservations(double[][] var1, double[] var2,
			double[] cond) throws Exception;

	/**
	 * As per {@link #setObservations(double[][], double[][], double[][])}
	 * but where dimensions of var1 and cond are 1 
	 * 
	 * @param var1 univariate observations for variable 1
	 * @param var2 observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond univariate observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void setObservations(double[] var1, double[][] var2,
			double[] cond) throws Exception;

	/**
	 * As per {@link #setObservations(double[][], double[][], double[][])}
	 * but where dimensions of var1 and var2 are 1 
	 * 
	 * @param var1 univariate observations for variable 1
	 * @param var2 univariate observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void setObservations(double[] var1, double[] var2,
			double[][] cond) throws Exception;

	/**
	 * As per {@link #setObservations(double[][], double[][], double[][])}
	 * but where dimension of cond is 1 
	 * 
	 * @param var1 observations for variable 1
	 * @param var2 observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond univariate observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void setObservations(double[][] var1, double[][] var2,
			double[] cond) throws Exception;

	/**
	 * As per {@link #setObservations(double[][], double[][], double[][])}
	 * but where dimension of var2 is 1 
	 * 
	 * @param var1 observations for variable 1
	 * @param var2 univariate observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void setObservations(double[][] var1, double[] var2,
			double[][] cond) throws Exception;

	/**
	 * As per {@link #setObservations(double[][], double[][], double[][])}
	 * but where dimension of var1 is 1 
	 * 
	 * @param var1 univariate observations for variable 1
	 * @param var2 observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void setObservations(double[] var1, double[][] var2,
			double[][] cond) throws Exception;

	/**
	 * Sets a single series from which to compute the PDF,
	 *  where all the various observations are valid.
	 * Cannot be called in conjunction with
	 * {@link #startAddObservations()} / {@link #addObservations(double[][], double[][], double[][])} or
	 * {@link #addObservations(double[][], double[][], double[][], int, int)} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * <p>The supplied series may be (multivariate) time-series or 
	 * simply a set of separate observations without a time interpretation.
	 * 
	 * @param var1 multivariate observations for variable 1
	 *  (first index is time or observation index, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond multivariate observations for the conditional
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param var1Valid series (with indices the same as var1)
	 *  indicating whether var1 at that point is valid.
	 * @param var2Valid series (with indices the same as var2)
	 *  indicating whether var2 at that point is valid.
	 * @param condValid series (with indices the same as cond)
	 *  indicating whether cond at that point is valid.
	 */
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[] var1Valid, boolean[] var2Valid,
			boolean[] condValid) throws Exception;

	/**
	 * Sets a single series from which to compute the PDF,
	 * where the observations are valid for all individual variables.
	 * Cannot be called in conjunction with
	 * {@link #startAddObservations()} / {@link #addObservations(double[][], double[][], double[][])} or
	 * {@link #addObservations(double[][], double[][], double[][], int, int)} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * <p>The supplied series may be (multivariate) time-series or 
	 * simply a set of separate observations without a time interpretation.
	 * 
	 * @param var1 multivariate observations for variable 1
	 *  (first index is time or observation index, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond multivariate observations for the conditional
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param var1Valid series (with indices the same as var1)
	 *  indicating whether each variable of var1 at that point is valid.
	 * @param var2Valid series (with indices the same as var2)
	 *  indicating whether each variable of var2 at that point is valid.
	 * @param condValid series (with indices the same as cond)
	 *  indicating whether each variable of cond at that point is valid.
	 */
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[][] var1Valid, boolean[][] var2Valid,
			boolean[][] condValid) throws Exception;

	/**
	 * Signal that we will add in the samples for computing the PDF 
	 * from several disjoint time-series or trials via calls to
	 * "addObservations" rather than "setObservations" type methods
	 * (defined by the child interfaces and classes).
	 *
	 */
	public void startAddObservations();
	
	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p>Note that the arrays must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param var1 multivariate observations for variable 1
	 *  (first index is time or observation index, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond multivariate observations for the conditional
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond) throws Exception;
	
	/**
	 * <p>As per {@link #addObservations(double[][], double[][], double[][])}
	 * but can only be used where dimensions of all
	 * variables have been set to 1.</p>
	 * 
	 * @param var1 univariate observations for variable 1
	 * @param var2 univariate observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond univariate observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void addObservations(double[] var1, double[] var2,
			double[] cond) throws Exception;
	
	/**
	 * <p>As per {@link #addObservations(double[][], double[][], double[][])}
	 * but can only be used where dimensions of
	 * var2 and cond have been set to 1.</p>
	 * 
	 * @param var1 observations for variable 1
	 * @param var2 univariate observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond univariate observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void addObservations(double[][] var1, double[] var2,
			double[] cond) throws Exception;

	/**
	 * <p>As per {@link #addObservations(double[][], double[][], double[][])}
	 * but can only be used where dimensions of
	 * var1 and cond have been set to 1.</p>
	 * 
	 * @param var1 univariate observations for variable 1
	 * @param var2 observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond univariate observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void addObservations(double[] var1, double[][] var2,
			double[] cond) throws Exception;

	/**
	 * <p>As per {@link #addObservations(double[][], double[][], double[][])}
	 * but can only be used where dimensions of
	 * var1 and var2 have been set to 1.</p>
	 * 
	 * @param var1 univariate observations for variable 1
	 * @param var2 univariate observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void addObservations(double[] var1, double[] var2,
			double[][] cond) throws Exception;

	/**
	 * <p>As per {@link #addObservations(double[][], double[][], double[][])}
	 * but can only be used where dimensions of
	 * var1 have been set to 1.</p>
	 * 
	 * @param var1 univariate observations for variable 1
	 * @param var2 observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void addObservations(double[] var1, double[][] var2,
			double[][] cond) throws Exception;

	/**
	 * <p>As per {@link #addObservations(double[][], double[][], double[][])}
	 * but can only be used where dimensions of
	 * var2 have been set to 1.</p>
	 * 
	 * @param var1 observations for variable 1
	 * @param var2 univariate observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void addObservations(double[][] var1, double[] var2,
			double[][] cond) throws Exception;

	/**
	 * <p>As per {@link #addObservations(double[][], double[][], double[][])}
	 * but can only be used where dimensions of
	 * cond have been set to 1.</p>
	 * 
	 * @param var1 observations for variable 1
	 * @param var2 observations for variable 2
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond univariate observations for the conditional
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void addObservations(double[][] var1, double[][] var2,
			double[] cond) throws Exception;

	/**
	 * <p>Adds a new sub-series of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p>Note that the arrays must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param var1 multivariate observations for variable 1
	 *  (first index is time or observation index, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond multivariate observations for the conditional
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use (including startTime)
	 * @throws Exception
	 */
	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond,
			int startTime, int numTimeSteps) throws Exception;

	/**
	 * <p>As per {@link #addObservations(double[][], double[][], duoble[][])};
     *  but also includes parameters to track which observation set
     *  the samples came from. Intended to only be used by other 
     *  estimator classes here and not by users directly.</p>
	 * 
	 * @param var1 multivariate observations for variable 1
	 *  (first index is time or observation index, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param cond multivariate observations for the conditional
	 *  (first index is time or observation index, second is variable number)
	 *  Length must match <code>var1</code>, and their indices must correspond.
	 * @param observationSetIndexToUse which set of observations these came fmor
	 * @param startTimeIndex which was the first time index of these
	 *  samples within that observation set.
	 * @throws Exception
	 */
	public void addObservationsTrackObservationIDs(double[][] var1, double[][] var2, double[][] cond,
			int observationSetIndexToUse, int startTimeIndex) throws Exception;

	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
	 * 
	 * @throws Exception 
	 */
	public void finaliseAddObservations() throws Exception;
	
	/**
	 * <p>Computes the local values of the conditional mutual information,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>If the samples were supplied via a single call such as
	 * {@link #setObservations(double[][], double[][], double[][])},
	 * then the return value is a single time-series of local
	 * channel measure values corresponding to these samples.</p>
	 * 
	 * <p>Otherwise where disjoint time-series observations were supplied using several 
	 *  calls such as {@link #addObservations(double[][], double[][], double[][])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single "time-series" return array.</p>
	 *  
	 * @return the "time-series" of local conditional MI values in either bits or nats
	 *  depending on the estimator.
	 * @throws Exception
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception;

	/**
	 * <p><b>Note</b> -- in contrast to the generic computeSignificance() method 
	 * described in the documentation below, this method for a conditional MI calculator currently fixes the relationship
	 * between variable 2 and the conditional, and shuffles
	 * variable 1 with respect to these.
	 * To shuffle variable 2 instead, call {@link #computeSignificance(int, int)}.
	 * </p>
	 * 
	 * @inheritDoc
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception;
	
	/**
	 * Defined as per {@link #computeSignificance(int)} except that this method allows
	 * the user to specify which variable is shuffled (whilst the other has its
	 * relationship with the conditional preserved).
	 * 
	 * @param variableToReorder which variable to shuffle:
	 * 	1 for variable 1, 2 for variable 2.
	 * @param numPermutationsToCheck number of surrogate samples to bootstrap
	 *  to generate the distribution.
	 * @return the distribution of channel measure scores under this null hypothesis.
	 * @see {@link #computeSignificance(int)}
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int variableToReorder,
			int numPermutationsToCheck) throws Exception;
	
	/**
	 * <p><b>Note</b> -- in contrast to the generic computeSignificance() method 
	 * described in the documentation below, this method for a conditional MI calculator currently fixes the relationship
	 * between variable 2 and the conditional, and shuffles
	 * variable 1 with respect to these.
	 * To shuffle variable 2 instead, call {@link #computeSignificance(int, int[][])}.
	 * </p>
	 * 
	 * @inheritDoc
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception;
	
	/**
	 * Defined as per {@link #computeSignificance(int[][])} except that this method allows
	 * the user to specify which variable is shuffled (whilst the other has its
	 * relationship with the conditional preserved).
	 * 
	 * @param variableToReorder which variable to shuffle:
	 * 	1 for variable 1, 2 for variable 2.
	 * @param newOrderings a specification of how to shuffle the values
	 *  of the variable specified by <code>variableToReorder</code>
	 *  to create the surrogates to generate the distribution with. The first
	 *  index is the permutation number (i.e. newOrderings.length is the number
	 *  of surrogate samples we use to bootstrap to generate the distribution here.)
	 *  Each array newOrderings[i] should be an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 * @return the distribution of channel measure scores under this null hypothesis.
	 * @see {@link #computeSignificance(int[][])} except that this method all 
	 * @throws Exception where the length of each permutation in newOrderings
	 *   is not equal to the number N samples that were previously supplied.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int[][] newOrderings) throws Exception;

	/**
	 * Compute the conditional mutual information if the given variable
	 *  were ordered as per the ordering
	 *  specified in newOrdering
	 * 
	 * @param variableToReorder which variable to shuffle:
	 * 	1 for variable 1, 2 for variable 2.
	 * @param newOrdering permutation of the indices for the given variable;
	 *  must be an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 * @return conditional MI under this reordering
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int variableToReorder, int[] newOrdering) throws Exception;

	/**
	 * Compute the local conditional MI values for each of the
	 * supplied samples in <code>states1</code>, <code>states2</code>
	 * and <code>condStates</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>states1</code>, <code>states2</code>
	 * and <code>condStates</code>
	 * (unless they were
	 * some of the previously supplied samples).</p>
	 * 
	 * @param states1 series of multivariate observations for variable 1
	 *  (first index is time or observation index, second is variable number)
	 * @param states2 series of multivariate observations for variable 2
	 *  (first index is time or observation index, second is variable number).
	 *  Length must match <code>states1</code>, and their indices must correspond.
	 * @param condStates series of multivariate observations for the conditional
	 *  (first index is time or observation index, second is variable number).
	 *  Length must match <code>states1</code>, and their indices must correspond.
	 * @return the series of local conditional MI values.
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double states1[][], double states2[][], double[][] condStates)
		throws Exception;

	/**
	 * <p>As per {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}
	 * but can only be used where dimensions of
	 * states1 and states2 have been set to 1.</p>
	 * 
	 * @param states1 series of univariate observations for variable 1
	 * @param states2 series of univariate observations for variable 2
	 *  Length must match <code>states1</code>, and their indices must correspond.
	 * @param condStates series of multivariate observations for the conditional
	 *  (first index is time or observation index, second is variable number).
	 *  Length must match <code>states1</code>, and their indices must correspond.
	 * @return the series of local conditional MI values.
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double states1[], double states2[], double[][] condStates)
		throws Exception;

	/**
	 * <p>As per {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}
	 * but can only be used all dimensions have been set to 1.</p>
	 * 
	 * @param states1 series of univariate observations for variable 1
	 * @param states2 series of univariate observations for variable 2.
	 *  Length must match <code>states1</code>, and their indices must correspond.
	 * @param condStates series of univariate observations for the conditional.
	 *  Length must match <code>states1</code>, and their indices must correspond.
	 * @return the series of local conditional MI values.
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double states1[], double states2[], double[] condStates)
		throws Exception;

	/**
	 * @throws Exception if the implementing class computes MI without
	 * explicit observations (e.g. see
	 * {@link infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian})
	 */
	@Override
	public int getNumObservations() throws Exception;
	
	/**
	 * Get whether the user has added more than one observation set
	 *  via the {@link #addObservations(double[][], double[][], double[][])}
	 *  and {@link #addObservations(double[][], double[][], double[][], int, int)}
	 *  methods.
	 * 
	 * @return whether the user has added more than one observation set.
	 */
	public boolean getAddedMoreThanOneObservationSet();

	/**
	 * Retrieve an array indicating which observation set each sample came from
	 * 
	 * @return array of integers
	 */
	public int[] getObservationSetIndices();

	/**
	 * Retrieve an array indicating which time index within its observation set that sample came from
	 * 
	 * @return array of integers
	 */
	public int[] getObservationTimePoints();
}
