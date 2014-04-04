package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * <p>Interface for multivariate implementations of the
 *    conditional mutual information.</p>
 *  
 * <p>
 * Intended usage of the child classes:
 * 	<ol>
 * 		<li>Construct</li>
 * 		<li>Set properties using {@link #setProperty(String, String)}</li>
 *		<li>{@link #initialise(int, int)} or {@link #initialise(int, int, double)}</li>
 * 		<li>Provide the observations to the calculator using:
 * 			{@link #setObservations(double[][], double[][])}, or
 * 			{@link #setCovariance(double[][])}, or
 * 			a sequence of:
 * 			{@link #startAddObservations()},
 *          multiple calls to either {@link #addObservations(double[][], double[][])}
 *          or {@link #addObservations(double[][], double[][], int, int)}, and then
 *          {@link #finaliseAddObservations()}.</li>
 * 		<li>Compute the required information-theoretic results, primarily:
 * 			{@link #computeAverageLocalOfObservations()} to return the average
 *          value based on the supplied observations; or other calls to compute
 *          local values or statistical significance.</li>
 * 	</ol>
 * </p>
 * 
 * @author Joseph Lizier, <a href="mailto:joseph.lizier at gmail.com">joseph.lizier at gmail.com</>
 *
 * @see "T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991)."
 */
public interface ConditionalMutualInfoCalculatorMultiVariate {

	/**
	 * Initialise the calculator
	 *
	 * @param var1Dimensions the number of joint variables in variable 1
	 * @param var2Dimensions the number of joint variables in variable 2
	 * @param condDimensions the number of joint variables in the conditional
	 */
	public void initialise(int var1Dimensions, int var2Dimensions, int condDimentions) throws Exception;

	/**
	 * Allows the user to set properties for the underlying calculator implementation
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[], double[])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param var1 multivariate observations for variable 1
	 *  (first index is time, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time, second is variable number)
	 * @param cond multivariate observations for the conditional
	 *  (first index is time, second is variable number)
	 * @throws Exception
	 */
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond) throws Exception;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param var1 multivariate observations for variable 1
	 *  (first index is time, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time, second is variable number)
	 * @param cond multivariate observations for the conditional
	 *  (first index is time, second is variable number)
	 * @param var1Valid time series (with time indices the same as var1)
	 *  indicating whether var1 at that point is valid.
	 * @param var2Valid time series (with time indices the same as var2)
	 *  indicating whether var2 at that point is valid.
	 * @param condValid time series (with time indices the same as cond)
	 *  indicating whether cond at that point is valid.
	 */
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[] var1Valid, boolean[] var2Valid,
			boolean[] condValid) throws Exception;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param var1 multivariate observations for variable 1
	 *  (first index is time, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time, second is variable number)
	 * @param cond multivariate observations for the conditional
	 *  (first index is time, second is variable number)
	 * @param var1Valid time series (with time indices the same as var1)
	 *  indicating whether each variable of var1 at that point is valid.
	 * @param var2Valid time series (with time indices the same as var2)
	 *  indicating whether each variable of var2 at that point is valid.
	 * @param condValid time series (with time indices the same as cond)
	 *  indicating whether each variable of cond at that point is valid.
	 */
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[][] var1Valid, boolean[][] var2Valid,
			boolean[][] condValid) throws Exception;

	/**
	 * Elect to add in the observations from several disjoint time series.
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
	 *  (first index is time, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time, second is variable number)
	 * @param cond multivariate observations for the conditional
	 *  (first index is time, second is variable number)
	 * @throws Exception
	 */
	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond) throws Exception;
	
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
	 *  (first index is time, second is variable number)
	 * @param var2 multivariate observations for variable 2
	 *  (first index is time, second is variable number)
	 * @param cond multivariate observations for the conditional
	 *  (first index is time, second is variable number)
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 * @throws Exception
	 */
	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond,
			int startTime, int numTimeSteps) throws Exception;

	/**
	 * Flag that the observations are complete, probability distribution functions can now be built.
	 * @throws Exception 
	 *
	 */
	public void finaliseAddObservations() throws Exception;
	
	/**
	 * 
	 * @return the average value of the conditional mutual information measure,
	 *  computed using all of the previously supplied observation sets.
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations() throws Exception;

	/**
	 * <p>Computes the local values of the conditional mutual information,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>If disjoint observations were supplied using several 
	 *  calls such as {@link ChannelCalculator#addObservations(double[], double[])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single return array,
	 *  though of course the time series for these disjoint observations were
	 *  not appended in computing the required PDFs).</p>
	 *  
	 * @return array of local values.
	 * @throws Exception
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception;

	/**
	 * <p>Compute the significance of obtaining the given average 
 	 *    measure from the given observations.</p>
	 * 
	 * <p>This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128.
	 * </p>
	 * 
	 * <p>Basically, we shuffle the observations of the named variable 
	 * against the other tuples.
	 * This keeps the marginal and joint PDFs of the unshuffled variables the same
	 *  but destroys any correlation between the named variable and the others.
	 * </p>
	 * 
	 * @param variableToReorder 1 for variable 1, 2 for variable 2
	 * @param numPermutationsToCheck number of new orderings of the source values to compare against
	 * @see "Chavez et. al., 'Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals', Journal of Neuroscience Methods 124 (2003) 113-128"
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int variableToReorder,
			int numPermutationsToCheck) throws Exception;
	
	/**
	 * <p>As per {@link #computeSignificance(int, int)} but supplies
	 *  the re-orderings of the observations of the named variable.</p>
	 * 
	 * @param variableToReorder 1 for variable 1, 2 for variable 2
	 * @param newOrderings first index is permutation number, i.e. newOrderings[i]
	 * 		is an array of 1 permutation of 0..n-1, where there were n observations.
	 * 		If the length of each permutation in newOrderings
	 * 		is not equal to numObservations, an Exception is thrown. 
	 * @return
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int[][] newOrderings) throws Exception;

	/**
	 * Compute the conditional mutual information if the given variable were ordered as per the ordering
	 *  specified in newOrdering
	 * 
	 * @param variableToReorder 1 for variable 1, 2 for variable 2
	 * @param newOrdering permutation of the indices for the given variable
	 * @return
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int variableToReorder, int[] newOrdering) throws Exception;

	/**
	 * Compute the local mutual information for the given states, using the
	 * PDFs from the previously supplied observations. 
	 * 
	 * @param states1
	 * @param states2
	 * @param condStates
	 * @return
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double states1[][], double states2[][], double[][] condStates)
		throws Exception;

	/**
	 * Set whether to print debug messages or not
	 * 
	 * @param debug whether to print debug messages or not
	 */
	public void setDebug(boolean debug);
	
	/**
	 * Get the last computed average of the measure
	 * 
	 * @return the last computed average
	 */
	public double getLastAverage();

	/**
	 * Get the number of observations that have been supplied for
	 *  computation of the PDFs
	 * 
	 * @return number of observations
	 * @throws Exception
	 */
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

}
