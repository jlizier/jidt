package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * <p>An abstract interface for calculators computing measures from a source to a destination.
 * </p>
 * 
 * <p>The interface is abstract because further specification is required
 * for addObservations or setObservations methods showing whether
 * this is a univariate or multivariate measure.</p>
 * 
 * 
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 *
 */
public abstract interface ChannelCalculatorCommon {

	/**
	 * Allows the user to set properties for the underlying calculator implementation
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;

	/**
	 * Elect to add in the observations from several disjoint time series.
	 *
	 */
	public void startAddObservations();
	
	/**
	 * Flag that the observations are complete, probability distribution functions can now be built.
	 * @throws Exception 
	 *
	 */
	public void finaliseAddObservations() throws Exception;
	
	/**
	 * 
	 * @return the average value of the implemented channel measure,
	 *  computed using all of the previously supplied observation sets.
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations() throws Exception;

	/**
	 * <p>Computes the local values of the implemented channel measure,
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
	 * <p>Compute the significance of obtaining the given average TE from the given observations.</p>
	 * 
	 * <p>This is as per Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128.
	 * </p>
	 * 
	 * <p>Basically, we shuffle the source observations against the destination tuples.
	 * This keeps the marginal PDFs the same (including the entropy rate of the destination)
	 *  but destroys any correlation between the source and state change of the destination.
	 * </p>
	 * 
	 * @param numPermutationsToCheck number of new orderings of the source values to compare against
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception;
	
	/**
	 * <p>As per {@link computeSignificance(int) computeSignificance()} but supplies
	 *  the re-orderings of the observations of the source variables.</p>
	 * 
	 * @param newOrderings first index is permutation number, i.e. newOrderings[i]
	 * 		is an array of 1 permutation of 0..n-1, where there were n observations.
	 * 		If the length of each permutation in newOrderings
	 * 		is not equal to numObservations, an Exception is thrown. 
	 * @return
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception;

	public void setDebug(boolean debug);
	
	public double getLastAverage();

	public int getNumObservations() throws Exception;

}
