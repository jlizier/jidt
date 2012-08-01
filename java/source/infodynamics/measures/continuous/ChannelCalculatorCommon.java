package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * An abstract  interface for calculators computing measures from a source to a destination.
 * 
 * 
 * @author Joseph Lizier, jlizier at gmail.com
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
	 *
	 */
	public void finaliseAddObservations();
	
	public double computeAverageLocalOfObservations() throws Exception;

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

	public int getNumObservations();

}
