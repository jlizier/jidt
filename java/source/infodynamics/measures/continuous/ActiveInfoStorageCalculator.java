package infodynamics.measures.continuous;

import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * Interface for calculators of the active information storage, as defined by
 * Lizier et al., 2012 (see below) 
 * 
 * @see J.T. Lizier, M. Prokopenko and A.Y. Zomaya, "Local measures of
 * information storage in complex distributed computation",
 * Information Sciences, vol. 208, pp. 39-54, 2012.
 * @see {@link http://dx.doi.org/10.1016/j.ins.2012.04.016}
 *
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 *
 */
public interface ActiveInfoStorageCalculator {

	/**
	 * Embedding length for the past history vector
	 */
	public static final String K_PROP_NAME = "k_HISTORY";
	/**
	 * Embedding delay for the past history vector
	 */
	public static final String TAU_PROP_NAME = "TAU";

	/**
	 * Initialise the calculator using the existing or default value of k
	 * 
	 */
	public void initialise() throws Exception;

	/**
	 * Initialise the calculator
	 * 
	 * @param k Length of past history to consider
	 */
	public void initialise(int k) throws Exception;
	
	/**
	 * Allows the user to set properties for the underlying calculator implementation
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;

	public void setObservations(double observations[]) throws Exception;
	
	/**
	 * Elect to add in the observations from several disjoint time series.
	 *
	 */
	public void startAddObservations();
	
	/**
	 * Add some more observations.
	 * Note that the array must not be over-written by the user
	 *  until after finaliseAddObservations() has been called.
	 * 
	 * @param observations
	 */
	public void addObservations(double[] observations) throws Exception;

	/**
	 * Add some more observations.
	 * 
	 * @param observations
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 */
	public void addObservations(double[] observations,
			int startTime, int numTimeSteps) throws Exception ;

	/**
	 * Flag that the observations are complete, probability distribution functions can now be built.
	 * @throws Exception 
	 *
	 */
	public void finaliseAddObservations() throws Exception;
	
	/**
	 * Sets the observations to compute the PDFs from.
	 * Cannot be called in conjunction with start/add/finaliseAddObservations.
	 * valid is a time series (with time indices the same as destination)
	 *  indicating whether the observation at that point is valid.
	 * 
	 * @param source observations for the source variable
	 * @param destValid
	 */
	public void setObservations(double[] observations,
			boolean[] valid) throws Exception;

	public double computeAverageLocalOfObservations() throws Exception;

	public double[] computeLocalOfPreviousObservations() throws Exception;

	public double[] computeLocalUsingPreviousObservations(double[] newObservations) throws Exception;

	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception;
	
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception;

	public void setDebug(boolean debug);
	
	public double getLastAverage();

	public int getNumObservations() throws Exception;
}
