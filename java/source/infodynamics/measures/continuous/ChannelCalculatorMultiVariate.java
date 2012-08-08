package infodynamics.measures.continuous;

/**
 * An interface for calculators computing measures from a source to a destination
 * for multiple joint variables in each (specific examples intended here
 * are the mutual information and transfer entropy).
 * 
 * 
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 *
 */
public interface ChannelCalculatorMultiVariate extends ChannelCalculatorCommon {

	/**
	 * Initialise the calculator
	 *
	 * @param sourceDimensions the number of joint variables in the source
	 * @param destDimensions the number of joint variables in the destination
	 */
	public void initialise(int sourceDimensions, int destDimensions) throws Exception;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param source multivariate observations for the source variable
	 *  (first index is time, second is variable number)
	 * @param destination multivariate observations for the destination variable
	 *  (first index is time, second is variable number)
	 * @throws Exception
	 */
	public void setObservations(double[][] source, double[][] destination) throws Exception;

	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * <p>Note that the arrays source and destination must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source multivariate observations for the source variable
	 *  (first index is time, second is variable number)
	 * @param destination multivariate observations for the destination variable
	 *  (first index is time, second is variable number)
	 * @throws Exception
	 */
	public void addObservations(double[][] source, double[][] destination) throws Exception;
	
	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * <p>Note that the arrays source and destination must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source multivariate observations for the source variable
	 *  (first index is time, second is variable number)
	 * @param destination multivariate observations for the destination variable
	 *  (first index is time, second is variable number)
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 * @throws Exception
	 */
	public void addObservations(double[][] source, double[][] destination,
			int startTime, int numTimeSteps) throws Exception;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param source multivariate observations for the source variable
	 *  (first index is time, second is variable number)
	 * @param destination multivariate observations for the destination variable
	 *  (first index is time, second is variable number)
	 * @param sourceValid time series (with time indices the same as source)
	 *  indicating whether the source at that point is valid.
	 * @param destValid time series (with time indices the same as destination)
	 *  indicating whether the destination at that point is valid.
	 */
	public void setObservations(double[][] source, double[][] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param source multivariate observations for the source variable
	 *  (first index is time, second is variable number)
	 * @param destination multivariate observations for the destination variable
	 *  (first index is time, second is variable number)
	 * @param sourceValid time series (with time indices the same as source)
	 *  indicating whether each variable of the source at that point is valid.
	 * @param destValid time series (with time indices the same as destination)
	 *  indicating whether each variable of the destination at that point is valid.
	 */
	public void setObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception;

}
