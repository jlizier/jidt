package infodynamics.measures.continuous;

/**
 * An interface for calculators computing measures from a source to a destination
 * for single variables in each.
 * 
 * 
 * @author Joseph Lizier, jlizier at gmail.com
 *
 */
public interface ChannelCalculator extends ChannelCalculatorCommon {

	public void initialise() throws Exception;

	public void setObservations(double source[], double destination[]) throws Exception;
	
	/**
	 * Add some more observations.
	 * Note that the arrays source and destination must not be over-written by the user
	 *  until after finaliseAddObservations() has been called.
	 * 
	 * @param source
	 * @param destination
	 */
	public void addObservations(double[] source, double[] destination) throws Exception;

	/**
	 * Add some more observations.
	 * 
	 * @param source
	 * @param destination
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 */
	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception ;

	/**
	 * Sets the observations to compute the PDFs from.
	 * Cannot be called in conjunction with start/add/finaliseAddObservations.
	 * destValid is a time series (with time indices the same as destination)
	 *  indicating whether the destination at that point is valid.
	 * sourceValid is the same for the source
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param sourceValid
	 * @param destValid
	 */
	public void setObservations(double[] source, double[] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception;

}
