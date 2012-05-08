package infodynamics.measures.continuous;

/**
 * An interface for calculators computing measures from a source to a destination
 * for multiple joint variables in each.
 * 
 * 
 * @author Joseph Lizier, jlizier at gmail.com
 *
 */
public interface ChannelCalculatorMultiVariate extends ChannelCalculatorCommon {

	/**
	 * Initialise the calculator
	 *
	 */
	public void initialise(int destDimensions, int sourceDimensions) throws Exception;

	public void setObservations(double[][] source, double[][] destination) throws Exception;

	public void addObservations(double[][] source, double[][] destination) throws Exception;
	
	public void addObservations(double[][] source, double[][] destination,
			int startTime, int numTimeSteps) throws Exception;

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
	public void setObservations(double[][] source, double[][] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception;

	/**
	 * Sets the observations to compute the PDFs from.
	 * Cannot be called in conjunction with start/add/finaliseAddObservations.
	 * destValid is a time series (with time indices the same as destination)
	 *  indicating whether the destination at that point is valid, with one
	 *  validity supplied for each sub-variable.
	 * sourceValid is the same for the source
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param sourceValid
	 * @param destValid
	 */
	public void setObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception;

}
