package infodynamics.measures.continuous;

public interface EntropyCalculator {

	/**
	 * Initialise the calculator using the existing or default value of calculator-specific parameters
	 * 
	 */
	public void initialise() throws Exception;
	
	/**
	 * Allows the user to set properties for the underlying calculator implementation
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;

	public void setObservations(double observations[]);
	
	public double computeAverageLocalOfObservations();

	public void setDebug(boolean debug);
}
