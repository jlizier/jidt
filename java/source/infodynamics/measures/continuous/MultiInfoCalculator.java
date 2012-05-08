package infodynamics.measures.continuous;

public interface MultiInfoCalculator {

	public void initialise(int dimensions);

	/**
	 * Allows the user to set properties for the underlying calculator implementation
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;
	
	public void setObservations(double observations[][]) throws Exception;
	
	public void startIndividualObservations();
	
	public void addObservation(double observation[]);
	
	public void endIndividualObservations() throws Exception;

	public double computeAverageLocalOfObservations() throws Exception;

	public double[] computeLocalOfPreviousObservations()
		throws Exception;
	
	public double[] computeLocalUsingPreviousObservations(double states[][])
		throws Exception;

	public void setDebug(boolean debug);
	
	public double getLastAverage();
}
