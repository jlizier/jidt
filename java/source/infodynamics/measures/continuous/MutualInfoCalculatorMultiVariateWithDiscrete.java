package infodynamics.measures.continuous;

import infodynamics.utils.MeasurementDistribution;

public interface MutualInfoCalculatorMultiVariateWithDiscrete {

	public void initialise(int dimensions, int base) throws Exception;
	
	public void setProperty(String propertyName, String propertyValue);
	
	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations) throws Exception;

	public double computeAverageLocalOfObservations() throws Exception;

	public double[] computeLocalUsingPreviousObservations(double[][] contStates, int[] discreteStates) throws Exception;

	public MeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception;

	public MeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception;

	public void setDebug(boolean debug);

	public double getLastAverage();
	
	public int getNumObservations();

}
