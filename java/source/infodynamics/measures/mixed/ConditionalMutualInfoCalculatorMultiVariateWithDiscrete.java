package infodynamics.measures.mixed;

import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * A conditional mutual information calculator between a joint set of continuous variables, 
 *  and a discrete variable, conditioned on another discrete variable. 
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com
 *
 */
public interface ConditionalMutualInfoCalculatorMultiVariateWithDiscrete {

	/**
	 * 
	 * 
	 * @param dimensions the number of joint continuous variables
	 * @param base the base of the discrete variable
	 * @param condBase the base of the discrete variable to condition on
	 * @throws Exception
	 */
	public void initialise(int dimensions, int base, int condBase) throws Exception;
	
	public void setProperty(String propertyName, String propertyValue);
	
	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations, int[] conditionedObservations) throws Exception;

	public double computeAverageLocalOfObservations() throws Exception;

	public double[] computeLocalUsingPreviousObservations(double[][] contStates,
			int[] discreteStates, int[] conditionedStates) throws Exception;

	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception;

	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception;

	public void setDebug(boolean debug);

	public double getLastAverage();
	
	public int getNumObservations();

}
