package infodynamics.measures.mixed;

import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * A conditional mutual information calculator between a joint set of continuous variables, 
 *  and a discrete variable, conditioned on another continuous variable vector. 
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com
 *
 */
public interface ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSource {

	/**
	 * 
	 * 
	 * @param dimensions the number of joint continuous variables
	 * @param base the base of the discrete variable
	 * @param dimensionsCond the number of joint continuous variables
	 *     to condition on
	 * @throws Exception
	 */
	public void initialise(int dimensions, int base, int dimensionsCond) throws Exception;
	
	public void setProperty(String propertyName, String propertyValue);
	
	public void startAddObservations();

	public void finaliseAddObservations() throws Exception;

	public void addObservations(double[][] continuousObservations,
			int[] discreteObservations, double[][] conditionedObservations) throws Exception;

	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations, double[][] conditionedObservations) throws Exception;

	public double computeAverageLocalOfObservations() throws Exception;

	public double[] computeLocalOfPreviousObservations() throws Exception;

	public double[] computeLocalUsingPreviousObservations(double[][] contStates,
			int[] discreteStates, double[][] conditionedStates) throws Exception;

	/**
	 * 
	 * @param reorderDiscreteVariable boolean for whether to reorder the discrete variable (true)
	 *    or the continuous variable (false)
	 * @param numPermutationsToCheck
	 * @return
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(boolean reorderDiscreteVariable, int numPermutationsToCheck) throws Exception;

	/**
	 * 
	 * @param reorderDiscreteVariable boolean for whether to reorder the discrete variable (true)
	 *    or the continuous variable (false)
	 * @param newOrderings
	 * @return
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(boolean reorderDiscreteVariable, int[][] newOrderings) throws Exception;

	public void setDebug(boolean debug);

	public double getLastAverage();
	
	public int getNumObservations();

}
