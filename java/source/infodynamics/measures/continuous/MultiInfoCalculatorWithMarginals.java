package infodynamics.measures.continuous;

/**
 * Adds methods for computing average and local marginal entropies, joint entropy, 
 *  and info distance.
 * 
 * @author Joseph Lizier
 *
 */
public interface MultiInfoCalculatorWithMarginals {

	public double computeAverageJointEntropy();

	public double computeAverageMarginalEntropy(int variableIndex);

	/**
	 * I'm not sure whether info distance is defined properly for multi-info in addition
	 *  to mutual info - I should check this.
	 * 
	 * @return
	 * @throws Exception
	 */
	public double computeAverageInfoDistanceOfObservations();

	public double[] computeLocalJointEntropyOfPreviousObservations()
		throws Exception;

	public double[] computeLocalJointEntropyUsingPreviousObservations(double states[][])
		throws Exception;

	public double[] computeLocalMarginalEntropyOfPreviousObservations(int variableIndex);

	public double[] computeLocalMarginalEntropyUsingPreviousObservations(double states[][], int variableIndex)
		throws Exception;

	/**
	 * I'm not sure whether info distance is defined properly for multi-info in addition
	 *  to mutual info - I should check this.
	 * 
	 * @return
	 * @throws Exception
	 */
	public double[] computeLocalInfoDistanceOfPreviousObservations()
		throws Exception;

	public double[] computeLocalInfoDistanceUsingPreviousObservations(double states[][])
		throws Exception;

}
