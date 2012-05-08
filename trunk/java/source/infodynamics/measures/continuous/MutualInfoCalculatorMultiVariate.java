package infodynamics.measures.continuous;

public interface MutualInfoCalculatorMultiVariate extends ChannelCalculatorMultiVariate {

	/**
	 * Time difference between data1 and data2 - assumed to be >= 0
	 */ 
	public static final String PROP_TIME_DIFF = "TIME_DIFF";

	/**
	 * Compute the mutual information if the second variable were ordered as per the ordering
	 *  specified in newOrdering
	 * 
	 * @param newOrdering
	 * @return
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[] newOrdering) throws Exception;

	public double[] computeLocalUsingPreviousObservations(double states1[][], double states2[][])
		throws Exception;
}
