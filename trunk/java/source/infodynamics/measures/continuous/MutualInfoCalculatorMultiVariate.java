package infodynamics.measures.continuous;

/**
 * <p>Interface for multivariate implementations of the mutual information.</p>
 *  
 * <p>
 * Intended usage of the child classes:
 * 	<ol>
 * 		<li>Construct</li>
 * 		<li>Set properties using {@link #setProperty(String, String)}</li>
 *		<li>{@link #initialise(int, int)} or {@link #initialise(int, int, double)}</li>
 * 		<li>Provide the observations to the calculator using:
 * 			{@link #setObservations(double[][], double[][])}, or
 * 			{@link #setCovariance(double[][])}, or
 * 			a sequence of:
 * 			{@link #startAddObservations()},
 *          multiple calls to either {@link #addObservations(double[][], double[][])}
 *          or {@link #addObservations(double[][], double[][], int, int)}, and then
 *          {@link #finaliseAddObservations()}.</li>
 * 		<li>Compute the required information-theoretic results, primarily:
 * 			{@link #computeAverageLocalOfObservations()} to return the average
 *          value based on the supplied observations; or other calls to compute
 *          local values or statistical significance.</li>
 * 	</ol>
 * </p>
 * 
 * @author Joseph Lizier, <a href="mailto:joseph.lizier at gmail.com">joseph.lizier at gmail.com</>
 *
 * @see "T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991)."
 */
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
