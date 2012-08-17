package infodynamics.measures.continuous;

/**
 * Interface for computing entropy for
 *  a multi-variate values (independent of underlying technique)
 * 
 * 
 * @author Joseph Lizier, joseph.lizier_at_gmail.com
 *
 */
public interface EntropyCalculatorMultiVariate {

	public void initialise(int dimensions);
	
	/**
	 * Supply the observations for which to compute the PDFs for the entropy
	 * 
	 * @param observations multivariate time series of observations; first index
	 *  is time step, second index is variable number (total should match dimensions
	 *  supplied to {@link #initialise(int)}
	 * @throws Exception if the dimensions of the observations do not match 
	 *  the expected value supplied in {@link #initialise(int)}; implementations
	 *  may throw other more specific exceptions also.
	 */
	public void setObservations(double observations[][]) throws Exception;
	
	public double computeAverageLocalOfObservations();
	
	public double[] computeLocalUsingPreviousObservations(double states[][]) throws Exception;

	public double[] computeLocalOfPreviousObservations() throws Exception;
	
	public double getLastAverage();
	
	public void setDebug(boolean debug);
	
}
