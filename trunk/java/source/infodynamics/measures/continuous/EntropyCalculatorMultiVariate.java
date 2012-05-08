package infodynamics.measures.continuous;

/**
 * Interface for computing entropy for
 *  a multi-variate values (independent of underlying technique)
 * 
 * 
 * @author Joseph Lizier
 *
 */
public interface EntropyCalculatorMultiVariate {

	public void initialise(int dimensions);
	
	public void setObservations(double observations[][]);
	
	public double computeAverageLocalOfObservations();
	
	public double[] computeLocalUsingPreviousObservations(double states[][]) throws Exception;

	public double[] computeLocalOfPreviousObservations();
	
	public double getLastAverage();
	
	public void setDebug(boolean debug);
	
}
