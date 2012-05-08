package infodynamics.measures.discrete;

/**
 * Interface to define adding observations and calculating
 * local and average values of info theoretic measures
 * for single agent metrics (e.g. entropy, active information).
 * Would ideally be an abstract class to be inherited from, but
 * it's more important for us to have inheritance from
 * ContextOfPastCalculator, and since java doesn't allow multiple
 * inheritance, one of them has to miss out.
 * 
 * @author Joseph Lizier
 *
 */
public interface SingleAgentMeasure {

	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states 1st index is time, 2nd index is agent number
	 */
	public void addObservations(int states[][]);

	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd index is agent number
	 * @param col index of agent
	 */
	public void addObservations(int states[][], int col);

	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 */
	public void addObservations(int states[][][]);

	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 * @param index1 first index of agent
	 * @param index2 index of agent in 2nd dimension
	 */
	public void addObservations(int states[][][], int index1, int index2);

	/**
	 * Returns the average information theoretic measure from
	 *  the observed values which have been passed in previously. 
	 *  
	 * Must set average, min and max
	 * 
	 * @return
	 */
	public double computeAverageLocalOfObservations();

	/**
	 * Computes local information theoretic measure for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method to be used for homogeneous agents only
	 *  
	 * Must set average, min and max
	 * 
	 * @param states 1st index is time, 2nd index is agent number
	 * @return
	 */
	public double[][] computeLocalFromPreviousObservations(int states[][]);

	/**
	 * Computes local information theoretic measure for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method is suitable for heterogeneous agents
	 *  
	 * Must set average, min and max
	 * 
	 * @param states 1st index is time, 2nd index is agent number
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][], int col);

	/**
	 * Computes local information theoretic measure for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method to be used for homogeneous agents only
	 *  
	 * Must set average, min and max
	 * 
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 * @return
	 */
	public double[][][] computeLocalFromPreviousObservations(int states[][][]);

	/**
	 * Computes local information theoretic measure for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method is suitable for heterogeneous agents
	 *  
	 * Must set average, min and max
	 * 
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][][], int index1, int index2);

	/**
	 * Standalone routine to 
	 * compute local information theoretic measure across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * 
	 * @param states 1st index is time, 2nd index is agent number
	 * @return
	 */
	public double[][] computeLocal(int states[][]);

	/**
	 * Standalone routine to 
	 * compute local information theoretic measure across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * 
	 * @param states 1st index is time, 2nd and 3rd index are agent number
	 * @return
	 */
	public double[][][] computeLocal(int states[][][]);

	/**
	 * Standalone routine to 
	 * compute average local information theoretic measure across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return the average
	 * This method to be called for homogeneous agents only
	 * 
	 * @param states -1st index is time, 2nd index is agent number
	 * @return
	 */
	public double computeAverageLocal(int states[][]);

	/**
	 * Standalone routine to 
	 * compute local information theoretic measure for one agent in a 2D spatiotemporal
	 *  array of the states of agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method should be used for heterogeneous agents
	 * 
	 * @param states 1st index is time, 2nd index is agent number
	 * @param col - column number of the agent in the states array
	 * @return
	 */
	public double[] computeLocalAtAgent(int states[][], int col);
	
	/**
	 * Standalone routine to 
	 * compute average local information theoretic measure 
	 * for a single agent
	 * Returns the average
	 * This method suitable for heterogeneous agents
	 * @param blocksize - Size of blocks to compute entropy over
	 * @param base - base of the states
	 * @param states 1st index is time, 2nd index is agent number
	 * @param col - column number of the agent in the states array
	 * 
	 * @return
	 */
	public double computeAverageLocalAtAgent(int states[][], int col);

}
