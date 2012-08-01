package infodynamics.measures.discrete;

import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * An interface for calculators computing measures from a source to a destination.
 * 
 * 
 * @author Joseph Lizier, jlizier at gmail.com
 *
 */
public interface ChannelCalculator {

	/**
	 * Initialise the calculator
	 *
	 */
	public void initialise();
	
	/**
	 * Add observations for the source and destination
	 * 
	 * @param states
	 * @param destIndex
	 * @param sourceIndex
	 */
	public void addObservations(int states[][], int destIndex, int sourceIndex);

	/**
	 * Add observations for the source and destination
	 * 
	 * @param dest
	 * @param source
	 */
	public void addObservations(int[] dest, int[] source);
	
	/**
	 * Compute the value of the measure
	 * 
	 * @return
	 */
	public double computeAverageLocalOfObservations();
	
	/**
	 * Compute the significance of the average value for the channel measure here
	 * 
	 * @param numPermutationsToCheck
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck);
}
