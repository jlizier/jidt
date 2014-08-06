/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
	 * @param sourceIndex
	 * @param destIndex
	 */
	public void addObservations(int states[][], int sourceIndex, int destIndex);

	/**
	 * Add observations for the source and destination
	 * @param source
	 * @param dest
	 */
	public void addObservations(int[] source, int[] dest);
	
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
