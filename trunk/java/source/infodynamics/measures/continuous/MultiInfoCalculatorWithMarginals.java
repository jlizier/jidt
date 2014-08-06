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
