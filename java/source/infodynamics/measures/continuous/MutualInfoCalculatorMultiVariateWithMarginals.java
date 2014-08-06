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

public interface MutualInfoCalculatorMultiVariateWithMarginals
	extends MutualInfoCalculatorMultiVariate {

	public double computeAverageJointEntropy();

	public double computeAverageEntropyOfObservation1();

	public double computeAverageEntropyOfObservation2();

	public double computeAverageInfoDistanceOfObservations();

	public double[] computeLocalJointEntropyOfPreviousObservations()
		throws Exception;

	public double[] computeLocalJointEntropyUsingPreviousObservations(double states1[][], double states2[][])
		throws Exception;

	public double[] computeLocalEntropy1OfPreviousObservations();

	public double[] computeLocalEntropy1UsingPreviousObservations(double states[][])
		throws Exception;

	public double[] computeLocalEntropy2OfPreviousObservations();

	public double[] computeLocalEntropy2UsingPreviousObservations(double states[][])
		throws Exception;

	public double[] computeLocalInfoDistanceOfPreviousObservations()
		throws Exception;

	public double[] computeLocalInfoDistanceUsingPreviousObservations(double states1[][], double states2[][])
		throws Exception;

}
