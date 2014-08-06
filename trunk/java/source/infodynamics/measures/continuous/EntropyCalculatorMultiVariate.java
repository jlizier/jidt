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
	 * Allows the user to set properties for the underlying calculator implementation
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;

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
	
	public int getNumObservations() throws Exception;
}
