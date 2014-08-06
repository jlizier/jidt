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

import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * <p>An abstract interface for calculators computing measures from a source to a destination.
 * </p>
 * 
 * <p>The interface is abstract because further specification is required
 * for addObservations or setObservations methods showing whether
 * this is a univariate or multivariate measure.</p>
 * 
 * 
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 *
 */
public abstract interface ChannelCalculatorCommon {

	/**
	 * Initialise the calculator for re-use with new observations.
	 * All parameters remain unchanged.
	 * 
	 * @throws Exception
	 */
	public void initialise() throws Exception;

	/**
	 * Allows the user to set properties for the underlying calculator implementation.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;

	/**
	 * Elect to add in the observations from several disjoint time series.
	 *
	 */
	public void startAddObservations();
	
	/**
	 * Flag that the observations are complete, probability distribution functions can now be built.
	 * @throws Exception 
	 *
	 */
	public void finaliseAddObservations() throws Exception;
	
	/**
	 * Return whether the user has added more than a single observation set via the
	 *  startAddObservations, addObservations, finaliseAddObservations sequence.
	 * 
	 * @return
	 */
	public boolean getAddedMoreThanOneObservationSet();
	
	/**
	 * 
	 * @return the average value of the implemented channel measure,
	 *  computed using all of the previously supplied observation sets.
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations() throws Exception;

	/**
	 * <p>Computes the local values of the implemented channel measure,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>If disjoint observations were supplied using several 
	 *  calls such as {@link ChannelCalculator#addObservations(double[], double[])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single return array,
	 *  though of course the time series for these disjoint observations were
	 *  not appended in computing the required PDFs).</p>
	 *  
	 * @return array of local values.
	 * @throws Exception
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception;

	/**
	 * <p>Compute the significance of obtaining the given average 
 	 *    measure from the given observations.</p>
	 * 
	 * <p>This is as per Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128.
	 * </p>
	 * 
	 * <p>Basically, we shuffle the source observations against the destination tuples.
	 * This keeps the marginal PDFs the same (including the entropy rate of the destination)
	 *  but destroys any correlation between the source and state change of the destination.
	 * </p>
	 * 
	 * @param numPermutationsToCheck number of new orderings of the source values to compare against
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception;
	
	/**
	 * <p>As per {@link computeSignificance(int) computeSignificance()} but supplies
	 *  the re-orderings of the observations of the source variables.</p>
	 * 
	 * @param newOrderings first index is permutation number, i.e. newOrderings[i]
	 * 		is an array of 1 permutation of 0..n-1, where there were n observations.
	 * 		If the length of each permutation in newOrderings
	 * 		is not equal to numObservations, an Exception is thrown. 
	 * @return
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception;

	/**
	 * Set whether to print debug messages or not
	 * 
	 * @param debug whether to print debug messages or not
	 */
	public void setDebug(boolean debug);
	
	/**
	 * Get the last computed average of the measure
	 * 
	 * @return the last computed average
	 */
	public double getLastAverage();

	/**
	 * Get the number of observations that have been supplied for
	 *  computation of the PDFs
	 * 
	 * @return number of observations
	 * @throws Exception
	 */
	public int getNumObservations() throws Exception;

}
