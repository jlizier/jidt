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

/**
 * Interface for calculators of information-theoretic measures
 * for single variables (e.g. entropy, active information storage).
 * The interface defines common operations such as
 * adding observations and calculating
 * local and average values, etc. 
 * 
 * <p>Usage is as per {@link InfoMeasureCalculatorDiscrete}, with
 * many methods for supplying observations and making 
 * calculations defined here.</p>
 * 
 * <p>It would ideally be an abstract class to be inherited from, but
 * it's more important for some of our calculators to have inheritance from
 * ContextOfPastCalculator, and since java doesn't allow multiple
 * inheritance, one of them has to miss out.
 * To get around this, we combine the two in
 * {@link UnivariateMeasureDiscreteInContextOfPastCalculator}.
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface UnivariateMeasureDiscrete {

	/**
	 * Property name for definite alphabet size. Attempting
	 * 	to use observations beyond this property will throw
	 * 	errors.
	 */
	public static final String ALPHABET_SIZE = "ALPHABET_SIZE";

	/**
	 * Property name for if the alphabet size is known.
	 */
	public static final String KNOWN_INTEGER_RANGE = "KNOWN_INTEGER_RANGE";

	/**
	 * Initialise the calculator with same or unknown alphabet size
	 */
	public void initialise();

	/**
	 * Sets the samples from which to compute the PDF for the entropy.
	 * Should only be called once, the last call contains the
	 * 	observations that are used (they are not accumulated). 
	 * 
	 * @param observations array of (univariate) samples
	 * @throws Exception
	 */
	public void setObservations(Object[] observations) throws Exception;

	/**
	 * Signal that we will add in the samples for computing the PDF 
	 * from several disjoint time-series or trials via calls to
	 * "addObservations" rather than "setObservations" type methods
	 * (defined by the child interfaces and classes).
	 */
	public void startAddObservations();

	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
	 * 
	 * @throws Exception when the estimator has no observations.
	 */
	public void finaliseAddObservations() throws Exception;
	
	/**
 	 * Add observations in to our estimates of the pdfs.
	 * Univariate.
	 *
	 * @param states series of samples
	 */
	public void addObservations(int[] states);
	
	
	/**
	 * Compute the average value of the measure
	 * from the previously-supplied samples.
	 * 
	 * Must set average, min and max
	 * 
	 * @return the estimate of the measure
	 */
	public double computeAverageLocalOfObservations();

	/**
	 * Computes local information theoretic measure for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 *  
	 * Must set average, min and max
	 * 
	 * @param states time series of samples
	 * @return time-series of local values (indexed as per states)
	 */
	public double[] computeLocalFromPreviousObservations(int states[]);
	
	/**
	 * Standalone routine to 
	 * compute the average information theoretic measure across a time-series
	 *  of states.
	 * Return the average.
	 * 
	 * @param states time series of samples
	 * @return average of the information-theoretic measure.
	 * @throws Exception if states parameter is empty.
	 */
	public double computeAverageLocal(int states[]) throws Exception;
	
	/**
	  * Add observations in to our estimates of the pdfs.
	  * This call suitable only for homogeneous agents, as all
	  *  agents will contribute to the PDFs.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 */
	@Deprecated
	public void addObservations(int states[][]) throws Exception;
	
	/**
	 * Initialise the calculator with (potentially) a new alphabet size
	 * 
	 * Deprecated, use initialise() and setProperty("ALPHABET_SIZE", "n") instead.
	 * @deprecated
	 * @param alphabetSize
	 */
	@Deprecated
	public void initialise(int alphabetSize);
}
