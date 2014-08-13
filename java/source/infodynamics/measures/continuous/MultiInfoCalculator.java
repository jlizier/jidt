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
 * <p>Interface for implementations of the multi-information or integration,
 * which may be applied to either multivariate or merely univariate
 * continuous data.
 * That is, it is applied to <code>double[][]</code> data, where the first index
 * is observation number or time, and the second is variable number.
 * See Tononi et al. below for the definition of multi-information/integration.</p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator using {@link #initialise(int)};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[][])}
 * 					for calculations based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to {@link #addObservations(double[][])} or
 * 							{@link #addObservation(double[])}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average multi-information:
 *					{@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local multi-information values for these samples:
 * 					{@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local multi-information values for a specific set of samples:
 * 					{@link #computeLocalUsingPreviousObservations(double[][])}</li>
 * 			</ul>
 * 		</li>
 * 		<li>
 * 		Return to step 2 or 3 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>"G. Tononi, O. Sporns, G. M. Edelman,
 *  <a href="http://dx.doi.org/10.1073/pnas.91.11.5033">"A measure for
 *  brain complexity:
 * 	relating functional segregation and integration in the nervous system"</a>
 *  Proceedings of the National Academy of Sciences, Vol. 91, No. 11.
 *  (1994), pp. 5033-5037.</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface MultiInfoCalculator {

	/**
	 * Initialise the calculator for (re-)use, with the existing
	 * (or default) values of parameters, with number of 
	 * joint variables specified.
	 * Clears an PDFs of previously supplied observations.
	 *
	 * @param dimensions the number of joint variables to consider
	 */
	public void initialise(int dimensions);

	/**
	 * Set properties for the calculator.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>No general properties are defined at the interface level here, i.e.
	 * there are only properties defined by child interfaces and classes.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;
	
	/**
	 * Sets a single series from which to compute the PDF for the multi-information.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * <p>The supplied series may be a time-series, or may be simply
	 * a set of separate observations
	 * without a time interpretation.</p>
	 * 
	 * <p>Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated).</p> 
	 * 
	 * @param observations series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 * @throws Exception
	 */
	public void setObservations(double observations[][]) throws Exception;
	
	/**
	 * Signal that we will add in the samples for computing the PDF 
	 * from several disjoint time-series or trials via calls to
	 * {@link #addObservation(double[])} or {@link #addObservations(double[][])}
	 * rather than {@link #setDebug(boolean)}.
	 */
	public void startAddObservations();
	
	/**
	 * <p>Adds a new (single) observation to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p>Note that the arrays must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param observation a single multivariate observation
	 *  (index is variable number)
	 */
	public void addObservation(double observation[]);
	
	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p>Note that the arrays must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param observations series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 */
	public void addObservations(double observations[][]);
	
	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
	 * 
	 * @throws Exception 
	 */
	public void finaliseAddObservations() throws Exception;

	/**
	 * Compute the multi-information from the previously-supplied samples.
	 * 
	 * @return the estimate of the multi-information
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations() throws Exception;

	/**
	 * <p>Computes the local values of the multi-information,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>If the samples were supplied via a single call i.e.
	 * {@link #setObservations(double[][])},
	 * then the return value is a single <i>time-series</i> of local
	 * channel measure values corresponding to these samples.</p>
	 * 
	 * <p>Otherwise where disjoint time-series observations were supplied using several 
	 *  calls such as {@link addObservations(double[][])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single "time-series" return array.</p>
	 *  
	 * @return the "time-series" of local multi-information values.
	 * @throws Exception
	 */
	public double[] computeLocalOfPreviousObservations()
		throws Exception;
	
	/**
	 * Compute the local multi-information values for each of the
	 * supplied samples in <code>states</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>states</code>
	 * (unless they were
	 * some of the previously supplied samples).</p>
	 * 
	 * @param states series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 * @return the series of local multi-information values.
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double states[][])
		throws Exception;

	/**
	 * Set or clear debug mode for extra debug printing to stdout
	 * 
	 * @param debug new setting for debug mode (on/off)
	 */
	public void setDebug(boolean debug);
	
	/**
	 * Return the multi-information last calculated in a call to {@link #computeAverageLocalOfObservations()}
	 * or {@link #computeLocalOfPreviousObservations()} after the previous
	 * {@link #initialise(int)} call.
	 * 
	 * @return the last computed multi-information value
	 */
	public double getLastAverage();
}
