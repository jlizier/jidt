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

package infodynamics.measures.spiking;

import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * <p>Interface for implementations of the <b>transfer entropy</b> (TE),
 * which may be applied to spiking time-series data.
 * That is, it is applied to <code>double[]</code> data, as an array
 * of time stamps at which spikes were recorded.
 * See Schreiber below for the definition of transfer entropy,
 * and Lizier et al. for the definition of local transfer entropy,
 * and (To be published) for how to measure transfer entropy on spike trains.
 * Specifically, this class implements the pairwise or <i>apparent</i>
 * transfer entropy; i.e. we compute the transfer that appears to
 *  come from a single source variable, without examining any other
 *  potential sources
 *  (see Lizier et al, PRE, 2008).</p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)}
 * 			e.g. including properties describing
 * 			the source and destination embedding;</li>
 *		<li>Initialise the calculator using
 *			{@link #initialise()} or {@link #initialise(int, int)};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[], double[])} 
 * 					for calculations based on single recordings, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to
 * 							{@link #addObservations(double[], double[])}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul></li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average TE: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local TE values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>the distribution of TE values under the null hypothesis
 * 					of no relationship between source and
 * 					destination values: {@link #computeSignificance(int)} or
 * 					{@link #computeSignificance(int[][])}.</li>
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
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 *  <li>To be published</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface TransferEntropyCalculatorSpiking {

	/**
	 * Property name to specify the destination history embedding length k
	 * (default value 1)
	 */
	public static final String K_PROP_NAME = "k_HISTORY";
	/**
	 * Property name for embedding length for the source past history vector
	 * (default value 1)
	 */
	public static final String L_PROP_NAME = "l_HISTORY";
	/* Could try to do this one later. I think we would just consider the history 
	 * of source as being up to this many units of time behind destination and
	 * no later. 
	 * 
	 * Property name for source-destination delay (default value is 0)
	 *
	public static final String DELAY_PROP_NAME = "DELAY";
	 */
	/**
	 * Property name for whether each series of time stamps of spikes is sorted
	 *  into temporal order. (default true)
	 */
	public static final String TIMESORTED_PROP_NAME = "TIME_SORTED";
	
	/**
	 * Initialise the calculator for re-use with new observations.
	 * All parameters remain unchanged.
	 * 
	 * @throws Exception
	 */
	public void initialise() throws Exception;
	
	/**
	 * Initialise the calculator for re-use with new observations.
	 * A new history embedding length k can be supplied here; all other parameters
	 * remain unchanged.
	 * 
	 * @param k destination history embedding length to be considered.
	 * @throws Exception
	 */
	public void initialise(int k) throws Exception;
	
	/**
	 * Initialise the calculator for re-use with new observations.
	 * New history embedding lengths k and l can be supplied here; all other parameters
	 * remain unchanged.
	 * 
	 * @param k destination history embedding length to be considered.
	 * @param l source history embedding length to be considered.
	 * @throws Exception
	 */
	public void initialise(int k, int l) throws Exception;

	/**
	 * Sets properties for the TE calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #K_PROP_NAME} -- destination history embedding length k
	 * 			(default value 1)</li>
	 * 		<li>{@link #L_PROP_NAME} -- embedding length for the source past history vector
	 * 			(default value 1)</li>
	 * 		<li>{@link #TIMESORTED_PROP_NAME} --  whether each series of time stamps of spikes is sorted
	 *  		into temporal order. (default "true")</li>
	 * </ul>
	 * <p><b>Note:</b> further properties may be defined by child classes.</p>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value, 
	 * or if the property is recognised but unsupported (eg some
	 * calculators do not support all of the embedding properties).
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;

	/**
	 * Get current property values for the calculator.
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, are the same as those for
	 * {@link #setProperty(String, String)}</p>
	 * 
	 * <p>Unknown property values are responded to with a null return value.</p>
	 * 
	 * @param propertyName name of the property
	 * @return current value of the property
	 * @throws Exception for invalid property values
	 */
	public String getProperty(String propertyName) throws Exception;
	
	/**
	 * Sets a single set of spiking observations from which to compute the PDF for transfer entropy.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * @param source series of time stamps of spikes for the source variable. 
	 * @param destination series of time stamps of spikes for the destination
	 *  variable. Length will generally be different to the <code>source</code>, 
	 *  unlike other transfer entropy implementations, e.g. for {@link infodynamics.measures.continuous.TransferEntropyCalculator}.
	 *   <code>source</code> and  <code>destination</code> must have the same reference time point,
	 *   and each array is assumed to be sorted into temporal order unless
	 *   the property {@link #TIMESORTED_PROP_NAME} has been set to false.
	 *  
	 * @throws Exception
	 */
	public void setObservations(double source[], double destination[]) throws Exception;
	
	/**
	 * Signal that we will add in the samples for computing the PDF 
	 * from several disjoint time-series or trials via calls to
	 * {@link #addObservations(double[], double[])} rather than 
	 * {@link #setObservations(double[], double[])} type methods
	 * (defined by the child interfaces and classes).
	 *
	 */
	public void startAddObservations();
	
	/**
	 * <p>Adds a new set of spiking observations to update the PDFs with.
	 * It is intended to be called multiple times, and must
	 * be called after {@link #startAddObservations()}. Call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append or overlay these observations to the previously
	 *  supplied observations, but treats them as independent trials - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * <p>Note that the arrays source and destination must not be over-written by the user
	 *  until after {@link #finaliseAddObservations()} has been called
	 *  (they are not copied by this method necessarily, the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source series of time stamps of spikes for the source variable. 
	 *   Will be returned in ascending sorted order.
	 * @param destination series of time stamps of spikes for the destination
	 *  variable. Length will generally be different to the <code>source</code>, 
	 *  unlike other transfer entropy implementations, e.g. for {@link infodynamics.measures.continuous.TransferEntropyCalculator}.
	 *   <code>source</code> and  <code>destination</code> must have the same reference time point,
	 *   and each array is assumed to be sorted into temporal order unless
	 *   the property {@link #TIMESORTED_PROP_NAME} has been set to false.
	 *   Will be returned in ascending sorted order.
	 * @throws Exception
	 */
	public void addObservations(double[] source, double[] destination) throws Exception;

	/**
	 * Signal that the observations are now all added via
	 * {@link #addObservations(double[], double[])}, PDFs can now be constructed.
	 * 
	 * @throws Exception 
	 */
	public void finaliseAddObservations() throws Exception;
	
	/**
	 * Query whether the user has added more than a single observation set via the
	 *  {@link #startAddObservations()}, "addObservations" (defined by child interfaces
	 *  and classes), {@link #finaliseAddObservations()} sequence.
	 * 
	 * @return true if more than a single observation set was supplied
	 */
	public boolean getAddedMoreThanOneObservationSet();
	
	/**
	 * Compute the TE from the previously-supplied samples.
	 * 
	 * @return the estimate of the channel measure
	 */
	public double computeAverageLocalOfObservations() throws Exception;

	/**
	 * This interface serves to indicate the return type of {@link TransferEntropyCalculator#computeLocalOfPreviousObservations()}
	 *  as each child implementation will return something specific
	 * 
	 * @author Joseph Lizier
	 *
	 */
	public interface SpikingLocalInformationValues {
		// Left empty intentionally
	}
	
	/**
	 * @return an object containing a representation of 
	 *   the of local TE values. The precise contents of this representation
	 *   will vary depending on the underlying implementation
	 */
	public SpikingLocalInformationValues computeLocalOfPreviousObservations() throws Exception;

	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck, double estimatedValue) throws Exception;

	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck,
								    double estimatedValue, long randomSeed) throws Exception;
	/**
	 * Set or clear debug mode for extra debug printing to stdout
	 * 
	 * @param debug new setting for debug mode (on/off)
	 */
	public void setDebug(boolean debug);
	
	/**
	 * Return the TE last calculated in a call to {@link #computeAverageLocalOfObservations()}
	 * or {@link #computeLocalOfPreviousObservations()} after the previous
	 * {@link #initialise()} call.
	 * 
	 * @return the last computed channel measure value
	 */
	public double getLastAverage();
}
