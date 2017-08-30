/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2017, Joseph T. Lizier, Ipek Oezdemir and Pedro Mediano
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
 * <p>Base interface for our information-theoretic calculators
 *  on continuous (double[]) data,
 *  providing common functionality
 *  for user-level measure classes.</p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * <ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>set properties via {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator using {@link #initialise()} or
 *			other initialise methods defined by child classes;</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			sets of "addObservations" methods defined by child classes, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average measure: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>or other quantities as defined by child classes.</li>
 * 			</ul>
 * 		</li>
 * 		<li>
 * 		Return to step 3 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface InfoMeasureCalculatorContinuous {

	/**
	 * Initialise the calculator for (re-)use, with the existing (or default) values of parameters.
	 * Clears any PDFs of previously supplied observations.
	 */
	public void initialise() throws Exception;
	
	/**
	 * Set properties for the underlying calculator implementation.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Property names are defined by the implementing class
	 * (generally in their documentation).
	 * Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
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
	 * Get the number of samples to be used for the PDFs here 
	 * which have been supplied by calls to
	 * "setObservations", "addObservations" etc.
	 * 	
	 * <p>Note that the number of samples may not be equal to the length of time-series
	 * supplied (e.g. for transfer entropy, where we need to accumulate
	 * a number of samples for the past history of the destination).
	 * </p>
	 * 
	 * @return the number of samples to be used for the PDFs
	 */
	public int getNumObservations() throws Exception;

	/**
	 * Compute the average value of the measure
	 * from the previously-supplied samples.
	 * 
	 * @return the estimate of the measure
	 *  (in bits or nats, depending on the estimator)
	 */
	public abstract double computeAverageLocalOfObservations() throws Exception;
	
	/**
	 * Return the measure last calculated in a call to
	 * {@link #computeAverageLocalOfObservations()}
	 * or related methods after the previous
	 * {@link #initialise()} call.
	 * 
	 * @return the last computed measure value
	 */
	public double getLastAverage();

	/**
	 * Set or clear debug mode for extra debug printing to stdout
	 * 
	 * @param debug new setting for debug mode (on/off)
	 */
	public void setDebug(boolean debug);
}
