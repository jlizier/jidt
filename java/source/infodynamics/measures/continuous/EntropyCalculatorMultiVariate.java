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
 * Interface for implementations of
 * entropy estimators on continuous multivariate data
 * (ie <code>double[][]</code> arrays
 * where first index is time or observation number, second is variable index).
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator {@link #initialise(int)};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      to set up the PDFs, using:
 * 			{@link #setObservations(double[][])};</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average entropy: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local entropy values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local entropy values for a specific set of samples:
 * 					{@link #computeLocalUsingPreviousObservations(double[])}.</li>
 * 			</ul>
 * 		</li>
 * 		<li>Return to step 2 or 3 to re-use the calculator on a new data set.</li>
 * 	</ol>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991).</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface EntropyCalculatorMultiVariate 
	extends EntropyCalculator {

	/**
	 * Property name for the number of dimensions
	 */
	public static final String NUM_DIMENSIONS_PROP_NAME = "NUM_DIMENSIONS";
	/**
	 * Property for whether we normalise the incoming observations to mean 0,
	 * standard deviation 1.
	 */
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	/**
	 * Property name for an amount of random Gaussian noise to be
	 *  added to the data (default is 1e-8, matching the MILCA toolkit).
	 */
	public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";
	/**
	 * Property name for the seed for the random number generator for noise to be
	 *  added to the data (default is no seed)
	 */
	public static final String PROP_NOISE_SEED = "NOISE_SEED";
	/**
	 * Property value to indicate no seed for the random number generator for noise to be
	 *  added to the data
	 */
	public static final String NOISE_NO_SEED_VALUE = "NONE";
	
	/**
	 * Set properties for the underlying calculator implementation.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Property names defined at the interface level, and what their
	 * values should represent, include:</p>
	 * <ul>
	 *  <li>{@link #NUM_DIMENSIONS_PROP_NAME} -- number of dimensions in the joint
	 * 		variable that we are computing the entropy of.</li>
	 * 	<li>{@link #NORMALISE_PROP_NAME} -- whether to normalise the incoming variable values
	 * 			to mean 0, standard deviation 1, or not (default false). Sets {@link #normalise}.</li>
	 *  <li>{@link #PROP_ADD_NOISE} -- a standard deviation for an amount of
	 *      random Gaussian noise to add to
	 *      each variable, to avoid having neighbourhoods with artificially
	 *      large counts. (We also accept "false" to indicate "0".)
	 *      The amount is added in after any normalisation,
	 *      so can be considered as a number of standard deviations of the data.
	 *      Default is 0 for most estimators; this is strongly recommended by
	 *      by Kraskov for the KSG method though, so for the Kozachenko estimator we 
	 *      use 1e-8 to match the MILCA toolkit (though note it adds in
	 *      a random amount of noise in [0,noiseLevel) ).</li>
	 *  <li>{@link #PROP_NOISE_SEED} -- a long value seed for the random noise generator or
	 *      the string {@link MutualInfoCalculatorMultiVariate#NOISE_NO_SEED_VALUE} for no seed (default)</li>
	 * </ul>
	 *  
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * <p>Note that implementing classes may defined additional properties.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception;
	
	/**
	 * Initialise the calculator for (re-)use, with the existing (or default) values
	 * of calculator-specific parameters.
	 * Clears any PDFs of previously supplied observations.
	 * 
	 * @param dimensions number of joint variables to be investigated
	 */
	public void initialise(int dimensions);
	
	/**
	 * Signal that we will add in the samples for computing the PDF 
	 * from several disjoint time-series or trials via calls to
	 * "addObservations" rather than "setObservations" type methods
	 * (defined by the child interfaces and classes).
	 */
	public void startAddObservations();
	
	/**
	 * Add more observations for which to compute the PDFs for the entropy.
	 * May be called multiple times between {@link #startAddObservations()} and
	 * {@link #finaliseAddObservations()}.
	 * 
	 * @param observations multivariate time series of observations; first index
	 *  is time step, second index is variable number (total should match dimensions
	 *  supplied to {@link #initialise(int)}
	 * @throws Exception if the dimensions of the observations do not match 
	 *  the expected value supplied in {@link #initialise(int)}; implementations
	 *  may throw other more specific exceptions also.
	 */
	public void addObservations(double[][] observations) throws Exception;

	/**
	 * Add more samples from which to compute the PDF for the entropy.
	 * Only allowed to be called when set up for dimension == 1.
	 * May be called multiple times between {@link #startAddObservations()} and
	 * {@link #finaliseAddObservations()}.
	 * 
	 * @param observations array of (univariate) samples
	 * @throws Exception if the expected dimensions were not 1.
	 */
	public void addObservations(double[] observations) throws Exception;

	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
	 * 
	 * @throws Exception 
	 */
	public void finaliseAddObservations() throws Exception;
	
	/**
	 * Set the observations for which to compute the PDFs for the entropy
	 * Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated).
	 * 
	 * @param observations multivariate time series of observations; first index
	 *  is time step, second index is variable number (total should match dimensions
	 *  supplied to {@link #initialise(int)}
	 * @throws Exception if the dimensions of the observations do not match 
	 *  the expected value supplied in {@link #initialise(int)}; implementations
	 *  may throw other more specific exceptions also.
	 */
	public void setObservations(double observations[][]) throws Exception;
	
	/**
	 * Compute the local entropy values for each of the
	 * supplied samples in <code>newObservations</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>newObservations</code> (unless they were
	 * some of the previously supplied samples).</p>
	 * 
	 * @param newObservations multivariate time-series for which to compute
	 * 	local entropy values (see {@link #setObservations(double[][])} for its format)
	 * @return time-series of local entropy values corresponding to each entry
	 * in <code>newObservations</code>
	 */
	public double[] computeLocalUsingPreviousObservations(double newObservations[][]) throws Exception;

}
