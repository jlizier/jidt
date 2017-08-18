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

import infodynamics.utils.EmpiricalNullDistributionComputer;

/**
 * Interface for implementations of
 * Predictive Information / Excess Entropy estimators on
 * continuous univariate data (ie double[] arrays).
 * 
 * <p>See definition of Predictive Information (PI) / Excess Entropy in the citations below.
 * Basically, PI is the mutual information between the past <i>state</i>
 * of a time-series process <i>X</i> and its future <i>state</i>.
 * The past <i>state</i> at time <code>n</code>
 * is represented by an embedding vector of <code>k</code> values from <code>X_n</code> backwards,
 * each separated by <code>\tau</code> steps, giving
 * <code><b>X^k_n</b> = [ X_{n-(k-1)\tau}, ... , X_{n-\tau}, X_n]</code>.
 * The future <i>state</i> from time <code>n+1</code>
 * is represented by an embedding vector of <code>k</code> values from <code>X_{n+1}</code> forwards,
 * each separated by <code>\tau</code> steps, giving
 * <code><b>X^{k+}_{n+1}</b> = [ X_{n+1}, X_{n+1+\tau}, ... , X_{n+1+(k-1)\tau}]</code>.
 * We call <code>k</code> the embedding dimension, and <code>\tau</code>
 * the embedding delay.
 * PI is then the mutual information between <b>X^k_n</b> and <b>X^{k+}_{n+1}</b>.</p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator using {@link #initialise()} or
 *			{@link #initialise(int)} or {@link #initialise(int, int)};
 *		</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[])} or
 * 					{@link #setObservations(double[], boolean[])} for calculations
 * 					based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to {@link #addObservations(double[])} or
 * 							{@link #addObservations(double[], int, int)}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul>
 * 		</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average PI: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local PI values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local PI values for a specific set of samples: {@link #computeLocalUsingPreviousObservations(double[])}</li>
 * 				<li>the distribution of PI values under the null hypothesis
 * 					of no relationship between past sequences in the series
 * 					and the next sequence: {@link #computeSignificance(int)} or
 * 					{@link #computeSignificance(int[][])}.</li>
 * 			</ul>
 * 		</li>
 * 		<li>
 * 		Return to step 2 or 3 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Bialek, W., Nemenman, I., and Tishby, N.,
 *  <a href="http://dx.doi.org/10.1016/S0378-4371(01)00444-7">
 * 	"Complexity through nonextensivity"</a>,
 *  Physica A, 302, 89-99. (2001).</li>
 * 	<li>J. P. Crutchfield, D. P. Feldman,
 *  <a href="http://dx.doi.org/10.1063/1.1530990">
 * 	"Regularities Unseen, Randomness Observed: Levels of Entropy Convergence"</a>,
 *  Chaos, Vol. 13, No. 1. (2003), pp. 25-54.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface PredictiveInfoCalculator
	extends InfoMeasureCalculatorContinuous, EmpiricalNullDistributionComputer {

	/**
	 * Property name for embedding length <code>k</code> of
	 * the past and future vectors (1 by default).
	 */
	public static final String K_EMBEDDING_PROP_NAME = "k_EMBEDDING";
	/**
	 * Alternative property name for embedding length <code>k</code> of
	 * the past and future vectors (1 by default);
	 * a second option, in order to provide compatibility with 
	 * {@link ActiveInfoStorageCalculator#K_PROP_NAME}
	 */
	public static final String K_PROP_NAME = "k_HISTORY";
	/**
	 * Property name for embedding delay <code>\tau</code> of the vectors
	 * (1 by default).
	 * The delay exists between points in the k-vectors
	 *  but there is always only one time step between the
	 *  past k-vector and the next k-vector.  
	 */
	public static final String TAU_PROP_NAME = "TAU";

	/**
	 * Initialise the calculator for (re-)use, with some parameters
	 * supplied here, and existing (or default) values of other parameters
	 * to be used.
	 * 
	 * @param k embedding length of past and future vectors
	 */
	public void initialise(int k) throws Exception;
	
	/**
	 * Initialise the calculator for (re-)use, with some parameters
	 * supplied here, and existing (or default) values of other parameters
	 * to be used.
	 * 
	 * @param k embedding length of past and future vectors
	 * @param tau embedding delay of past and future vectors
	 */
	public void initialise(int k, int tau) throws Exception;

	/**
	 * Set properties for the underlying calculator implementation.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Property names defined at the interface level, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 	<li>{@link #K_PROP_NAME} or {@link #K__EMBEDDING_PROP_NAME} --
	 * 		embedding length <code>k</code> of
	 * 		the past and future vectors</li>
	 *  <li>{@link #TAU_PROP_NAME} -- embedding delay between each of the
	 *     <code>k</code> points in the vectors.</li>
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
	 * Sets a single time-series from which to compute the PDF for the PI.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * @param observations time-series array of (univariate) samples,
	 *  where the array index is time.
	 */
	public void setObservations(double observations[]) throws Exception;
	
	/**
	 * Signal that we will add in the samples for computing the PDF 
	 * from several disjoint time-series or trials via calls to
	 * {@link #addObservations(double[])} etc.
	 *
	 */
	public void startAddObservations();
	
	/**
	 * Add more time-series for the computation of the PDF.
	 * The array observations must not be over-written by the user
	 *  until after finaliseAddObservations() has been called.
	 * 
	 * @param observations time-series array of (univariate) samples,
	 *  where the array index is time.
	 */
	public void addObservations(double[] observations) throws Exception;

	/**
	 * Add more time-series for the computation of the PDF, using
	 * only a sub-series of <code>observations</code>.
	 * The array observations must not be over-written by the user
	 *  until after finaliseAddObservations() has been called.
	 * 
	 * @param observations time-series array of (univariate) samples,
	 *  where the array index is time.
	 * @param startTime first time index to extract samples from
	 * @param numTimeSteps number of time steps to extract starting from startTime
	 */
	public void addObservations(double[] observations,
			int startTime, int numTimeSteps) throws Exception ;

	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
	 * 
	 * @throws Exception 
	 */
	public void finaliseAddObservations() throws Exception;
	
	/**
	 * Sets a single time-series from which to compute the PDF for the AIS,
	 * subject to the validity of each sample in that series.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 *  
	 * @param observations time-series array of (univariate) samples,
	 *  where the array index is time.
	 * @param valid a time series (with indices the same as observations)
	 *  indicating whether the entry in observations at that index is valid; 
	 *  we only take vectors as samples to add to the observation set where
	 *  all points in the time series (even between points in 
	 *  the embedded k-vector with embedding delays) are valid.
	 */
	public void setObservations(double[] observations,
			boolean[] valid) throws Exception;

	/**
	 * Compute the local PI values for each of the
	 * previously-supplied samples.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations.</p>
	 * 
	 * <p>If the samples were supplied via a single call such as
	 * {@link #setObservations(double[])},
	 * then the return value is a single time-series of local
	 * PI values corresponding to these samples (with the first, and last,
	 * <code>(k-1)*tau + 1</code>
	 * values set to 0 since PI is undefined there).</p>
	 * 
	 * <p>Otherwise where disjoint time-series observations were supplied using several 
	 *  calls such as {@link addObservations(double[])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single "time-series" return array (without any <code>(k-1)*tau + 1</code>
	 *  leading or lagging 0 values).</p>
	 *  
	 * @return the "time-series" of local AIS values.
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception;

	/**
	 * Compute the local PI values for each of the
	 * supplied samples in <code>newObservations</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>newObservations</code> (unless they were
	 * some of the previously supplied samples).</p>
	 * 
	 * @param newObservations time-series for which to compute local PI values
	 * @return time-series of local PI values corresponding to newObservations
	 *   (the first, and last, <code>(k-1)*tau + 1</code>
	 *    values are set to 0 since PI is undefined there)
	 * @throws Exception for an invalid input array, e.g. of length
	 *  less than <code>2*k</code>
	 */
	public double[] computeLocalUsingPreviousObservations(double[] newObservations) throws Exception;

}
