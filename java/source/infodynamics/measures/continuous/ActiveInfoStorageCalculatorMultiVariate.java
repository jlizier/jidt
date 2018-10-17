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
 * Interface for implementations of Active Information Storage estimators on
 * multivariate continuous time-series data. That is, it is applied to
 * <code>double[][]</code> data, indexed first by time then by variable number.
 * 
 * <p>See definition of Active Information Storage (AIS) by Lizier et al. below.
 * Basically, AIS is the mutual information between the past <i>state</i>
 * of a time-series process <i>X</i> and its next value. The past <i>state</i> at time <code>n</code>
 * is represented by an embedding vector of <code>k</code> values from <code>X_n</code> backwards,
 * each separated by <code>\tau</code> steps, giving
 * <code><b>X^k_n</b> = [ X_{n-(k-1)\tau}, ... , X_{n-\tau}, X_n]</code>.
 * We call <code>k</code> the embedding dimension, and <code>\tau</code>
 * the embedding delay.
 * AIS is then the mutual information between <b>X^k_n</b> and X_{n+1}.</p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator using {@link #initialise(int)} or
 *			{@link #initialise(int, int)} or {@link #initialise(int, int, int)};
 *		</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[][])} for calculations
 * 					based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to {@link #addObservations(double[][])} or
 * 							{@link #addObservations(double[][], int, int)}, then</li>
 * 						<li>{@link #finaliseAddObservations()}; OR</li>
 * 					</ol></li>
 * 				<li>univariate calls as supported by {@link ActiveInfoStorageCalculator} </li>
 * 			</ul>
 * 		</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average AIS: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local AIS values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local AIS values for a specific set of samples: {@link #computeLocalUsingPreviousObservations(double[])}</li>
 * 				<li>the distribution of AIS values under the null hypothesis
 * 					of no relationship between past sequences in the series
 * 					and the next value: {@link #computeSignificance(int)} or
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
 * 	<li>J.T. Lizier, M. Prokopenko and A.Y. Zomaya,
 * 		<a href="http://dx.doi.org/10.1016/j.ins.2012.04.016">
 * 		"Local measures of information storage in complex distributed computation"</a>,
 * 		Information Sciences, vol. 208, pp. 39-54, 2012.</li>
 * </ul>
 * 
 * @author Pedro AM Mediano (<a href="pmediano at imperial.ac.uk">email</a>,
 * <a href="https://www.doc.ic.ac.uk/~pam213/">www</a>)
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 *
 * @see ActiveInfoStorageCalculator
 */
public interface ActiveInfoStorageCalculatorMultiVariate extends ActiveInfoStorageCalculator {

	/**
	 * Property name for the number of dimensions in the multivariate data.
	 */
	public static final String PROP_DIMENSIONS = "DIMENSIONS";

	/**
	 * Initialise the calculator for re-use with new observations.
	 * History length k, source and destination dimensions are
	 * specified here; all other parameters remain unchanged. 
	 * 
	 * @param dimensions number of joint variables in the system.
	 * @param k history embedding length to be considered.
	 * @param tau history embedding delay to be considered.
	 * @throws Exception
	 */
	public void initialise(int dimensions, int k, int tau) throws Exception;

	/**
	 * <p>Sets a single set of <b>univariate</b> observations to compute the PDFs from.
	 * Can only be called on this multivariate calculator if the dimension of the system
	 * is 1, otherwise throws an exception</p>
	 *  
	 * {@inheritDoc}
	 * 
	 * @throws Exception if initialised dimensions were not 1
	 */
	@Override
	public void setObservations(double observations[]) throws Exception;
	
	/**
	 * As per {@link #setObservations(double[])} only with multivariate 
	 *  observations
	 * 
	 * @param observations time-series array of (multivariate) samples,
	 *  where the first index is time.
	 * @throws Exception
	 */
	public void setObservations(double[][] observations) throws Exception;
	
	/**
	 * <p>Adds a new set of <b>univariate</b> observations to compute the PDFs from.
	 * Can only be called on this multivariate calculator if the dimension of the system
	 * is 1, otherwise throws an exception</p>
	 * 
	 * {@inheritDoc}
	 * 
	 * @throws Exception if initialised dimensions were not 1
	 */
	@Override
	public void addObservations(double[] observations) throws Exception;

	/**
	 * As per {@link #addObservations(double[])} only with multivariate 
	 *  observations
	 * 
	 * @param observations time-series array of (multivariate) samples,
	 *  where the first index is time.
	 * @throws Exception
	 */
	public void addObservations(double[][] observations) throws Exception;
	
	/**
	 * This method will return an exception if this class was not initialised 
	 *  for univariate data only.
	 *  
	 * {@inheritDoc}
	 * 
	 * @throws Exception if initialised dimensions were not 1
	 */
	@Override
	public void addObservations(double[] observations,
			int startTime, int numTimeSteps) throws Exception;

	/**
	 * As per {@link #addObservations(double[], int, int)} only with multivariate 
	 *  observations
	 * 
	 * @param observations time-series array of (multivariate) samples,
	 *  where the first index is time.
	 * @throws Exception
	 */
	public void addObservations(double[][] observations, int startTime,
			int numTimeSteps) throws Exception;
	
	/**
	 * This method will return an exception if this class was not initialised 
	 *  for univariate data only.
	 *  
	 * {@inheritDoc}
	 * 
	 * @throws Exception if initialised dimensions were not 1
	 */
	@Override
	public void setObservations(double[] observations,
			boolean[] valid) throws Exception;

	/**
	 * As per {@link #setObservations(double[], boolean[])} only with multivariate 
	 *  observations
	 * 
	 * @param observations time-series array of (multivariate) samples,
	 *  where the first index is time.
	 * @throws Exception
	 */
	public void setObservations(double[][] observations,
			boolean[] valid) throws Exception;
	
	/**
	 * This method will return an exception if this class was not initialised 
	 *  for univariate data only.
	 *  
	 * {@inheritDoc}
	 * 
	 * @throws Exception if initialised dimensions were not 1
	 */
	@Override
	public void addObservations(double[] observations, boolean[] valid) throws Exception;

	/**
	 * As per {@link #addObservations(double[], boolean[])} only with multivariate 
	 *  observations
	 * 
	 * @param observations time-series array of (multivariate) samples,
	 *  where the first index is time.
	 * @throws Exception
	 */
	public void addObservations(double[][] observations, boolean[] valid) throws Exception;

	/**
	 * This method will return an exception if this class was not initialised 
	 *  for univariate data only.
	 *  
	 * {@inheritDoc}
	 * 
	 * @throws Exception if initialised dimensions were not 1
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(double[] newObservations) throws Exception;

	/**
	 * As per {@link #computeLocalUsingPreviousObservations(double[])} only with multivariate 
	 *  observations
	 * 
	 * @param newObservations time-series array of (multivariate) samples,
	 *  where the first index is time.
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] newObservations) throws Exception;
}
