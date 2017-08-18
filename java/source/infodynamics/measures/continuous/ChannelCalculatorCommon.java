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
 * A basic interface for calculators computing measures on a <i>channel</i> from a
 * source to a destination time-series (i.e. mutual information and transfer entropy).
 * In the following, we refer to the abstract measure computed by this calculator
 * as the <i>"channel measure"</i>.
 * 
 * <p>The interface is <b>abstract</b> because further specification is required
 * for addObservations or setObservations methods showing whether
 * this is a univariate ({@link ChannelCalculator}) or multivariate 
 * {@link ChannelCalculatorMultiVariate} measure.</p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator using {@link #initialise()} or
 *			other initialise methods defined by child classes;
 *		</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>a "setObservations" method (defined by child classes)
 * 				for calculations based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to "addObservations" methods
 * 						(defined by children), then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul>
 * 		</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average of the channel measure: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local channel measure values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local channel measure values for a specific set of samples: "computeLocalUsingPreviousObservations"
 * 					(to be defined by child classes);</li>
 * 				<li>the distribution of the channel measure values under the null hypothesis
 * 					of no relationship between source and destination
 * 					values: {@link #computeSignificance(int)} or
 * 					{@link #computeSignificance(int[][])}.</li>
 * 			</ul>
 * 		</li>
 * 		<li>
 * 		Return to step 2 or 3 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see ChannelCalculator
 * @see ChannelCalculatorMultiVariate
 *
 */
public abstract interface ChannelCalculatorCommon
	extends InfoMeasureCalculatorContinuous, EmpiricalNullDistributionComputer {

	/**
	 * Signal that we will add in the samples for computing the PDF 
	 * from several disjoint time-series or trials via calls to
	 * "addObservations" rather than "setObservations" type methods
	 * (defined by the child interfaces and classes).
	 *
	 */
	public void startAddObservations();
	
	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
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
	 * <p>Computes the local values of the implemented channel measure,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>If the samples were supplied via a single call such as
	 * {@link #setObservations(double[])},
	 * then the return value is a single time-series of local
	 * channel measure values corresponding to these samples.</p>
	 * 
	 * <p>Otherwise where disjoint time-series observations were supplied using several 
	 *  calls such as {@link addObservations(double[])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single "time-series" return array.</p>
	 *  
	 * @return the "time-series" of local channel measure values.
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception;

}
