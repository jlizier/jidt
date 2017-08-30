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

package infodynamics.utils;

/**
 * Calculators implementing this interface must provide
 *  {@link #computeSignificance(int)} and {@link #computeSignificance(int[][])}
 *  methods to compute
 *  the statistical significance of their measurement, returning an empirically
 *  determined distribution {@link EmpiricalMeasurementDistribution}.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface EmpiricalNullDistributionComputer {

	/**
	 * Generate an <b>empirical</b> (bootstrapped) distribution of what the given measure would look like,
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination value (possibly
	 * in the context of a conditional value).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * conditional MI and TE.
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using addObservations methods etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p><b>References:</b><br/>
	 *  <ul>
	 *   <li>J.T. Lizier, "JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems", 2014.</li>
     * </ul>
	 * 
	 * @param numPermutationsToCheck number of surrogate samples to bootstrap
	 *  to generate the distribution.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @return the empirical distribution of measure scores under this null hypothesis.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception;
	
	/**
	 * Generate an <b>empirical</b> (bootstrapped) distribution of what the given measure would look like,
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination value (possibly
	 * in the context of a conditional value).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done MI,
	 * conditional MI and TE.
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using addObservations methods etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p>This method (in contrast to {@link #computeSignificance(int)})
	 * allows the user to specify how to construct the surrogates,
	 * such that repeatable results may be obtained.</p>
	 * 
	 * <p><b>References:</b><br/>
	 *  <ul>
	 *   <li>J.T. Lizier, "JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems", 2014.</li>
     * </ul>
	 * 
	 * @param newOrderings a specification of how to shuffle the next values
	 *  to create the surrogates to generate the distribution with. The first
	 *  index is the permutation number (i.e. newOrderings.length is the number
	 *  of surrogate samples we use to bootstrap to generate the distribution here.)
	 *  Each array newOrderings[i] should be an array of length N (where
	 *  would be the value returned by a call to getNumObservations() for the given measure)
	 *  containing a permutation of the values in 0..(N-1).
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @return the empirical distribution of measure scores under this null hypothesis.
	 * @throws Exception where e.g. the newOrderings don't supply arrays of the correct
	 *   length matching the number of observations that we have.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception;
	
}
