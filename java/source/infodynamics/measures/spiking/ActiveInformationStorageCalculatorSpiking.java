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
 * Interface for implementations of Active Information Storage (AIS) calculators
 * on spike trains or other event-based data.
 * 
 * <p>AIS is the mutual information between a process' past and its next state.
 * For spiking processes, this measures how much information about future spiking
 * is contained in the past spiking patterns.</p>
 * 
 * <p>Usage is intended to follow this paradigm:</p>
 * <ol>
 *  <li>Construct the calculator via an implementing class constructor;</li>
 *  <li>Set properties using {@link #setProperty(String, String)};</li>
 *  <li>Initialize the calculator using {@link #initialise()} or
 *      {@link #initialise(int)};</li>
 *  <li>Provide the observations/samples for the calculator to set up the PDFs,
 *      using one or more calls to sets of {@link #setObservations(double[])} or
 *      {@link #addObservations(double[])} or
 *      {@link #startAddObservations()}, followed by multiple calls to
 *      {@link #addObservations(double[])} and then
 *      {@link #finaliseAddObservations()};</li>
 *  <li>Compute the required quantities, being one or more of:
 *      <ul>
 *          <li>the average AIS: {@link #computeAverageLocalOfObservations()};</li>
 *          <li>the local AIS values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 *
 *          <li>the distribution of AIS values under the null hypothesis
 *              of no relationship between past and future:
 *              {@link #computeSignificance(int)} or
 *              {@link #computeSignificance(int, double)}.</li>
 *      </ul>
 *  </li>
 * </ol>
 * 
 * @author Michael Fang (fangmichael33@gmail.com)
 *
 */
public interface ActiveInformationStorageCalculatorSpiking {



    /**
	 * Property name for a comma-separated list of integers representing 
	 * the past intervals to use in the embedding.
	 */
	public static final String PAST_INTERVALS_PROP_NAME = "PAST_INTERVALS";



    /**
     * Initialise the calculator for (re-)use, with the existing
     * (or default) values of parameters.
     * 
     * @throws Exception
     */
    public void initialise() throws Exception;
    
    /**
     * Initialise the calculator for (re-)use, with some parameters
     * supplied here rather than in later method calls.
     * 
     * @param k Length of past history to consider (i.e. embedding length)
     * @throws Exception
     */
    public void initialise(int k) throws Exception;
    
    /**
     * Set properties for the calculator.
     * New property values are not guaranteed to take effect until the next call
     * to an initialise method. 
     * 
     * <p>Valid property names, and what their
     * values should represent, include:</p>
     * <ul>
     *  <li>{@link #PAST_INTERVALS_PROP_NAME} -- a comma-separated list
     *      of integers representing the past intervals to use in the embedding.</li>
     *  <li>Any other properties defined by the implementing class.</li>
     * </ul>
     * 
     * @param propertyName name of the property to set
     * @param propertyValue value of the property to set
     * @throws Exception for invalid property values
     */
    public void setProperty(String propertyName, String propertyValue) throws Exception;
    
    /**
     * Get property values for the calculator.
     * 
     * <p>Valid property names, and what their
     * values should represent, are the same as those for
     * {@link #setProperty(String, String)}</p>
     * 
     * @param propertyName name of the property
     * @return the value of the property
     * @throws Exception for invalid property values
     */
    public String getProperty(String propertyName) throws Exception;
    
    /**
     * Sets a single set of observations for the calculator to use.
     * Cannot be called once {@link #startAddObservations()} has been called,
     * and cannot be called after {@link #addObservations(double[])} has been
     * called.
     * 
     * @param observations time-series array of spike times
     * @throws Exception
     */
    public void setObservations(double[] observations) throws Exception;
    
    /**
     * Signal that we will add in the observations for calculating the PDFs
     * from several disjoint time-series or trials.
     * 
     * @throws Exception
     */
    public void startAddObservations() throws Exception;
    
    /**
     * Add observations for the PDFs for a single time-series.
     * 
     * @param observations time-series array of spike times
     * @throws Exception
     */
    public void addObservations(double[] observations) throws Exception;
    
    /**
     * Signal that we have finished adding in the observations.
     * 
     * @throws Exception
     */
    public void finaliseAddObservations() throws Exception;
    
    /**
     * Returns whether more than one time-series has been added
     * to the calculator (either via {@link #setObservations(double[])}
     * or via {@link #addObservations(double[])})
     * 
     * @return true if more than one time-series has been supplied
     */
    public boolean getAddedMoreThanOneObservationSet();
    
    /**
     * Compute the average AIS from the previously-supplied samples.
     * 
     * @return the estimate of the AIS
     * @throws Exception
     */
    public double computeAverageLocalOfObservations() throws Exception;
    
    /**
     * This interface serves to indicate the return type of {@link #computeLocalOfPreviousObservations()}
     * as each child implementation will return something specific
     */
    public interface SpikingLocalInformationValues {
        // Left empty intentionally
    }
    
    /**
     * Compute the local AIS values for the previously-supplied samples.
     * 
     * @return an object containing a representation of the local AIS values
     * @throws Exception
     */
    public SpikingLocalInformationValues computeLocalOfPreviousObservations() throws Exception;
    
    public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck, double estimatedValue) throws Exception;
    
    public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck, double estimatedValue, long randomSeed) throws Exception;
    
    /**
     * Set whether to display debug messages or not.
     * 
     * @param debug display debug messages if true
     */
    public void setDebug(boolean debug);
    
    /**
     * Return the AIS last calculated in a call to {@link #computeAverageLocalOfObservations()}
     * or {@link #computeLocalOfPreviousObservations()} after the previous
     * {@link #initialise()} call.
     * 
     * @return the last computed average AIS value
     */
    public double getLastAverage();
} 