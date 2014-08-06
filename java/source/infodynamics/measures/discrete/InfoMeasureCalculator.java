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
 * <p>Info theoretic measure calculator base class, providing common functionality
 *  for user-level measure classes.</p>
 * 
 * <p>Usage of child classes is intended to follow this general pattern:
 * <ol>
 * 	<li>Construct;</li>
 *  <li>{@link #initialise()};</li>
 *  <li>Then, either of the following:
 *  	<ol>
 *  		<li>Continuous accumulation of observations before computing, via:
 * 				<ol>
 * 					<li>calling "addObservations()" methods of children
 * 					several times over;</li>
 *  				<li>Compute required quantities, using 
 *  				{@link #computeAverageLocalOfObservations()} or 
 *		        	 "computeLocalUsinPreviousObservations()".</li>
 *         		</ol></li>
 *  		<li>Standalone computation from a single set of observations; call:
 *   			"computeLocal()" or "computeAverageLocal()".</li>
 * 		</ol>
 * </ol></p>
 * 
 * <p>Note: various functions referred to above (e.g. "addObservations()")
 * are not specified here, so that this class can be a superclass
 * for both univariate, pairwise and multivariate methods.</p>
 * 
 * @author Joseph Lizier
 * joseph.lizier at gmail.com
 * http://lizier.me/joseph/
 *
 */
public abstract class InfoMeasureCalculator {

	protected double average = 0.0;
	protected double max = 0.0;
	protected double min = 0.0;
	protected double std = 0.0;
	protected int observations = 0;
	protected int base = 0; // number of individual states. Need initialised to 0 for changedSizes

	protected double log_base = 0;
	protected double log_2 = Math.log(2.0);
	protected boolean power_of_2_base = false;
	protected int log_2_base = 0;
	
	protected boolean debug = false;
	
	/**
	 * 
	 * @param blocksize
	 * @param base
	 */
	protected InfoMeasureCalculator(int base) {
		
		this.base = base;
		log_base = Math.log(base);

		if (base < 2) {
			throw new RuntimeException("Can't calculate info theoretic measures for base " + base);
		}
		
		// Check if we've got a power of 2
		power_of_2_base = isPowerOf2(base);
		if (power_of_2_base) {
			log_2_base = (int) Math.round(Math.log(base) / Math.log(2));
		}
	}
	
	/**
	 * Initialise calculator, preparing to take observation sets in
	 * Should be called prior to any of the addObservations() methods.
	 * You can reinitialise without needing to create a new object.
	 */
	public void initialise(){
		average = 0.0;
		max = 0.0;
		min = 0.0;
		std = 0.0;
		observations = 0;
	}
	
	public final double getLastAverage() {
		return average;
	}

	/**
	 * Not declaring this final so that separable calculator
	 *  can throw an exception on it since it does not support it
	 * 
	 * @return
	 */
	public double getLastMax() {
		return max;
	}

	/**
	 * Not declaring this final so that separable calculator
	 *  can throw an exception on it since it does not support it
	 * 
	 * @return
	 */
	public double getLastMin() {
		return min;
	}
	
	public final double getLastStd() {
		return std;
	}
	
	public final int getNumObservations() {
		return observations;
	}
	
	public final static boolean isPowerOf2(int num) {
		int bits = 0;
		int shiftedValue = num;
		for (int b = 0; b < Integer.SIZE; b ++) {
			if ((shiftedValue & 0x01) > 0) {
				// LSB is a 1
				bits++;
				if (bits > 1) {
					// Encountered more than 1 bit set to 1.
					// Is not a power of 2
					return false;
				}
			}
			// Shift a new bit down into the LSB
			shiftedValue = shiftedValue >> 1;
		}
		// Post: num has either 1 bit set to 1 or 0 bits set to 1
		if (bits == 0) {
			return false;
		}
		return true;
	}
	
	/**
	 * Compute the average value of the measure from the previously supplied
	 * observations
	 * 
	 * @return average value
	 */
	public abstract double computeAverageLocalOfObservations();
	
	/**
	 * @param debug the debug status to set
	 */
	public void setDebug(boolean debug) {
		this.debug = debug;
	}
}
