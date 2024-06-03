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

package infodynamics.measures.continuous.kernel;

import infodynamics.measures.continuous.EntropyCalculatorMultiVariate;
import infodynamics.measures.continuous.EntropyCalculatorMultiVariateCommon;

/**
 * <p>Computes the differential entropy of a given set of observations
 *  (implementing {@link EntropyCalculatorMultiVariate}, using box-kernel estimation.
 *  For details on box-kernel estimation, see Kantz and Schreiber (below).</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link EntropyCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #EntropyCalculatorMultiVariateKernel()}.</li>
 *  <li>Further properties are available, see {@link #setProperty(String, String)};</li>
 *  <li>An additional {@link #initialise(int, double)} option;</li>
 * </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>H. Kantz and T. Schreiber, "Nonlinear Time Series Analysis".
 *   Cambridge, MA: Cambridge University Press, 1997.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EntropyCalculatorMultiVariateKernel
	extends EntropyCalculatorMultiVariateCommon 
	implements EntropyCalculatorMultiVariate
	{

	private KernelEstimatorMultiVariate mvke = null;
	
	/**
	 * Default value for kernel width
	 */
	private static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	private double kernelWidth = DEFAULT_EPSILON;
	/**
	 * Property name for the kernel width
	 */
	public static final String KERNEL_WIDTH_PROP_NAME = "KERNEL_WIDTH";
	/**
	 * Legacy property name for the kernel width
	 */
	public static final String EPSILON_PROP_NAME = "EPSILON";

	/**
	 * Construct an instance
	 */
	public EntropyCalculatorMultiVariateKernel() {
		super();
		mvke = new KernelEstimatorMultiVariate();
		mvke.setDebug(debug);
		mvke.setNormalise(normalise);
	}

	public void initialise(int dimensions) {
		initialise(dimensions, kernelWidth);
	}

	/**
	 * Initialise the calculator for (re-)use, with a specific 
	 * number of joint variables, specific kernel width,
	 * and existing (or default) values of other parameters.
	 * Clears an PDFs of previously supplied observations.
	 *
	 * @param dimensions number of joint variables
	 * @param kernelWidth if {@link #NORMALISE_PROP_NAME} property has
	 *  been set, then this kernel width corresponds to the number of
	 *  standard deviations from the mean (otherwise it is an absolute value)
	 */
	public void initialise(int dimensions, double kernelWidth) {
		super.initialise(dimensions);
		this.kernelWidth = kernelWidth;
		mvke.initialise(dimensions, kernelWidth);
	}

	@Override
	public void finaliseAddObservations() throws Exception {
		super.finaliseAddObservations();
		mvke.setObservations(observations);
	}

	@Override
	public double computeAverageLocalOfObservations() {
		if (isComputed) {
			return lastAverage;
		}
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvke.getProbability(observations[b]);
			double cont = Math.log(prob);
			entropy -= cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + (-cont/Math.log(2.0)) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		lastAverage = entropy / (double) totalObservations / Math.log(2.0);
		isComputed = true;
		return lastAverage;
	}

	@Override
	public double[] computeLocalOfPreviousObservations() {
		return computeLocalUsingPreviousObservations(true, observations);
	}

	@Override
	public double[] computeLocalUsingPreviousObservations(double states[][]) {
		return computeLocalUsingPreviousObservations(false, states);
	}
	
	/**
	 * Internal method for computing locals of a set of observations
	 *  based on existing PDFs.
	 * 
	 * @param isPreviousObservations whether we're using the points
	 *   that the PDFs were computed from or not.
	 * @param states samples to compute the local joint entropies on.
	 * @return
	 */
	protected double[] computeLocalUsingPreviousObservations(
			boolean isPreviousObservations, double states[][]) {
		
		double entropy = 0.0;
		double[] localEntropy = new double[states.length];
		for (int b = 0; b < states.length; b++) {
			double prob = mvke.getProbability(states[b]);
			double cont = -Math.log(prob);
			localEntropy[b] = cont / Math.log(2.0);
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		entropy = entropy / (double) totalObservations / Math.log(2.0); // Don't use /= as I'm not sure of the order of operations it applies.
		if (isPreviousObservations) {
			lastAverage = entropy;
			isComputed = true;
		}
		return localEntropy;
	}

	@Override
	public void setDebug(boolean debug) {
		super.setDebug(debug);
		mvke.setDebug(debug);
	}

	/**
	 * <p>Set properties for the kernel entropy calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #KERNEL_WIDTH_PROP_NAME} (legacy value is {@link #EPSILON_PROP_NAME}) --
	 * 			kernel width to be used in the calculation. If {@link #normalise} is set,
	 * 		    then this is a number of standard deviations; otherwise it
	 * 			is an absolute value. Default is {@link #DEFAULT_KERNEL_WIDTH}.</li>
	 *  	<li>any valid properties for {@link EntropyCalculatorMultiVariate#setProperty(String, String)}.</li>
	 * </ul> 
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;

		// TODO If we implement a dynamic correlation exclusion property,
		//  then we will need to call getProbability(double, int) instead of
		//  just getProbability(double) above.
		
		if (propertyName.equalsIgnoreCase(KERNEL_WIDTH_PROP_NAME) ||
				propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			kernelWidth = Double.parseDouble(propertyValue);
		} else {
			// No property was set
			propertySet = false;
			// try the superclass:
			super.setProperty(propertyName, propertyValue);
			if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
				// Handle an additional step for this one:
				mvke.setNormalise(normalise);
			}
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(KERNEL_WIDTH_PROP_NAME) ||
				propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			return Double.toString(kernelWidth);
		} else {
		      // try the superclass:
		      return super.getProperty(propertyName);
		}
	}

}
