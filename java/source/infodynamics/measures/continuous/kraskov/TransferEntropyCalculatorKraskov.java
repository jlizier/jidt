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

package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.TransferEntropyCalculator;
import infodynamics.measures.continuous.TransferEntropyCalculatorViaCondMutualInfo;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential transfer entropy (TE) between two univariate
 *  <code>double[]</code> time-series of observations
 *  (implementing {@link TransferEntropyCalculator}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see references below).
 *  This estimator is realised here by plugging in
 *  a {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}
 *  as the calculator into the parent class {@link TransferEntropyCalculatorViaCondMutualInfo}.</p>
 *  
 * <p>Crucially, the calculation is performed by examining
 * neighbours in the full joint space (as specified by Frenzel and Pompe,
 * and Gomez-Herrero et al.)
 * rather than two MI calculators.</p>
 * 
 * <p>Usage is as per the paradigm outlined for {@link TransferEntropyCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step is either a simple call to {@link #TransferEntropyCalculatorKraskov()},
 *      or else specifies which KSG algorithm to implement via
 *      {@link #TransferEntropyCalculatorKraskov(String)};</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties defined for both
 * 		{@link TransferEntropyCalculator#setProperty(String, String)} and
 *      {@link ConditionalMutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)}
 *      as outlined
 *      in {@link TransferEntropyCalculatorViaCondMutualInfo#setProperty(String, String)});
 *      as well as for {@link #PROP_KRASKOV_ALG_NUM}.
 *      Embedding parameters may be automatically determined as outlined in 
 *      {@link TransferEntropyCalculatorViaCondMutualInfo#setProperty(String, String)}.</li>
 *      </li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 * 	<li>Frenzel and Pompe, <a href="http://dx.doi.org/10.1103/physrevlett.99.204101">
 * 	"Partial Mutual Information for Coupling Analysis of Multivariate Time Series"</a>,
 * 	Physical Review Letters, <b>99</b>, p. 204101+ (2007).</li>
 * 	<li>G. Gomez-Herrero, W. Wu, K. Rutanen, M. C. Soriano, G. Pipa, and R. Vicente,
 * 	<a href="http://arxiv.org/abs/1008.0539">
 * 	"Assessing coupling dynamics from an ensemble of time series"</a>,
 * 	arXiv:1008.0539 (2010).</li>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see TransferEntropyCalculator
 * @see ConditionalMutualInfoCalculatorMultiVariateKraskov
 */
public class TransferEntropyCalculatorKraskov
	extends TransferEntropyCalculatorViaCondMutualInfo {

	/**
	 * Class name for KSG conditional MI estimator via KSG algorithm 1
	 */
	public static final String COND_MI_CALCULATOR_KRASKOV1 = ConditionalMutualInfoCalculatorMultiVariateKraskov1.class.getName();
	/**
	 * Class name for KSG conditional MI estimator via KSG algorithm 2
	 */
	public static final String COND_MI_CALCULATOR_KRASKOV2 = ConditionalMutualInfoCalculatorMultiVariateKraskov2.class.getName();
	
	/**
	 * Property for setting which underlying Kraskov-Grassberger algorithm to use  (1 or 2).
	 * Will only be applied at the next initialisation.
	 */
	public final static String PROP_KRASKOV_ALG_NUM = "ALG_NUM";
	
	/**
	 * Which Kraskov algorithm number we are using
	 */
	protected int kraskovAlgorithmNumber = 1;
	protected boolean algChanged = false;

	/**
	 * Creates a new instance of the Kraskov-estimate style transfer entropy calculator
	 * 
	 * Uses algorithm 1 by default, as per Gomez-Herro et al.
	 * 
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public TransferEntropyCalculatorKraskov() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(COND_MI_CALCULATOR_KRASKOV1);
		kraskovAlgorithmNumber = 1;
	}

	/**
	 * Creates a new instance of the Kraskov-Grassberger style transfer entropy calculator,
	 *  with the supplied conditional MI calculator name
	 * 
	 * @param calculatorName fully qualified name of the underlying MI class.
	 *    Must be {@link #COND_MI_CALCULATOR_KRASKOV1} or {@link #COND_MI_CALCULATOR_KRASKOV2}
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public TransferEntropyCalculatorKraskov(String calculatorName) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(calculatorName);
		// Now check that it was one of our Kraskov-Grassberger calculators:
		if (calculatorName.equalsIgnoreCase(COND_MI_CALCULATOR_KRASKOV1)) {
			kraskovAlgorithmNumber = 1;
		} else if (calculatorName.equalsIgnoreCase(COND_MI_CALCULATOR_KRASKOV2)) {
			kraskovAlgorithmNumber = 2;
		} else {
			throw new ClassNotFoundException("Must be an underlying Kraskov-Grassberger conditional MI calculator");
		}
	}

	/**
	 * Creates a new instance of the Kraskov-Grassberger style transfer entropy calculator,
	 *  with the CMI Kraskov calculator using the given Kraskov algorithm number
	 * 
	 * @param calculatorName fully qualified name of the underlying MI class.
	 *    Must be {@link #COND_MI_CALCULATOR_KRASKOV1} or {@link #COND_MI_CALCULATOR_KRASKOV2}
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public TransferEntropyCalculatorKraskov(int algorithm) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super((algorithm == 1) ? COND_MI_CALCULATOR_KRASKOV1 : COND_MI_CALCULATOR_KRASKOV2);
		if ((algorithm != 1) && (algorithm != 2)) {
			throw new ClassNotFoundException("Algorithm must be 1 or 2");
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.TransferEntropyCalculatorViaCondMutualInfo#initialise(int, int, int, int, int)
	 */
	@Override
	public void initialise(int k, int k_tau, int l, int l_tau, int delay)
			throws Exception {
		if (algChanged) {
			// The algorithm number was changed in a setProperties call:
			String newCalcName = COND_MI_CALCULATOR_KRASKOV1;
			if (kraskovAlgorithmNumber == 2) {
				newCalcName = COND_MI_CALCULATOR_KRASKOV2;
			}
			@SuppressWarnings("unchecked")
			Class<ConditionalMutualInfoCalculatorMultiVariate> condMiClass = 
					(Class<ConditionalMutualInfoCalculatorMultiVariate>) Class.forName(newCalcName);
			ConditionalMutualInfoCalculatorMultiVariate newCondMiCalc = condMiClass.newInstance();
			construct(newCondMiCalc);
			// Set the properties for the Kraskov MI calculators (may pass in properties for our super class
			//  as well, but they should be ignored)
			for (String key : props.keySet()) {
				newCondMiCalc.setProperty(key, props.get(key));
			}
			algChanged = false;
		}
		
		super.initialise(k, k_tau, l, l_tau, delay);
	}

	/**
	 * Sets properties for the TE calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #PROP_KRASKOV_ALG_NUM} -- which Kraskov algorithm number to use (1 or 2).</li>
	 * 		<li>Any properties accepted by {@link TransferEntropyCalculatorViaCondMutualInfo#setProperty(String, String)}</li>
	 * 		<li>Or properties accepted by the underlying
	 * 		{@link ConditionalMutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)} implementation.</li>
	 * </ul>
	 * <p><b>Note:</b> further properties may be defined by child classes.</p>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value).
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		if (propertyName.equalsIgnoreCase(PROP_KRASKOV_ALG_NUM)) {
			int previousAlgNumber = kraskovAlgorithmNumber;
			kraskovAlgorithmNumber = Integer.parseInt(propertyValue);
			if ((kraskovAlgorithmNumber != 1) && (kraskovAlgorithmNumber != 2)) {
				throw new Exception("Kraskov algorithm number (" + kraskovAlgorithmNumber
						+ ") must be either 1 or 2");
			}
			if (kraskovAlgorithmNumber != previousAlgNumber) {
				algChanged = true;
			}
			if (debug) {
				System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
						" to " + propertyValue);
			}
		} else {
			// Assume it was a property for the parent class or underlying conditional MI calculator
			super.setProperty(propertyName, propertyValue);
		}
	}

	@Override
	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(PROP_KRASKOV_ALG_NUM)) {
			return Integer.toString(kraskovAlgorithmNumber);
		} else if (propertyName.equalsIgnoreCase(PROP_RAGWITZ_NUM_NNS)) {
			// Need to deal with this one here instead of in super class, 
			//  since if it's not set especially by the user we'll default to the
			//  kNNs used in the KSG method:
			if (ragwitz_num_nns_set) {
				return Integer.toString(ragwitz_num_nns);
			} else {
				return condMiCalc.getProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K);
			}
		} else {
			// Assume it was a property for the parent class or underlying conditional MI calculator
			return super.getProperty(propertyName);
		}
	}
	
	/**
     * Debug method to return the k nearest neighbour distances that 
     *  would be utilised for each sample point here.
     * Note that this is specifically the max-norm across the source-target-targetPast variables, which is used for each
     *  range search in algorithm 1 (although algorithm 2 would use the max distance
     *  for each variable within the kNNs in their separate range searches). 
     *  
     * @param startTimePoint
     * @param numTimePoints
     * @return
     * @throws Exception
     */
	public double[] kNNDistances(int startTimePoint, int numTimePoints) throws Exception {
		// Defer the call to the underlying KSG CMI estimator
		return ((ConditionalMutualInfoCalculatorMultiVariateKraskov) condMiCalc).kNNDistances(startTimePoint, numTimePoints);
	}
	
	/**
	 * Debug method to return the k nearest neighbour distances that
	 *  would be utilised in {@link #computeLocalUsingPreviousObservations(double[], double[])}
	 *  for a cross TE.
     * Note that this is specifically the max-norm across the source-target-targetPast variables, which is used for each
     *  range search in algorithm 1 (although algorithm 2 would use the max distance
     *  for each variable within the kNNs in their separate range searches). 
	 *  
	 * @param startTimePoint
	 * @param numTimePoints
	 * @param newSourceObservations
	 * @param newDestObservations
	 * @return
	 * @throws Exception
	 */
	public double[] kNNDistancesForNewSamples(int startTimePoint, int numTimePoints,
			double[] newSourceObservations, double[] newDestObservations) throws Exception {
		if (newSourceObservations.length != newDestObservations.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					newSourceObservations.length, newDestObservations.length));
		}
		if (newDestObservations.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to compute for here
			return new double[newDestObservations.length];
		}
		// Now embed as per computeLocalUsingPreviousObservations() in super:
		double[][] newDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newDestObservations, k, k_tau,
						startTimeForFirstDestEmbedding,
						newDestObservations.length - startTimeForFirstDestEmbedding - 1);
		double[][] newDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(newDestObservations, 1,
						startTimeForFirstDestEmbedding + 1,
						newDestObservations.length - startTimeForFirstDestEmbedding - 1);
		double[][] newSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newSourceObservations, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						newSourceObservations.length - startTimeForFirstDestEmbedding - 1);
		return ((ConditionalMutualInfoCalculatorMultiVariateKraskov) condMiCalc).
				kNNDistancesForNewSamples(startTimePoint, numTimePoints, newSourcePastVectors, newDestNextVectors, newDestPastVectors);
	}

}
