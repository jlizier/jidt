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

import java.util.Hashtable;

import infodynamics.measures.continuous.ActiveInfoStorageCalculatorMultiVariate;
import infodynamics.measures.continuous.ActiveInfoStorageCalculatorMultiVariateViaMutualInfo;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.MatrixUtils;

/**
 * An Active Information Storage (AIS) calculator for multivariate time-series
 * data (implementing {@link ActiveInfoStorageCalculatorMultiVariate})
 * which is affected using a 
 * Kraskov-Stoegbauer-Grassberger (KSG) Mutual Information (MI) calculator
 * ({@link MutualInfoCalculatorMultiVariateKraskov}) to make the calculations.
 * 
 * <p>
 * That is, this class implements an AIS calculator using the KSG nearest-neighbour approach.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateKraskov}
 * as the calculator into the parent class {@link ActiveInfoStorageCalculatorMultiVariateViaMutualInfo}.
 * </p> 
 * 
 * <p>Usage is as per the paradigm outlined for {@link ActiveInfoStorageCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step is either a simple call to {@link #ActiveInfoStorageCalculatorKraskov()},
 *      or else specifies which KSG algorithm to implement via {@link #ActiveInfoStorageCalculatorMultiVariateKraskov(int)}
 *      or {@link #ActiveInfoStorageCalculatorMultiVariateKraskov(String)};</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties for
 *      {@link MutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)}
 *      (except {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} as outlined
 *      in {@link ActiveInfoStorageCalculatorMultiVariateViaMutualInfo#setProperty(String, String)}).
 *      Embedding parameters may be automatically determined as per the Ragwitz criteria
 *      by setting the property {@link #PROP_AUTO_EMBED_METHOD} to {@link #AUTO_EMBED_METHOD_RAGWITZ},
 *      or as per the max. bias-corrected AIS criteria by 
 *      setting the property {@link #PROP_AUTO_EMBED_METHOD} to {@link #AUTO_EMBED_METHOD_MAX_CORR_AIS}
 *      (plus additional parameter settings for these).</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>J.T. Lizier, M. Prokopenko and A.Y. Zomaya,
 * 		<a href="http://dx.doi.org/10.1016/j.ins.2012.04.016">
 * 		"Local measures of information storage in complex distributed computation"</a>,
 * 		Information Sciences, vol. 208, pp. 39-54, 2012.</li>
 * 	<li>Ragwitz and Kantz, "Markov models from data by simple nonlinear time series
 *  	predictors in delay embedding spaces", Physical Review E, vol 65, 056201 (2002).</li>
 *  <li>J. Garland, R. G. James, E. Bradley, <a href="http://dx.doi.org/10.1103/physreve.93.022221">
 *  	"Leveraging information storage to select forecast-optimal parameters for delay-coordinate reconstructions"</a>,
 *  	Physical Review E, Vol. 93 (2016), 022221, doi:</li>
 * </ul>
 * 
 * @author Pedro AM Mediano (<a href="pmediano at imperial.ac.uk">email</a>,
 * <a href="https://www.doc.ic.ac.uk/~pam213/">www</a>)
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 *
 * @see ActiveInfoStorageCalculatorMultiVariate
 * @see ActiveInfoStorageCalculatorMultiVariateViaMutualInfo
 * @see MutualInfoCalculatorMultiVariateKraskov
 *
 */
public class ActiveInfoStorageCalculatorMultiVariateKraskov
	extends ActiveInfoStorageCalculatorMultiVariateViaMutualInfo {
	
	/**
	 * Class name for KSG MI estimator via KSG algorithm 1
	 */
	public static final String MI_CALCULATOR_KRASKOV1 = MutualInfoCalculatorMultiVariateKraskov1.class.getName();
	/**
	 * Class name for KSG MI estimator via KSG algorithm 2
	 */
	public static final String MI_CALCULATOR_KRASKOV2 = MutualInfoCalculatorMultiVariateKraskov2.class.getName();
	/**
	 * Property for setting which underlying Kraskov-Grassberger algorithm to use.
	 * Will only be applied at the next initialisation.
	 */
	public final static String PROP_KRASKOV_ALG_NUM = "ALG_NUM";
	
	/**
	 * Which Kraskov algorithm number we are using
	 */
	protected int kraskovAlgorithmNumber = 2;
	protected boolean algChanged = false;
	/**
	 * Storage for the properties ready to pass onto the underlying conditional MI calculators should they change 
	 */
	protected Hashtable<String,String> props;

	/**
	 * Property name for the auto-embedding method. Defaults to {@link #AUTO_EMBED_METHOD_NONE}
	 * Other valid values are {@link #AUTO_EMBED_METHOD_RAGWITZ} or
	 * {@link #AUTO_EMBED_METHOD_MAX_CORR_AIS}
	 */
	public static final String PROP_AUTO_EMBED_METHOD = "AUTO_EMBED_METHOD";
	/**
	 * Valid value for the property {@link #PROP_AUTO_EMBED_METHOD} indicating that
	 *  no auto embedding should be done (i.e. to use manually supplied parameters)
	 */
	public static final String AUTO_EMBED_METHOD_NONE = "NONE";
	/**
	 * Valid value for the property {@link #PROP_AUTO_EMBED_METHOD} indicating that
	 *  the Ragwitz optimisation technique should be used for automatic embedding
	 */
	public static final String AUTO_EMBED_METHOD_RAGWITZ = "RAGWITZ";
	/**
	 * Valid value for the property {@link #PROP_AUTO_EMBED_METHOD} indicating that
	 *  the automatic embedding should be done by maximising the bias corrected
	 *  AIS (as per Garland et al. in the references above).
	 */
	public static final String AUTO_EMBED_METHOD_MAX_CORR_AIS = "MAX_CORR_AIS";
	/**
	 * Internal variable tracking what type of auto embedding (if any)
	 *  we are using
	 */
	protected String autoEmbeddingMethod = AUTO_EMBED_METHOD_NONE;
	
	/**
	 * Property name for maximum k (embedding length) for the auto-embedding search. Default to 1
	 */
	public static final String PROP_K_SEARCH_MAX = "AUTO_EMBED_K_SEARCH_MAX";
	/**
	 * Internal variable for storing the maximum embedding length to search up to for
	 *  automating the parameters.
	 */
	protected int k_search_max = 1;

	/**
	 * Property name for maximum tau (embedding delay) for the auto-embedding search. Default to 1
	 */
	public static final String PROP_TAU_SEARCH_MAX = "AUTO_EMBED_TAU_SEARCH_MAX";
	/**
	 * Internal variable for storing the maximum embedding delay to search up to for
	 *  automating the parameters.
	 */
	protected int tau_search_max = 1;

	/**
	 * Property name for the number of nearest neighbours to use for the auto-embedding search (Ragwitz criteria).
	 * Defaults to match the value in use for {@link MutualInfoCalculatorMultiVariateKraskov#PROP_K}
	 */
	public static final String PROP_RAGWITZ_NUM_NNS = "AUTO_EMBED_RAGWITZ_NUM_NNS";
	/**
	 * Internal variable for storing the number of nearest neighbours to use for the
	 *  auto embedding search (Ragwitz criteria)
	 */
	protected int ragwitz_num_nns = 1;
	/** 
	 * Internal variable to track whether the property {@link #PROP_RAGWITZ_NUM_NNS} has been
	 * set yet
	 */
	protected boolean ragwitz_num_nns_set = false;
	
	/**
	 * Creates a new instance of the Kraskov-Stoegbauer-Grassberger style AIS calculator.
	 * Uses algorithm 2 by default.
	 * 
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorMultiVariateKraskov() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(MI_CALCULATOR_KRASKOV2);
		props = new Hashtable<String,String>();
	}

	/**
	 * Creates a new instance of the Kraskov-Stoegbauer-Grassberger estimator for AIS,
	 *  with the supplied MI calculator name.
	 * 
	 * @param calculatorName fully qualified name of the underlying MI class.
	 *    Must be equal to {@link #MI_CALCULATOR_KRASKOV1} or {@link #MI_CALCULATOR_KRASKOV2}
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorMultiVariateKraskov(String calculatorName) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(calculatorName);
		// Now check that it was one of our Kraskov-Grassberger calculators:
		if (!calculatorName.equalsIgnoreCase(MI_CALCULATOR_KRASKOV1) &&
				!calculatorName.equalsIgnoreCase(MI_CALCULATOR_KRASKOV2)) {
			throw new ClassNotFoundException("Must be an underlying Kraskov-Stoegbauer-Grassberger calculator");
		}
		props = new Hashtable<String,String>();
	}

	/**
	 * Creates a new instance of the Kraskov-Stoegbauer-Grassberger estimator for AIS,
	 *  with the supplied Kraskov-Stoegbauer-Grassberger MI algorithm number
	 * 
	 * @param algorithm must be either 1 or 2 for the first or second KSG algorithm
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorMultiVariateKraskov(int algorithm) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(algorithm == 1 ? MI_CALCULATOR_KRASKOV1 : MI_CALCULATOR_KRASKOV2);
		if ((algorithm != 1) && (algorithm != 2)) {
			throw new ClassNotFoundException("Algorithm must be 1 or 2");
		}
		props = new Hashtable<String,String>();
	}

  /**
   * {@inheritDoc}
   */
  @Override
  public void initialise(int dimensions, int k, int tau) throws Exception {
		if (algChanged) {
			// The algorithm number was changed in a setProperties call:
			String newCalcName = MI_CALCULATOR_KRASKOV1;
			if (kraskovAlgorithmNumber == 2) {
				newCalcName = MI_CALCULATOR_KRASKOV2;
			}
			@SuppressWarnings("unchecked")
			Class<MutualInfoCalculatorMultiVariate> miClass = 
					(Class<MutualInfoCalculatorMultiVariate>) Class.forName(newCalcName);
			MutualInfoCalculatorMultiVariate newMiCalc = miClass.newInstance();
			construct(newMiCalc);
			// Set the properties for the Kraskov MI calculators (may pass in properties for our super class
			//  as well, but they should be ignored)
			for (String key : props.keySet()) {
				newMiCalc.setProperty(key, props.get(key));
			}
			algChanged = false;
		}
		
		super.initialise(dimensions, k, tau);
  }


	/**
	 * Sets properties for the AIS calculator.
   *
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #PROP_AUTO_EMBED_METHOD} -- method by which the calculator
	 * 		automatically determines the embedding history length ({@link #K_PROP_NAME})
	 * 		and embedding delay ({@link #TAU_PROP_NAME}). Default is {@link #AUTO_EMBED_METHOD_NONE} meaning
	 * 		values are set manually; other accepted values include: {@link #AUTO_EMBED_METHOD_RAGWITZ} for use
	 * 		of the Ragwitz criteria (searching up to {@link #PROP_K_SEARCH_MAX} and 
	 * 		{@link #PROP_TAU_SEARCH_MAX}). Use of any value other than {@link #AUTO_EMBED_METHOD_NONE}
	 * 		will lead to any previous settings for k and tau (via e.g. {@link #initialise(int, int)} or
	 * 		auto-embedding during previous calculations) will be overwritten after observations
	 * 		are supplied.</li>
	 * 		<li>{@link #PROP_K_SEARCH_MAX} -- maximum embedded history length to search
	 * 		up to if automatically determining the embedding parameters (as set by
	 * 		{@link #PROP_AUTO_EMBED_METHOD}); default is 1</li>
	 * 		<li>{@link #PROP_TAU_SEARCH_MAX} -- maximum embedded history length to search
	 * 		up to if automatically determining the embedding parameters (as set by
	 * 		{@link #PROP_AUTO_EMBED_METHOD}); default is 1</li>
	 * 		<li>{@link #PROP_RAGWITZ_NUM_NNS} -- number of nearest neighbours to use
	 * 		in the auto-embedding if the property {@link #PROP_AUTO_EMBED_METHOD}
	 * 		has been set to {@link #AUTO_EMBED_METHOD_RAGWITZ}. Defaults to the property value
	 *      set for {@link MutualInfoCalculatorMultiVariateKraskov.PROP_K}</li>
	 * 		<li>{@link #PROP_KRASKOV_ALG_NUM} -- which KSG algorithm to use, 1 or 2.
	 * 		<li>Any properties accepted by {@link super#setProperty(String, String)}</li>
	 * 		<li>Or properties accepted by the underlying
	 * 		{@link MutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)} implementation.</li>
	 * </ul>
	 * <p>One should set {@link MutualInfoCalculatorMultiVariateKraskov#PROP_K} here, the number
	 *  of neighbouring points one should count up to in determining the joint kernel size.</p> 
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value).
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_AUTO_EMBED_METHOD)) {
			// New method set for determining the embedding parameters
			autoEmbeddingMethod = propertyValue;
		} else if (propertyName.equalsIgnoreCase(PROP_K_SEARCH_MAX)) {
			// Set max embedding history length for auto determination of embedding
			k_search_max = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_TAU_SEARCH_MAX)) {
			// Set maximum embedding delay for auto determination of embedding
			tau_search_max = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_RAGWITZ_NUM_NNS)) {
			// Set the number of nearest neighbours to use in case of Ragwitz auto embedding:
			ragwitz_num_nns = Integer.parseInt(propertyValue);
			ragwitz_num_nns_set = true;
		} else if (propertyName.equalsIgnoreCase(PROP_KRASKOV_ALG_NUM)) {
			int previousAlgNumber = kraskovAlgorithmNumber;
			kraskovAlgorithmNumber = Integer.parseInt(propertyValue);
			if ((kraskovAlgorithmNumber != 1) && (kraskovAlgorithmNumber != 2)) {
				throw new Exception("Kraskov algorithm number (" + kraskovAlgorithmNumber
						+ ") must be either 1 or 2");
			}
			if (kraskovAlgorithmNumber != previousAlgNumber) {
				algChanged = true;
			}
		} else {
      propertySet = false;
			// Assume it was a property for the parent class or underlying MI calculator
			super.setProperty(propertyName, propertyValue);
			props.put(propertyName, propertyValue); // This will keep properties for the super class as well as the MI calculator, but this is ok
		}
    if (debug && propertySet) {
      System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
          " to " + propertyValue);
    }
	}

  // TODO: PEDRO: should I implement a getProperty here?
	// Joe: Yes. But look into this in wider context of what else we could cut from this class --
	//  there seems to be an awful lot that's re-implementing variables and code that we 
	//  inherit already ...
	@Override
	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(PROP_KRASKOV_ALG_NUM)) {
			return Integer.toString(kraskovAlgorithmNumber);
		} else {
			// Assume it was a property for the parent class or underlying conditional MI calculator
			return super.getProperty(propertyName);
		}
	}

	@Override
	public void preFinaliseAddObservations() throws Exception {
		// Automatically determine the embedding parameters for the given time series
		
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_NONE)) {
			return;
		}
		// Else we need to auto embed
		
		double bestPredictionError = Double.POSITIVE_INFINITY;
		int k_candidate_best = 1;
		int tau_candidate_best = 1;
		
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ)) {
			if (debug) {
				System.out.printf("Beginning Ragwitz auto-embedding with k_max=%d, tau_max=%d\n",
						k_search_max, tau_search_max);
			}
			
			for (int k_candidate = 1; k_candidate <= k_search_max; k_candidate++) {
				for (int tau_candidate = 1; tau_candidate <= tau_search_max; tau_candidate++) {
					// Use our internal MI calculator in case it has any particular 
					//  properties we need to have been set already
					miCalc.initialise(k_candidate*dimensions, dimensions);
					miCalc.startAddObservations();
					// Send all of the observations through:
          if (dimensions == 1) {

            for (double[] observations : vectorOfObservationTimeSeries) {
              double[][] currentDestPastVectors = 
                  MatrixUtils.makeDelayEmbeddingVector(observations, k_candidate,
                      tau_candidate, (k_candidate-1)*tau_candidate,
                      observations.length - (k_candidate-1)*tau_candidate - 1);
              double[][] currentDestNextVectors =
                  MatrixUtils.makeDelayEmbeddingVector(observations, 1,
                      (k_candidate-1)*tau_candidate + 1,
                      observations.length - (k_candidate-1)*tau_candidate - 1);
              miCalc.addObservations(currentDestPastVectors, currentDestNextVectors);
            }

          } else {

            for (double[][] observations : vectorOfMultiVariateObservationTimeSeries) {
              double[][] currentDestPastVectors = 
                  MatrixUtils.makeDelayEmbeddingVector(observations, k_candidate,
                      tau_candidate, (k_candidate-1)*tau_candidate,
                      observations.length - (k_candidate-1)*tau_candidate - 1);
              double[][] currentDestNextVectors =
                  MatrixUtils.makeDelayEmbeddingVector(observations, 1,
                      (k_candidate-1)*tau_candidate + 1,
                      observations.length - (k_candidate-1)*tau_candidate - 1);
              miCalc.addObservations(currentDestPastVectors, currentDestNextVectors);
            }

          }
					miCalc.finaliseAddObservations();
					// Now grab the prediction errors of the next value from the required number of
					// nearest neighbours of the previous state. Since the prediction error is a sum
          // of squared terms, we can also sum the prediction error across dimensions.
					double[] predictionError;
					if (ragwitz_num_nns_set) {
						predictionError = 
							((MutualInfoCalculatorMultiVariateKraskov) miCalc).
								computePredictionErrorsFromObservations(false, ragwitz_num_nns);
					} else {
						predictionError = 
								((MutualInfoCalculatorMultiVariateKraskov) miCalc).
									computePredictionErrorsFromObservations(false);
					}
					if (debug) {
						System.out.printf("Embedding prediction error (dim=%d) for k=%d,tau=%d is %.3f\n",
								predictionError.length, k_candidate, tau_candidate,
								MatrixUtils.sum(predictionError) / (double) miCalc.getNumObservations());
					}
					if ((MatrixUtils.sum(predictionError) / (double) miCalc.getNumObservations())
							< bestPredictionError) {
						// This parameter setting is the best so far:
						// (Note division by number of observations to normalise
						//  for less observations for larger k and tau)
						bestPredictionError = MatrixUtils.sum(predictionError) / (double) miCalc.getNumObservations();
						k_candidate_best = k_candidate;
						tau_candidate_best = tau_candidate;
					}
					if (k_candidate == 1) {
						// tau is irrelevant, so no point testing other values
						break;
					}
				}
			}
		} else if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS)) {
			double bestAIS = Double.NEGATIVE_INFINITY;
			if (debug) {
				System.out.printf("Beginning max bias corrected AIS auto-embedding with k_max=%d, tau_max=%d\n",
						k_search_max, tau_search_max);
			}
			
			for (int k_candidate = 1; k_candidate <= k_search_max; k_candidate++) {
				for (int tau_candidate = 1; tau_candidate <= tau_search_max; tau_candidate++) {
					// Use our internal MI calculator in case it has any particular 
					//  properties we need to have been set already
					miCalc.initialise(k_candidate*dimensions, dimensions);
					miCalc.startAddObservations();
					// Send all of the observations through:
          if (dimensions == 1) {

            for (double[] observations : vectorOfObservationTimeSeries) {
              double[][] currentDestPastVectors = 
                  MatrixUtils.makeDelayEmbeddingVector(observations, k_candidate,
                      tau_candidate, (k_candidate-1)*tau_candidate,
                      observations.length - (k_candidate-1)*tau_candidate - 1);
              double[][] currentDestNextVectors =
                  MatrixUtils.makeDelayEmbeddingVector(observations, 1,
                      (k_candidate-1)*tau_candidate + 1,
                      observations.length - (k_candidate-1)*tau_candidate - 1);
              miCalc.addObservations(currentDestPastVectors, currentDestNextVectors);
            }

          } else {

            for (double[][] observations : vectorOfMultiVariateObservationTimeSeries) {
              double[][] currentDestPastVectors = 
                  MatrixUtils.makeDelayEmbeddingVector(observations, k_candidate,
                      tau_candidate, (k_candidate-1)*tau_candidate,
                      observations.length - (k_candidate-1)*tau_candidate - 1);
              double[][] currentDestNextVectors =
                  MatrixUtils.makeDelayEmbeddingVector(observations, 1,
                      (k_candidate-1)*tau_candidate + 1,
                      observations.length - (k_candidate-1)*tau_candidate - 1);
              miCalc.addObservations(currentDestPastVectors, currentDestNextVectors);
            }

          }
					miCalc.finaliseAddObservations();
					// Now grab the AIS estimate here -- it's already bias corrected for 
					double thisAIS = miCalc.computeAverageLocalOfObservations();
					if (debug) {
						System.out.printf("AIS for k=%d,tau=%d is %.3f\n",
								k_candidate, tau_candidate, thisAIS);
					}
					if (thisAIS > bestAIS) {
						// This parameter setting is the best so far:
						bestAIS = thisAIS;
						k_candidate_best = k_candidate;
						tau_candidate_best = tau_candidate;
					}
					if (k_candidate == 1) {
						// tau is irrelevant, so no point testing other values
						break;
					}
				}
			}
		} else {
			throw new RuntimeException("Unexpected value " + autoEmbeddingMethod +
					" for property " + PROP_AUTO_EMBED_METHOD);
		}
		// Make sure the embedding length and delay are set here
		k = k_candidate_best;
		tau = tau_candidate_best;
		if (debug) {
			System.out.printf("Embedding parameters set to k=%d,tau=%d (for prediction error %.3f)\n",
				k, tau, bestPredictionError);
		}
	}

}

