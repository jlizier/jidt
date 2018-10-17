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

import java.util.Iterator;
import java.util.Vector;

import infodynamics.measures.continuous.kraskov.ActiveInfoStorageCalculatorKraskov;
import infodynamics.measures.continuous.kraskov.ActiveInfoStorageCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;

/**
 * An Active Information Storage (AIS) calculator (implementing {@link ActiveInfoStorageCalculator})
 * which is affected using a 
 * given Mutual Information (MI) calculator (implementing
 * {@link MutualInfoCalculatorMultiVariate}) to make the calculations.
 * 
 * <p>Usage is as per the paradigm outlined for {@link ActiveInfoStorageCalculator},
 * except that in the constructor(s) for this class the implementation for
 * a {@link MutualInfoCalculatorMultiVariate} must be supplied.
 * Further properties may be set on this class via {@link #setProperty(String, String)}
 * (including auto-embedding parameters) as described in the javadocs for that method.
 * </p>
 * 
 * <p>This class <i>may</i> be used directly, however users are advised that
 * several child classes are available which already plug-in the various MI estimators
 * to provide AIS calculators (taking specific caution associated with
 * each type of estimator):</p>
 * <ul>
 * 	<li>{@link infodynamics.measures.continuous.gaussian.ActiveInfoStorageCalculatorGaussian}</li>
 * 	<li>{@link infodynamics.measures.continuous.kernel.ActiveInfoStorageCalculatorKernel}</li>
 * 	<li>{@link infodynamics.measures.continuous.kraskov.ActiveInfoStorageCalculatorKraskov}</li>
 * </ul>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>J.T. Lizier, M. Prokopenko and A.Y. Zomaya,
 * 		<a href="http://dx.doi.org/10.1016/j.ins.2012.04.016">
 * 		"Local measures of information storage in complex distributed computation"</a>,
 * 		Information Sciences, vol. 208, pp. 39-54, 2012.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see ActiveInfoStorageCalculator
 */
public class ActiveInfoStorageCalculatorViaMutualInfo implements
		ActiveInfoStorageCalculator {

	/**
	 * The underlying mutual information calculator
	 */
	protected MutualInfoCalculatorMultiVariate miCalc;
	/**
	 * Length of past history to consider (embedding length)
	 */
	protected int k = 1;
	/**
	 * Embedding delay to use between elements of the embeding vector.
	 * We're hard-coding a delay of 1 between the history vector and the next 
	 *  observation however.
	 */
	protected int tau = 1;
	/**
	 * Whether debug mode is on
	 */
	protected boolean debug = false;
	
	/**
	 * Storage for observations supplied via {@link #addObservations(double[])}
	 * type calls
	 */
	protected Vector<double[]> vectorOfObservationTimeSeries;
	/**
	 * Storage for validity arrays for supplied observations.
	 * Entries are null where the whole corresponding observation time-series is valid
	 */
	protected Vector<boolean[]> vectorOfValidityOfObservations;

	/**
	 * Property name for the auto-embedding method. Defaults to {@link #AUTO_EMBED_METHOD_NONE}.
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
	 * Construct using an instantiation of the named MI calculator
	 * 
	 * @param miCalculatorClassName fully qualified class name of the MI calculator to instantiate
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 */
	public ActiveInfoStorageCalculatorViaMutualInfo(String miCalculatorClassName)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		@SuppressWarnings("unchecked")
		Class<MutualInfoCalculatorMultiVariate> miClass = 
				(Class<MutualInfoCalculatorMultiVariate>) Class.forName(miCalculatorClassName);
		MutualInfoCalculatorMultiVariate miCalc = miClass.newInstance();
		construct(miCalc);
	}
	
	/**
	 * Construct using an instantiation of the given MI class
	 * 
	 * @param miCalcClass Class of the MI calculator to instantiate and use
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 */
	protected ActiveInfoStorageCalculatorViaMutualInfo(Class<MutualInfoCalculatorMultiVariate> miCalcClass)
			throws InstantiationException, IllegalAccessException {
		MutualInfoCalculatorMultiVariate miCalc = miCalcClass.newInstance();
		construct(miCalc);
	}
	
	/**
	 * Construct using the given (constructed but not initialised)
	 * MI calculator.
	 * 
	 * @param miCalc MI calculator which is already constructed but
	 *  there has not been a call to its {@link MutualInfoCalculatorMultiVariate#initialise()}
	 *  method yet 
	 */
	protected ActiveInfoStorageCalculatorViaMutualInfo(MutualInfoCalculatorMultiVariate miCalc) {
		construct(miCalc);
	}

	/**
	 * Internal routine to execute common code for constructing an instance
	 * 
	 * @param miCalc
	 */
	protected void construct(MutualInfoCalculatorMultiVariate miCalc) {
		this.miCalc = miCalc;
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#initialise()
	 */
	@Override
	public void initialise() throws Exception {
		initialise(k, tau); // Initialise with current value of k
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#initialise(int)
	 */
	@Override
	public void initialise(int k) throws Exception {
		initialise(k, tau);
	}

	/**
	 * {@inheritDoc}
	 * 
	 * <p>All child classes <b>must</b> call this routine on this as the super class
	 *  once they have finished executing their specialised code
	 *  for their {@link #initialise()} implementations.
	 * </p>
	 * 
	 */
	@Override
	public void initialise(int k, int tau) throws Exception {
		this.k = k;
		this.tau = tau;
		vectorOfObservationTimeSeries = null;
		vectorOfValidityOfObservations = null;
	}

	/**
	 * Set properties for the underlying calculator implementation.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Allowable property names include:</p>
	 * <ul>
	 * 	<li>Those defined for the {@link ActiveInfoStorageCalculator} interface
	 *      (i.e. {@link #K_PROP_NAME} or {@link #TAU_PROP_NAME})</li>
	 * 	<li>{@link #PROP_AUTO_EMBED_METHOD} -- method by which the calculator
	 * 	   automatically determines the embedding history length ({@link #K_PROP_NAME})
	 * 	   and embedding delay ({@link #TAU_PROP_NAME}). Default is {@link #AUTO_EMBED_METHOD_NONE} meaning
	 * 	   values are set manually; other accepted values include: {@link #AUTO_EMBED_METHOD_RAGWITZ} for use
	 * 		of the Ragwitz criteria and {@link #AUTO_EMBED_METHOD_MAX_CORR_AIS} for using
	 * 		the maz bias-corrected AIS criteria (both searching up to {@link #PROP_K_SEARCH_MAX} and 
	 * 		{@link #PROP_TAU_SEARCH_MAX}, as outlined by Garland et al. in the references list above).
	 * 		Use of any value other than {@link #AUTO_EMBED_METHOD_NONE}
	 * 		will lead to any previous settings for k and tau (via e.g. {@link #initialise(int, int)} or
	 * 		auto-embedding during previous calculations) will be overwritten after observations
	 * 		are supplied.</li>
	 * 	<li>{@link #PROP_K_SEARCH_MAX} -- maximum embedded history length to search
	 * 		up to if automatically determining the embedding parameters (as set by
	 * 		{@link #PROP_AUTO_EMBED_METHOD}); default is 1</li>
	 * 	<li>{@link #PROP_TAU_SEARCH_MAX} -- maximum embedded history length to search
	 * 		up to if automatically determining the embedding parameters (as set by
	 * 		{@link #PROP_AUTO_EMBED_METHOD}); default is 1</li>
	 * 	<li>{@link #PROP_RAGWITZ_NUM_NNS} -- number of nearest neighbours to use
	 * 		in the auto-embedding if the property {@link #PROP_AUTO_EMBED_METHOD}
	 * 		has been set to {@link #AUTO_EMBED_METHOD_RAGWITZ}. Defaults to the property value
	 *      set for {@link MutualInfoCalculatorMultiVariateKraskov.PROP_K}</li>
	 *  <li>Any properties defined for the underlying
	 *     {@link MutualInfoCalculatorMultiVariate#setProperty(String, String)} implementation,
	 *     <b>however</b> the user is <b>not</b> allowed to set the property 
	 *     {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} here.
	 *     This would set a time difference from the history vector to the next
	 *     step, which we currently do not allow.
	 *     (If we change our mind one day and allow it, we could implement
	 *     it simply by letting the time diff property be set here).</li>
	 * </ul>
	 *  
	 * <p>Note that implementing classes may defined additional properties.</p>
	 *
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		if (propertyName.equalsIgnoreCase(MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF)) {
			throw new Exception("Cannot set " + MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF
					+ " property on the ActiveInfoStorageCalculator");
		}
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(TAU_PROP_NAME)) {
			tau = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_AUTO_EMBED_METHOD)) {
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
		} else {
			// No property was set on this class, assume it is for the underlying
			//  MI calculator
			miCalc.setProperty(propertyName, propertyValue);
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public String getProperty(String propertyName)
			throws Exception {
		
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			return Integer.toString(k);
		} else if (propertyName.equalsIgnoreCase(TAU_PROP_NAME)) {
			return Integer.toString(tau);
		} else if (propertyName.equalsIgnoreCase(PROP_AUTO_EMBED_METHOD)) {
			return autoEmbeddingMethod;
		} else if (propertyName.equalsIgnoreCase(PROP_K_SEARCH_MAX)) {
			return Integer.toString(k_search_max);
		} else if (propertyName.equalsIgnoreCase(PROP_TAU_SEARCH_MAX)) {
			return Integer.toString(tau_search_max);
		} else if (propertyName.equalsIgnoreCase(PROP_RAGWITZ_NUM_NNS)) {
			if (ragwitz_num_nns_set) {
				return Integer.toString(ragwitz_num_nns);
			} else {
				return miCalc.getProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_K);
			}
		} else {
			// No property was set on this class, assume it is for the underlying
			//  MI calculator, even if it is for 
			//  MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF which
			//  is not a valid property for the AIS calculator:
			return miCalc.getProperty(propertyName);
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#setObservations(double[])
	 */
	@Override
	public void setObservations(double[] observations) throws Exception {
		startAddObservations();
		addObservations(observations);
		finaliseAddObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#startAddObservations()
	 */
	@Override
	public void startAddObservations() {
		vectorOfObservationTimeSeries = new Vector<double[]>();
		vectorOfValidityOfObservations = new Vector<boolean[]>();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#addObservations(double[])
	 */
	@Override
	public void addObservations(double[] observations) throws Exception {
		// Store these observations in our vector for now
		vectorOfObservationTimeSeries.add(observations);
		vectorOfValidityOfObservations.add(null); // All observations were valid
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#addObservations(double[])
	 */
	@Override
	public void addObservations(double[] observations, boolean[] valid) throws Exception {
		// Store these observations in our vector for now
		vectorOfObservationTimeSeries.add(observations);
		vectorOfValidityOfObservations.add(valid); // All observations were valid
	}

	/**
	 * Protected method to internally parse and submit observations through
	 *  to the supplied MI calculator with the given embedding parameters
	 * 
	 * @param miCalc_in_use MI calculator to supply 
	 * @param k_in_use k embedding dimension to use
	 * @param tau_in_use tau embedding delay to use
	 * @param observations time series of observations
	 * @throws Exception
	 */
	protected void addObservationsWithGivenParams(MutualInfoCalculatorMultiVariate miCalc_in_use,
			int k_in_use, int tau_in_use, double[] observations) throws Exception {
		if (observations.length - (k_in_use-1)*tau_in_use - 1 <= 0) {
			// There are no observations to add here
			// Don't throw an exception, do nothing since more observations
			//  can be added later.
			return;
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(observations, k_in_use, tau_in_use,
						(k_in_use-1)*tau_in_use, observations.length - (k_in_use-1)*tau_in_use - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(observations, 1, (k_in_use-1)*tau_in_use + 1,
						observations.length - (k_in_use-1)*tau_in_use - 1);
		miCalc_in_use.addObservations(currentDestPastVectors, currentDestNextVectors);
	}
	
	/**
	 * Protected method to internally parse and submit observations through
	 *  to the supplied MI calculator with the given embedding parameters.
	 * This is done given a time-series of booleans indicating whether each entry
	 *  is valid
	 * 
	 * @param miCalc_in_use MI calculator to supply 
	 * @param k_in_use k embedding dimension to use
	 * @param tau_in_use tau embedding delay to use
	 * @param observations time series of observations
	 * @param valid a time series (with indices the same as observations) indicating
	 *  whether the entry in observations at that index is valid; we only take vectors
	 *  as samples to add to the observation set where all points in the time series
	 *  (even between points in the embedded k-vector with embedding delays) are valid.
	 * @throws Exception
	 */
	protected void addObservationsWithGivenParams(MutualInfoCalculatorMultiVariate miCalc_in_use,
			int k_in_use, int tau_in_use, double[] observations, boolean[] valid) throws Exception {
		
		// compute the start and end times using our determined embedding parameters:
		Vector<int[]> startAndEndTimePairs = computeStartAndEndTimePairs(k_in_use, tau_in_use, valid);
		
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservationsWithGivenParams(miCalc_in_use, k_in_use, tau_in_use,
					MatrixUtils.select(observations, startTime, endTime - startTime + 1));
		}
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#addObservations(double[], int, int)
	 */
	@Override
	public void addObservations(double[] observations, int startTime,
			int numTimeSteps) throws Exception {
		addObservations(MatrixUtils.select(observations, startTime, numTimeSteps));
	}

	/**
	 * Hook in case child implementations need to perform any processing on the 
	 *  observation time series prior to their being processed and supplied
	 *  to the underlying MI calculator.
	 * Primarily this is to allow the child implementation to automatically determine
	 *  embedding parameters if desired, and a default implementation is provided
	 *  for this for the main two auto-embedding methods.
	 * Child implementations do not need to override this default empty implementation
	 *  if no new functionality is required.
	 */
	protected void preFinaliseAddObservations() throws Exception {
		// Automatically determine the embedding parameters for the given time series
		
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_NONE)) {
			return;
		}
		// Else we need to auto embed
		
		// TODO Could make sure the rest of the code could handle k=0
		//  as default if nothing can improve on this, though 
		//  I think I prefer k=1 to stay as default.
		int k_candidate_best = 1;
		int tau_candidate_best = 1;
		
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ)) {
			double bestPredictionError = Double.POSITIVE_INFINITY;
			if (debug) {
				System.out.printf("Beginning Ragwitz auto-embedding with k_max=%d, tau_max=%d\n",
						k_search_max, tau_search_max);
			}
			
			for (int k_candidate = 1; k_candidate <= k_search_max; k_candidate++) {
				for (int tau_candidate = 1; tau_candidate <= tau_search_max; tau_candidate++) {
					try {
						// Use a KSG MI calculator, which can do Ragwitz fairly easily.
						MutualInfoCalculatorMultiVariateKraskov miCalcKraskov;
						if ((this instanceof ActiveInfoStorageCalculatorKraskov) || 
							(this instanceof ActiveInfoStorageCalculatorMultiVariateKraskov)) {
							// Use our internal MI calculator in case it has any particular 
							//  properties we need to have been set already
							miCalcKraskov = (MutualInfoCalculatorMultiVariateKraskov) miCalc;
						} else {
							// We'll create one to use, but we won't give the user the opportunity to set most of the properties
							//  on it, just the number of nearest neighbours. Leave NORM_TYPE etc as default.
							miCalcKraskov = new MutualInfoCalculatorMultiVariateKraskov1();
						}
						prepareMICalculator(miCalcKraskov, k_candidate, tau_candidate);
						// Now grab the prediction errors of the next value from the required number of
						// nearest neighbours of the previous state: (array is of only one term)
						double[] predictionError;
						if (ragwitz_num_nns_set) {
							predictionError = 
								miCalcKraskov.computePredictionErrorsFromObservations(false, ragwitz_num_nns);
						} else {
							predictionError = 
								miCalcKraskov.computePredictionErrorsFromObservations(false);
						}
						if (debug) {
							System.out.printf("Embedding prediction error (dim=%d) for k=%d,tau=%d is %.3f\n",
									predictionError.length, k_candidate, tau_candidate,
									predictionError[0] / (double) miCalcKraskov.getNumObservations());
						}
						if ((predictionError[0] / (double) miCalcKraskov.getNumObservations())
								< bestPredictionError) {
							// This parameter setting is the best so far:
							// (Note division by number of observations to normalise
							//  for less observations for larger k and tau)
							bestPredictionError = predictionError[0] / (double) miCalcKraskov.getNumObservations();
							k_candidate_best = k_candidate;
							tau_candidate_best = tau_candidate;
						}
						if (k_candidate == 1) {
							// tau is irrelevant, so no point testing other values
							break;
						}
					} catch (Exception ex) {
						throw new Exception("Exception encountered in attempting auto-embedding, evaluating candidates k=" + k_candidate +
								", tau=" + tau_candidate, ex);
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
					try {
						// Use our internal MI calculator in case it has any particular 
						//  properties we need to have been set already
						prepareMICalculator(miCalc, k_candidate, tau_candidate);
						// Now grab the AIS estimate here
						double thisAIS = miCalc.computeAverageLocalOfObservations();
						thisAIS -= computeAdditionalBiasToRemove();
						if (debug) {
							System.out.printf("AIS (bias corrected) for k=%d,tau=%d (%d samples) is %.5f\n",
									k_candidate, tau_candidate, miCalc.getNumObservations(), thisAIS);
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
					} catch (Exception ex) {
						throw new Exception("Exception encountered in attempting auto-embedding, evaluating candidates k=" + k_candidate +
								", tau=" + tau_candidate, ex);
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
			System.out.printf("Embedding parameters set to k=%d,tau=%d\n",
				k, tau);
		}
	}
	
	/**
	 * Internal method to compute any additional bias correction in the underlying calculator
	 * during auto-embedding in {@link #preFinaliseAddObservations()} if required.
	 * 
	 * @return additional bias correction to remove (will be zero if assumed to be already bias corrected).
	 * @throws Exception 
	 */
	protected double computeAdditionalBiasToRemove() throws Exception {
		// Default implementation does nothing
		return 0;
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		// Auto embed if required
		preFinaliseAddObservations();
		
		prepareMICalculator(miCalc, k, tau);
		
		vectorOfObservationTimeSeries = null; // No longer required
		vectorOfValidityOfObservations = null;
	}

	/**
	 * Prepare the given pre-instantiated (and properties supplied)
	 *  Mutual information calculator with this data set,
	 *  using the embedding parameters supplied.
	 * This may be used in the final calculation, or by the auto-embedding 
	 *  procedures, hence the use of method arguments rather than
	 *  using the member variables directly.
	 * 
	 * @param miCalc_in_use MI calculator to supply 
	 * @param k_in_use k embedding dimension to use
	 * @param tau_in_use tau embedding delay to use
	 * @throws Exception
	 */
	protected void prepareMICalculator(MutualInfoCalculatorMultiVariate miCalc_in_use,
			int k_in_use, int tau_in_use) throws Exception {
		
		// Initialise the MI calculator, including any auto-embedding length
		miCalc_in_use.initialise(k_in_use, 1);
		miCalc_in_use.startAddObservations();
		// Send all of the observations through:
		Iterator<boolean[]> validityIterator = vectorOfValidityOfObservations.iterator();
		for (double[] observations : vectorOfObservationTimeSeries) {
			boolean[] validity = validityIterator.next();
			if (validity == null) {
				// Add the whole time-series
				addObservationsWithGivenParams(miCalc_in_use, k_in_use,
						tau_in_use, observations);
			} else {
				addObservationsWithGivenParams(miCalc_in_use, k_in_use,
						tau_in_use, observations, validity);
			}
		}
		// TODO do we need to throw an exception if there are no observations to add?
		miCalc_in_use.finaliseAddObservations();
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#setObservations(double[], boolean[])
	 */
	@Override
	public void setObservations(double[] observations, boolean[] valid)
			throws Exception {
		startAddObservations();
		// Add these observations and the indication of their validity
		vectorOfObservationTimeSeries.add(observations);
		vectorOfValidityOfObservations.add(valid);
		finaliseAddObservations();
	}

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid series of observations.
	 * 
	 * <p>This method is made public so it can be used if one wants to compute the number of
	 *  observations prior to making a call to {@link #setObservations(double[], boolean[])}.</p>
	 * 
	 * <p>Functions as per {@link #computeStartAndEndTimePairs(int, int, boolean[])}
	 * with k and tau set to the property values of this calculator.</p>
	 * 
	 * @param valid a time series (with indices the same as observations)
	 *  indicating whether the entry in observations at that index is valid; 
	 *  we only take vectors as samples to add to the observation set where
	 *  all points in the time series (even between points in 
	 *  the embedded k-vector with embedding delays) are valid.
	 * @return a vector for start and end time pairs of valid series
	 *  of observations (as defined by <code>valid</code>).
	 */
	public Vector<int[]> computeStartAndEndTimePairs(boolean[] valid) {
		return computeStartAndEndTimePairs(k, tau, valid);
	}

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid series of observations.
	 * 
	 * <p>This method is made public so it can be used if one wants to compute the number of
	 *  observations prior to making a call to {@link #setObservations(double[], boolean[])}.</p>
	 * 
	 * @param k_in_use the k embedding dimension parameter to use here
	 * @param tau_in_use the tau embedding delay parameter to use here
	 * @param valid a time series (with indices the same as observations)
	 *  indicating whether the entry in observations at that index is valid; 
	 *  we only take vectors as samples to add to the observation set where
	 *  all points in the time series (even between points in 
	 *  the embedded k-vector with embedding delays) are valid.
	 * @return a vector for start and end time pairs of valid series
	 *  of observations (as defined by <code>valid</code>).
	 */
	public Vector<int[]> computeStartAndEndTimePairs(int k_in_use, int tau_in_use, boolean[] valid) {
		// Scan along the data avoiding invalid values
		int startTime = 0;
		int endTime = 0;
		boolean lookingForStart = true;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();
		for (int t = 0; t < valid.length; t++) {
			if (lookingForStart) {
				// Precondition: startTime holds a candidate start time
				if (valid[t]) {
					// This point is OK at the destination
					if (t - startTime < (k_in_use-1)*tau_in_use+1) {
						// We're still checking the past history only, so
						continue;
					} else {
						// We've got the full past history ok
						// set a candidate endTime
						endTime = t;
						lookingForStart = false;
						if (t == valid.length - 1) {
							// we need to terminate now
							int[] timePair = new int[2];
							timePair[0] = startTime;
							timePair[1] = endTime;
							startAndEndTimePairs.add(timePair);
							// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
						}
					}
				} else {
					// We need to keep looking.
					// Move the potential start time to the next point
					startTime = t + 1;
				}
			} else {
				// Precondition: startTime holds the start time for this set, 
				//  endTime holds a candidate end time
				// Check if we can include the current time step
				boolean terminateSequence = false;
				if (valid[t]) {
					// We can extend
					endTime = t;
				} else {
					terminateSequence = true;
				}
				if (t == valid.length - 1) {
					// we need to terminate the sequence anyway
					terminateSequence = true;
				}
				if (terminateSequence) {
					// This section is done
					int[] timePair = new int[2];
					timePair[0] = startTime;
					timePair[1] = endTime;
					startAndEndTimePairs.add(timePair);
					// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
					lookingForStart = true;
					startTime = t + 1;
				}
			}
		}
		return startAndEndTimePairs;
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeAverageLocalOfObservations()
	 */
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		return miCalc.computeAverageLocalOfObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeLocalOfPreviousObservations()
	 */
	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] local = miCalc.computeLocalOfPreviousObservations();
		if (!miCalc.getAddedMoreThanOneObservationSet()) {
			double[] localsToReturn = new double[local.length + (k-1)*tau + 1];
			System.arraycopy(local, 0, localsToReturn, (k-1)*tau + 1, local.length);
			return localsToReturn;
		} else {
			return local;
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeLocalUsingPreviousObservations(double[])
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(double[] newObservations) throws Exception {
		if (newObservations.length - (k-1)*tau - 1 <= 0) {
			// There are no observations to compute for here
			return new double[newObservations.length];
		}
		double[][] newDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newObservations, k, tau, (k-1)*tau, newObservations.length - (k-1)*tau - 1);
		double[][] newDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(newObservations, 1, (k-1)*tau + 1, newObservations.length - (k-1)*tau - 1);
		double[] local = miCalc.computeLocalUsingPreviousObservations(newDestPastVectors, newDestNextVectors);
		// Pad the front of the array with zeros where local AIS isn't defined:
		double[] localsToReturn = new double[local.length + (k-1)*tau + 1];
		System.arraycopy(local, 0, localsToReturn, (k-1)*tau + 1, local.length);
		return localsToReturn;

	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeSignificance(int)
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		return miCalc.computeSignificance(numPermutationsToCheck);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeSignificance(int[][])
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		return miCalc.computeSignificance(newOrderings);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#setDebug(boolean)
	 */
	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
		miCalc.setDebug(debug);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#getLastAverage()
	 */
	@Override
	public double getLastAverage() {
		return miCalc.getLastAverage();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#getNumObservations()
	 */
	@Override
	public int getNumObservations() throws Exception {
		return miCalc.getNumObservations();
	}
}
