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

import infodynamics.measures.continuous.gaussian.ActiveInfoStorageCalculatorGaussian;
import infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian;
import infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian;
import infodynamics.measures.continuous.kraskov.ActiveInfoStorageCalculatorKraskov;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov1;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov2;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;

import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

/**
 * A Transfer Entropy (TE) calculator (implementing {@link TransferEntropyCalculator})
 * which is affected using a 
 * given Conditional Mutual Information (MI) calculator (implementing
 * {@link ConditionalMutualInfoCalculatorMultiVariate}) to make the calculations.
 * 
 * <p>Usage is as per the paradigm outlined for {@link TransferEntropyCalculator},
 * except that in the constructor(s) for this class the implementation for
 * a {@link ConditionalMutualInfoCalculatorMultiVariate} must be supplied.
 * </p>
 * 
 * <p>This class <i>may</i> be used directly, however users are advised that
 * several child classes are available which already plug-in the various
 * conditional MI estimators
 * to provide TE calculators (taking specific caution associated with
 * each type of estimator):</p>
 * <ul>
 * 	<li>{@link infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian}</li>
 * 	<li>{@link infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov}</li>
 * </ul>
 * 
 * <p>Embedding parameters may be automatically determined as per the Ragwitz criteria
 * by setting the property {@link #PROP_AUTO_EMBED_METHOD} to {@link #AUTO_EMBED_METHOD_RAGWITZ}
 * or {@link #AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY},
 * or as per the max. bias-corrected AIS criteria by 
 * setting the property {@link #PROP_AUTO_EMBED_METHOD} to {@link #AUTO_EMBED_METHOD_MAX_CORR_AIS}
 * (as per Garland et al. in the reference list)
 * plus additional parameter settings for these.
 * </p>
 *      
 * TODO Delete TransferEntropyCalculatorCommon once we've switched everything over to use this?
 * Might be useful to leave it after all, and move common functionality from here to there.
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 * 	<li>Ragwitz and Kantz, "Markov models from data by simple nonlinear time series
 *  	predictors in delay embedding spaces", Physical Review E, vol 65, 056201 (2002).</li>
 *  <li>J. Garland, R. G. James, E. Bradley, <a href="http://dx.doi.org/10.1103/physreve.93.022221">
 *  	"Leveraging information storage to select forecast-optimal parameters for delay-coordinate reconstructions"</a>,
 *  	Physical Review E, Vol. 93 (2016), 022221, doi:</li>
 * </ul>
 *  
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 */
public class TransferEntropyCalculatorViaCondMutualInfo implements
		TransferEntropyCalculator {

	/**
	 * Underlying conditional mutual information calculator
	 */
	protected ConditionalMutualInfoCalculatorMultiVariate condMiCalc;
	/**
	 * Length of past destination history to consider (embedding length)
	 */
	protected int k = 1;
	/**
	 * Embedding delay to use between elements of the destination embeding vector.
	 * We're hard-coding a delay of 1 between the history vector and the next 
	 *  observation however.
	 */
	protected int k_tau = 1;
	/**
	 * Length of past source history to consider (embedding length)
	 */
	protected int l = 1;
	/**
	 * Embedding delay to use between elements of the source embeding vector.
	 */
	protected int l_tau = 1;
	/**
	 * Source-destination next observation delay
	 */
	protected int delay = 1;
	/**
	 * Time index of the last point in the destination embedding of the first
	 *  (destination past, source past, destination next) tuple than can be 
	 *  taken from any set of time-series observations. 
	 */
	protected int startTimeForFirstDestEmbedding;

	/**
	 * Whether we're in debugging mode
	 */
	protected boolean debug = false;

	/**
	 * Storage for source observations supplied via {@link #addObservations(double[], double[])} etc.
	 */
	protected Vector<double[]> vectorOfSourceTimeSeries;

	/**
	 * Storage for destination observations supplied via {@link #addObservations(double[], double[])} etc.
	 */
	protected Vector<double[]> vectorOfDestinationTimeSeries;

	/**
	 * Storage for validity arrays for supplied source observations.
	 * Entries are null where the whole corresponding observation time-series is valid
	 */
	protected Vector<boolean[]> vectorOfValidityOfSource;

	/**
	 * Storage for validity arrays for supplied destination observations.
	 * Entries are null where the whole corresponding observation time-series is valid
	 */
	protected Vector<boolean[]> vectorOfValidityOfDestination;

	/**
	 * Array of the number of observations added by each separate call to
	 * {@link #addObservations(double[], double[])}, in order of which those calls
	 * were made. This is returned by {@link #getSeparateNumObservations()}
	 */
	protected int[] separateNumObservations;
	
	/**
	 * Property name for specifying which (if any) auto-embedding method to use.
	 * Valid values include {@link #AUTO_EMBED_METHOD_RAGWITZ}, {@link #AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY},
	 * {@link #AUTO_EMBED_METHOD_MAX_CORR_AIS}, {@link #AUTO_EMBED_METHOD_MAX_CORR_AIS_DEST_ONLY},
	 * {@link #AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE} and {@link #AUTO_EMBED_METHOD_NONE}.
	 * Defaults to {@link #AUTO_EMBED_METHOD_NONE}
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
	 *  for both source and destination time-series
	 */
	public static final String AUTO_EMBED_METHOD_RAGWITZ = "RAGWITZ";
	/**
	 * Valid value for the property {@link #PROP_AUTO_EMBED_METHOD} indicating that
	 *  the Ragwitz optimisation technique should be used for automatic embedding
	 *  for the destination time-series only
	 */
	public static final String AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY = "RAGWITZ_DEST_ONLY";
	/**
	 * Valid value for the property {@link #PROP_AUTO_EMBED_METHOD} indicating that
	 *  the automatic embedding should be done by maximising the bias corrected
	 *  AIS, for both source and destination time series
	 */
	public static final String AUTO_EMBED_METHOD_MAX_CORR_AIS = "MAX_CORR_AIS";
	/**
	 * Valid value for the property {@link #PROP_AUTO_EMBED_METHOD} indicating that
	 *  the automatic embedding should be done by maximising the bias corrected
	 *  AIS for the target and subsequently maximising the bias-corrected TE over source embeddings,
	 *  given a fixed source-target delay.
	 */
	public static final String AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE = "MAX_CORR_AIS_AND_TE";
	/**
	 * Valid value for the property {@link #PROP_AUTO_EMBED_METHOD} indicating that
	 *  the automatic embedding should be done by maximising the bias corrected
	 *  AIS, for destination time series only
	 */
	public static final String AUTO_EMBED_METHOD_MAX_CORR_AIS_DEST_ONLY = "MAX_CORR_AIS_DEST_ONLY";
	/**
	 * Internal variable tracking what type of auto embedding (if any)
	 *  we are using
	 */
	protected String autoEmbeddingMethod = AUTO_EMBED_METHOD_NONE;
	
	/**
	 * Property name for maximum embedding lengths (i.e. k for destination, and l for source if we're auto-embedding
	 *  the source as well) for the auto-embedding search. Defaults to 1
	 */
	public static final String PROP_K_SEARCH_MAX = "AUTO_EMBED_K_SEARCH_MAX";
	/**
	 * Internal variable for storing the maximum embedding length to search up to for
	 *  automating the parameters.
	 */
	protected int k_search_max = 1;

	/**
	 * Property name for maximum embedding delay (i.e. k_tau for destination, and l_tau for source if we're auto-embedding
	 *   the source as well) for the auto-embedding search. Defaults to 1
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
	protected int ragwitz_num_nns = 4;
	/** 
	 * Internal variable to track whether the property {@link #PROP_RAGWITZ_NUM_NNS} has been
	 * set yet
	 */
	protected boolean ragwitz_num_nns_set = false;

	/**
	 * Storage for the properties ready to pass onto the underlying conditional MI calculators should they change 
	 */
	protected Hashtable<String,String> props = new Hashtable<String,String>();
	
	/**
	 * Construct a transfer entropy calculator using an instance of
	 * condMiCalculatorClassName as the underlying conditional mutual information calculator.
	 * 
	 * @param condMiCalculatorClassName fully qualified name of the class which must implement
	 * 	{@link ConditionalMutualInfoCalculatorMultiVariate}
	 * @throws InstantiationException if the given class cannot be instantiated
	 * @throws IllegalAccessException if illegal access occurs while trying to create an instance
	 *   of the class
	 * @throws ClassNotFoundException if the given class is not found
	 */
	public TransferEntropyCalculatorViaCondMutualInfo(String condMiCalculatorClassName)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		@SuppressWarnings("unchecked")
		Class<ConditionalMutualInfoCalculatorMultiVariate> condMiClass = 
				(Class<ConditionalMutualInfoCalculatorMultiVariate>) Class.forName(condMiCalculatorClassName);
		ConditionalMutualInfoCalculatorMultiVariate condMiCalc = condMiClass.newInstance();
		construct(condMiCalc);
	}

	/**
	 * Construct a transfer entropy calculator using an instance of
	 * condMiCalcClass as the underlying conditional mutual information calculator.
	 * 
	 * @param condMiCalcClass the class which must implement
	 * 	{@link ConditionalMutualInfoCalculatorMultiVariate}
	 * @throws InstantiationException if the given class cannot be instantiated
	 * @throws IllegalAccessException if illegal access occurs while trying to create an instance
	 *   of the class
	 * @throws ClassNotFoundException if the given class is not found
	 */
	public TransferEntropyCalculatorViaCondMutualInfo(Class<ConditionalMutualInfoCalculatorMultiVariate> condMiCalcClass)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		ConditionalMutualInfoCalculatorMultiVariate condMiCalc = condMiCalcClass.newInstance();
		construct(condMiCalc);
	}

	/**
	 * Construct this calculator by passing in a constructed but not initialised
	 * underlying Conditional Mutual information calculator.
	 * 
	 * @param condMiCalc An instantiated conditional mutual information calculator.
	 * @throws Exception if the supplied calculator has not yet been instantiated.
	 */
	public TransferEntropyCalculatorViaCondMutualInfo(ConditionalMutualInfoCalculatorMultiVariate condMiCalc) throws Exception {
		if (condMiCalc == null) {
			throw new Exception("Conditional MI calculator used to construct ConditionalTransferEntropyCalculatorViaCondMutualInfo " +
					" must have already been instantiated.");
		}
		construct(condMiCalc);
	}
	
	/**
	 * Internal method to set the conditional mutual information calculator.
	 * Can be overridden if anything else needs to be done with it by the child classes.
	 * 
	 * @param condMiCalc
	 */
	protected void construct(ConditionalMutualInfoCalculatorMultiVariate condMiCalc) {
		this.condMiCalc = condMiCalc;
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorCommon#initialise()
	 */
	@Override
	public void initialise() throws Exception {
		initialise(k, k_tau, l, l_tau, delay);
	}
	
	@Override
	public void initialise(int k) throws Exception {
		initialise(k, k_tau, l, l_tau, delay);
	}
	
	/**
	 * Initialise the calculator for re-use with new observations.
	 * New embedding parameters and source-destination delay
	 * may be supplied here; all other parameters
	 * remain unchanged.
	 * 
	 * @param k embedding length of destination past history to consider
	 * @param k_tau embedding delay for the destination variable
	 * @param l embedding length of source past history to consider
	 * @param l_tau embedding delay for the source variable
	 * @param delay time lag between last element of source and destination next value
	 */
	public void initialise(int k, int k_tau, int l, int l_tau, int delay) throws Exception {
		if (delay < 0) {
			throw new Exception("Cannot compute TE with source-destination delay < 0");
		}
		this.k = k;
		this.k_tau = k_tau;
		this.l = l;
		this.l_tau = l_tau;
		this.delay = delay;
		
		startTimeForFirstDestEmbedding =
				computeStartTimeForFirstDestEmbedding(k, k_tau, l, l_tau, delay);

		vectorOfSourceTimeSeries = null;
		vectorOfDestinationTimeSeries = null;
		vectorOfValidityOfSource = null;
		vectorOfValidityOfDestination = null;
		separateNumObservations = new int[] {};
	}

	/**
	 * Protected internal method to 
	 * set the point at which we can start taking observations from in any
	 * addObservations call.
	 * 
	 * User supplied parameters for k etc, so this can be used not only to
	 *  set the definitive startTimeForFirstDestEmbedding but also when
	 *  we're searching the parameter space in auto-embedding.
	 * 
	 */
	protected static int computeStartTimeForFirstDestEmbedding(
			int k_in_use, int k_tau_in_use, int l_in_use, int l_tau_in_use, int delay_in_use) {
		// These two integers represent the last
		// point of the destination embedding, in the cases where the destination
		// embedding itself determines where we can start taking observations, or
		// the case where the source embedding plus delay is longer and so determines
		// where we can start taking observations.
		int startTimeBasedOnDestPast = (k_in_use-1)*k_tau_in_use;
		int startTimeBasedOnSourcePast = (l_in_use-1)*l_tau_in_use + delay_in_use - 1;
		return Math.max(startTimeBasedOnDestPast, startTimeBasedOnSourcePast);
	}
	
	/**
	 * Sets properties for the TE calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #PROP_AUTO_EMBED_METHOD} -- method by which the calculator
	 * 		automatically determines the embedding history length ({@link #K_PROP_NAME})
	 * 		and embedding delay ({@link #TAU_PROP_NAME}) for destination and potentially source.
	 * 		Default is {@link #AUTO_EMBED_METHOD_NONE} meaning
	 * 		values are set manually; other accepted values include: {@link #AUTO_EMBED_METHOD_RAGWITZ} for use
	 * 		of the Ragwitz criteria for both source and destination (searching up to {@link #PROP_K_SEARCH_MAX} and 
	 * 		{@link #PROP_TAU_SEARCH_MAX}), and {@link #AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY} for use
	 * 		of the Ragwitz criteria for the destination only;
	 * 		{@link #AUTO_EMBED_METHOD_MAX_CORR_AIS} for use of the max bias corrected AIS criteria
	 * 		for both source and destination (searching up to {@link #PROP_K_SEARCH_MAX} and 
	 * 		{@link #PROP_TAU_SEARCH_MAX}), {@link #AUTO_EMBED_METHOD_MAX_CORR_AIS_DEST_ONLY} for use of
	 * 		this criteria for the destination only and {@link #AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE} for
	 * 		use of this criteria for the target, plus the max bias corrected TE for source embeddings.
	 * 		Use of any value other than {@link #AUTO_EMBED_METHOD_NONE}
	 * 		will lead to previous settings for embedding lengths and delays (via e.g. {@link #initialise(int, int)} or
	 * 		auto-embedding during previous calculations) for the destination and perhaps source to
	 * 		be overwritten after observations are supplied.</li>
	 * 		<li>{@link #PROP_K_SEARCH_MAX} -- maximum embedded history length to search
	 * 		up to if automatically determining the embedding parameters (as set by
	 * 		{@link #PROP_AUTO_EMBED_METHOD}) for the time-series to be embedded; default is 1</li>
	 * 		<li>{@link #PROP_TAU_SEARCH_MAX} -- maximum embedded history length to search
	 * 		up to if automatically determining the embedding parameters (as set by
	 * 		{@link #PROP_AUTO_EMBED_METHOD}) for the time-series to be embedded; default is 1</li>
	 * 		<li>{@link #PROP_RAGWITZ_NUM_NNS} -- number of nearest neighbours to use
	 * 		in the auto-embedding if the property {@link #PROP_AUTO_EMBED_METHOD}
	 * 		has been set to {@link #AUTO_EMBED_METHOD_RAGWITZ} or {@link #AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY}.
	 * 		Defaults to the property value
	 *      set for {@link ConditionalMutualInfoCalculatorMultiVariateKraskov#PROP_K}</li>
	 * 		<li>Any properties accepted by {@link TransferEntropyCalculator#setProperty(String, String)}</li>
	 * 		<li>Or properties accepted by the underlying
	 * 		{@link ConditionalMutualInfoCalculatorMultiVariate#setProperty(String, String)} implementation.</li>
	 * </ul>
	 * <p><b>Note:</b> further properties may be defined by child classes.</p>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value.
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(K_TAU_PROP_NAME)) {
			k_tau = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			l = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(L_TAU_PROP_NAME)) {
			l_tau = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(DELAY_PROP_NAME)) {
			delay = Integer.parseInt(propertyValue);
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
			//  conditional MI calculator
			condMiCalc.setProperty(propertyName, propertyValue);
			props.put(propertyName, propertyValue);
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			return Integer.toString(k);
		} else if (propertyName.equalsIgnoreCase(K_TAU_PROP_NAME)) {
			return Integer.toString(k_tau);
		} else if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			return Integer.toString(l);
		} else if (propertyName.equalsIgnoreCase(L_TAU_PROP_NAME)) {
			return Integer.toString(l_tau);
		} else if (propertyName.equalsIgnoreCase(DELAY_PROP_NAME)) {
			return Integer.toString(delay);
		} else if (propertyName.equalsIgnoreCase(PROP_AUTO_EMBED_METHOD)) {
			return autoEmbeddingMethod;
		} else if (propertyName.equalsIgnoreCase(PROP_K_SEARCH_MAX)) {
			return Integer.toString(k_search_max);
		} else if (propertyName.equalsIgnoreCase(PROP_TAU_SEARCH_MAX)) {
			return Integer.toString(tau_search_max);
		} else if (propertyName.equalsIgnoreCase(PROP_RAGWITZ_NUM_NNS)) {
			// if we're using a KSG estimator this will have been handled
			//  by the child class. Else it doesn't matter whether it has been
			//  explicitly set or not by the user, we'll return the
			//  current value held here:
			return Integer.toString(ragwitz_num_nns);
		} else {
			// No property matches for this class, assume it is for the underlying
			//  conditional MI calculator
			return condMiCalc.getProperty(propertyName);
		}
	}
	
	@Override
	public void setObservations(double[] source, double[] destination) throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}

	@Override
	public void startAddObservations() {
		vectorOfSourceTimeSeries = new Vector<double[]>();
		vectorOfDestinationTimeSeries = new Vector<double[]>();
		vectorOfValidityOfSource = new Vector<boolean[]>();
		vectorOfValidityOfDestination = new Vector<boolean[]>();
	}
	
	@Override
	public void addObservations(double[] source, double[] destination)
			throws Exception {
		// Store these observations in our vectors for now
		vectorOfSourceTimeSeries.add(source);
		vectorOfDestinationTimeSeries.add(destination);
		vectorOfValidityOfSource.add(null); // All observations were valid
		vectorOfValidityOfDestination.add(null); // All observations were valid
	}

	/**
	 * Adds this set of observations to compute the PDFs from, but
	 * only where these observations are indicated to be valid.
	 * 
	 * @param source time series of observations for the source variable.
	 * @param destination time series of observations for the destination variable.
	 * @param sourceValid array (with indices the same as source) indicating whether the source at that index is valid.
	 * @param destValid array (with indices the same as destination) indicating whether the destination at that index is valid.
	 * @throws Exception
	 */
	public void addObservations(double[] source, double[] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		// Store these observations in our vectors for now
		vectorOfSourceTimeSeries.add(source);
		vectorOfDestinationTimeSeries.add(destination);
		vectorOfValidityOfSource.add(sourceValid); // All observations were valid
		vectorOfValidityOfDestination.add(destValid); // All observations were valid
	}

	/**
	 * Protected method to internally parse and submit observations through
	 *  to the given conditional MI calculator with specific embedding parameter settings
	 *  supplied.
	 * This may be used in the final calculation, or by the auto-embedding 
	 *  procedures, hence the use of static method and arguments rather than
	 *  using any member variables directly.
	 * 
	 * @param condMiCalc_in_use conditional MI calculator to use 
	 * @param k_in_use k embedding dimension to use for target
	 * @param k_tau_in_use target tau embedding delay to use
	 * @param l_in_use l embedding dimension to use for source
	 * @param l_tau_in_use source tau embedding delay to use
	 * @param delay_in_use source-target delay to use
	 * @param source time series of source observations
	 * @param destination time series of destination observations
	 * @return the number of observations added
	 * @throws Exception
	 */
	protected static int addObservationsWithGivenParams(
			ConditionalMutualInfoCalculatorMultiVariate condMiCalc_in_use,
			int k_in_use, int k_tau_in_use, int l_in_use, int l_tau_in_use,
			int delay_in_use,
			double[] source, double[] destination) throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		int startTimeForFirstDestEmbedding_in_use =
				computeStartTimeForFirstDestEmbedding(k_in_use, k_tau_in_use,
						l_in_use, l_tau_in_use, delay_in_use);
		if (source.length < startTimeForFirstDestEmbedding_in_use + 2) {
			// There are no observations to add here, the time series is too short
			// Don't throw an exception, do nothing since more observations
			//  can be added later.
			return 0;
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(destination, k_in_use, k_tau_in_use,
						startTimeForFirstDestEmbedding_in_use,
						destination.length - startTimeForFirstDestEmbedding_in_use - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(destination, 1,
						startTimeForFirstDestEmbedding_in_use + 1,
						destination.length - startTimeForFirstDestEmbedding_in_use - 1);
		double[][] currentSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(source, l_in_use, l_tau_in_use,
						startTimeForFirstDestEmbedding_in_use + 1 - delay_in_use,
						source.length - startTimeForFirstDestEmbedding_in_use - 1);
		condMiCalc_in_use.addObservations(currentSourcePastVectors, currentDestNextVectors, currentDestPastVectors);
		return destination.length - startTimeForFirstDestEmbedding_in_use - 1;
	}

	/**
	 * Protected method to internally parse and submit observations through
	 *  to the given conditional MI calculator with specific embedding parameter settings
	 *  supplied.
	 * This is done given time-series of booleans indicating whether each entry
	 *  is valid
	 * This may be used in the final calculation, or by the auto-embedding 
	 *  procedures, hence the use of static method and arguments rather than
	 *  using any member variables directly.
	 * 
	 * @param condMiCalc_in_use conditional MI calculator to use
	 * @param k_in_use k embedding dimension to use for target
	 * @param k_tau_in_use target tau embedding delay to use
	 * @param l_in_use l embedding dimension to use for source
	 * @param l_tau_in_use source tau embedding delay to use
	 * @param delay_in_use source-target delay to use
	 * @param source time series of source observations
	 * @param destination time series of destination observations
	 * @param source time series of source observations
	 * @param destination time series of destination observations
	 * @param sourceValid array (with indices the same as source) indicating whether
	 * 	the source at that index is valid.
	 * @param destValid array (with indices the same as destination) indicating whether
	 * 	the destination at that index is valid.
	 * @return total number of observations added 
	 * @throws Exception
	 */
	protected static int addObservationsWithGivenParams(
			ConditionalMutualInfoCalculatorMultiVariate condMiCalc_in_use,
			int k_in_use, int k_tau_in_use, int l_in_use, int l_tau_in_use,
			int delay_in_use,
			double[] source, double[] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		
		// Compute the start and end time pairs using our embedding parameters:
		Vector<int[]> startAndEndTimePairs =
				computeStartAndEndTimePairs(k_in_use, k_tau_in_use, l_in_use,
						l_tau_in_use, delay_in_use, sourceValid, destValid);
		
		int totalObservationsAdded = 0;
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			totalObservationsAdded += addObservationsWithGivenParams(
					condMiCalc_in_use, k_in_use, k_tau_in_use, l_in_use,
					l_tau_in_use, delay_in_use,
					MatrixUtils.select(source, startTime, endTime - startTime + 1),
					MatrixUtils.select(destination, startTime, endTime - startTime + 1));
		}
		return totalObservationsAdded;
	}

	@Override
	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length < startTime + numTimeSteps) {
			// There are not enough observations given the arguments here
			throw new Exception("Not enough observations to set here given startTime and numTimeSteps parameters");
		}
		addObservations(MatrixUtils.select(source, startTime, numTimeSteps),
					    MatrixUtils.select(destination, startTime, numTimeSteps));
	}

	/**
	 * Hook in case this or child implementations need to perform any processing on the 
	 *  observation time series prior to their being processed and supplied
	 *  to the underlying MI calculator.
	 * Primarily this is to allow this or child implementation to automatically determine
	 *  embedding parameters if desired.
	 * Child implementations do not need to override this default empty implementation
	 *  if no new functionality is required.
	 */
	protected void preFinaliseAddObservations() throws Exception {
		// Automatically determine the embedding parameters for the given time series
		
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_NONE)) {
			return;
		}
		// Else we need to auto embed
		
		// Check user has set a valid embedding method:
		if (!(autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS_DEST_ONLY) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE))) {
			throw new Exception("Invalid auto-embed method: " + autoEmbeddingMethod);
		}
		
		// Use an AIS calculator to embed both time-series individually:
		ActiveInfoStorageCalculator aisCalc;
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY)) {
			// We're doing Ragwitz auto-embedding:
			// Use a KSG estimator, as Ragwitz is easy with these
			aisCalc = new ActiveInfoStorageCalculatorKraskov();
		} else {
			// We're doing max bias-corrected AIS embedding, grab an instance of 
			//  the corresponding calculator for the estimator type.
			// This is not as nicely object-oriented as I would like, but I can't attach an
			//  to this superclass to grab a relevant AIS calculator without making it abstract,
			//  and I didn't want this class to be abstract.
			if (condMiCalc instanceof ConditionalMutualInfoCalculatorMultiVariateGaussian) {
				aisCalc = new ActiveInfoStorageCalculatorGaussian();
				String numCorrectingSurrogates = getProperty(TransferEntropyCalculatorGaussian.PROP_MAX_CORR_NUM_SURROGATES);
				if (numCorrectingSurrogates != null) {
					aisCalc.setProperty(
							ActiveInfoStorageCalculatorGaussian.PROP_MAX_CORR_AIS_NUM_SURROGATES,
								numCorrectingSurrogates);
				}
			} else if (condMiCalc instanceof ConditionalMutualInfoCalculatorMultiVariateKraskov1) {
				aisCalc = new ActiveInfoStorageCalculatorKraskov(1);
			} else if (condMiCalc instanceof ConditionalMutualInfoCalculatorMultiVariateKraskov2) {
				aisCalc = new ActiveInfoStorageCalculatorKraskov(2);
			// Add these lines in once we have a CMI kernel calculator:
			// } else if (condMiCalc instanceof ConditionalMutualInfoCalculatorMultiVariateKernel) {
			//	aisCalc = new ActiveInfoStorageCalculatorKernel();
			} else {
				throw new RuntimeException("Invalid CMI type found during auto-embedding: " + condMiCalc.getClass().getName());
			}
		}
		// Set the properties for the underlying MI calculator here to match our
		//  properties for our underlying CMI calculator:
		for (String key : props.keySet()) {
			aisCalc.setProperty(key, props.get(key));
		}
		// Set the auto-embedding properties as we require:
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY)) {
			// We're doing Ragwitz auto-embedding
			aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_AUTO_EMBED_METHOD,
					ActiveInfoStorageCalculatorViaMutualInfo.AUTO_EMBED_METHOD_RAGWITZ);
			// In case !ragwitz_num_nns_set and our condMiCalc has a different default number of
			//  kNNs for Kraskov search than miCalc, we had best supply the number directly here:
			aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_RAGWITZ_NUM_NNS,
						getProperty(PROP_RAGWITZ_NUM_NNS));
		} else {
			// We're doing max bias-corrected AIS embedding:
			aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_AUTO_EMBED_METHOD,
					ActiveInfoStorageCalculatorViaMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS);
		}
		aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_K_SEARCH_MAX,
				Integer.toString(k_search_max));
		aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_TAU_SEARCH_MAX,
				Integer.toString(tau_search_max));
		
		// Embed the destination:
		if (debug) {
			System.out.println("Starting embedding of destination:");
		}
		prepareAISCalculator(aisCalc, vectorOfDestinationTimeSeries, vectorOfValidityOfDestination);
		// Set the auto-embedding parameters for the destination:
		k = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.K_PROP_NAME));
		k_tau = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.TAU_PROP_NAME));
		if (debug) {
			System.out.printf("Embedding parameters for destination set to k=%d,k_tau=%d\n",
				k, k_tau);
		}
	
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ) ||
				autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS)) {
			// Embed the source also:
			if (debug) {
				System.out.println("Starting embedding of source:");
			}
			prepareAISCalculator(aisCalc, vectorOfSourceTimeSeries, vectorOfValidityOfSource);
			// Set the auto-embedding parameters for the source:
			l = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.K_PROP_NAME));
			l_tau = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.TAU_PROP_NAME));
			if (debug) {
				System.out.printf("Embedding parameters for source set to l=%d,l_tau=%d\n",
					l, l_tau);
			}

		} else if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE)) {
			if (debug) {
				System.out.println("Starting embedding of source:");
			}

			double bestTE = Double.NEGATIVE_INFINITY;
			int l_candidate_best = 1;
			int l_tau_candidate_best = 1;

			// Iterate over all possible source embeddings
			for (int l_candidate = 1; l_candidate <= k_search_max; l_candidate++) {
				for (int l_tau_candidate = 1; l_tau_candidate <= tau_search_max; l_tau_candidate++) {

					// Use our internal CMI calculator in case it has any particular 
					//  properties we need to have been set already
					prepareCMICalculator(condMiCalc, k, k_tau, l_candidate, l_tau_candidate, delay);
					double thisTE = condMiCalc.computeAverageLocalOfObservations();
					thisTE -= computeAdditionalBiasToRemove();
					
					if (debug) {
						System.out.printf("TE for l=%d, l_tau=%d is %.3f\n",
								l_candidate, l_tau_candidate, thisTE);
					}

					if (thisTE > bestTE) {
						// This parameter setting is the best so far:
						bestTE = thisTE;
						l_candidate_best = l_candidate;
						l_tau_candidate_best = l_tau_candidate;
					}
					if (l_candidate == 1) {
						// tau is irrelevant, so no point testing other values
						break;
					}
				}
			}
			l = l_candidate_best;
			l_tau = l_tau_candidate_best;
			if (debug) {
				System.out.printf("Embedding parameters for source set to l=%d,l_tau=%d\n",
					l, l_tau);
			}
		}

		// Now that embedding parameters are finalised:
		startTimeForFirstDestEmbedding =
				computeStartTimeForFirstDestEmbedding(k, k_tau, l, l_tau, delay);
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

	@Override
	public void finaliseAddObservations() throws Exception {
		// Auto embed if required
		preFinaliseAddObservations();

		separateNumObservations = prepareCMICalculator(condMiCalc, k, k_tau, l, l_tau, delay);
		
		vectorOfSourceTimeSeries = null; // No longer required
		vectorOfDestinationTimeSeries = null; // No longer required
		vectorOfValidityOfSource = null;
		vectorOfValidityOfDestination = null;
	}

	/**
	 * Prepare the given pre-instantiated (and properties supplied)
	 *  Active information storage calculator with the given data set,
	 *  for a calculation of auto-embedding parameters.
	 * 
	 * @param aisCalc_in_use AIS calculator to supply
	 * @param setOfTimeSeriesSamples set of time series samples for the calculation
	 * @param setOfValidities set of time series of validity indications. Each can be a null array if all are valid
	 * @throws Exception
	 */
	protected static void prepareAISCalculator(ActiveInfoStorageCalculator aisCalc,
			Vector<double[]> setOfTimeSeriesSamples, Vector<boolean[]> setOfValidities) 
					throws Exception {
		
		aisCalc.initialise();
		aisCalc.startAddObservations();
		Iterator<boolean[]> validityIterator = setOfValidities.iterator();
		for (double[] timeSeries : setOfTimeSeriesSamples) {
			boolean[] validity = validityIterator.next();
			if (validity == null) {
				aisCalc.addObservations(timeSeries);
			} else {
				aisCalc.addObservations(timeSeries, validity);
			}
		}
		aisCalc.finaliseAddObservations();
	}
	
	/**
	 * Prepare the given pre-instantiated (and properties supplied)
	 *  Conditional mutual information calculator with this data set,
	 *  using the embedding parameters supplied.
	 * This may be used in the final calculation, or by the auto-embedding 
	 *  procedures, hence the use of method arguments rather than
	 *  using the member variables directly.
	 * 
	 * @param condMiCalc_in_use conditional MI calculator to use
	 * @param k_in_use k embedding dimension to use for target
	 * @param k_tau_in_use target tau embedding delay to use
	 * @param l_in_use l embedding dimension to use for source
	 * @param l_tau_in_use source tau embedding delay to use
	 * @param delay_in_use source-target delay to use
	 * @return integer array of number of samples added for each time series pair in the sample set.
	 * @throws Exception
	 */
	protected int[] prepareCMICalculator(
			ConditionalMutualInfoCalculatorMultiVariate condMiCalc_in_use,
			int k_in_use, int k_tau_in_use, int l_in_use, int l_tau_in_use,
			int delay_in_use) throws Exception {

		// Initialise the conditional MI calculator, including any auto-embedding length
		condMiCalc_in_use.initialise(l_in_use, 1, k_in_use);
		condMiCalc_in_use.startAddObservations();
		// Send all of the observations through:
		Iterator<double[]> destIterator = vectorOfDestinationTimeSeries.iterator();
		Iterator<boolean[]> sourceValidityIterator = vectorOfValidityOfSource.iterator();
		Iterator<boolean[]> destValidityIterator = vectorOfValidityOfDestination.iterator();
		int[] separateNumObservationsArray = new int[vectorOfDestinationTimeSeries.size()];
		int setNum = 0;
		for (double[] source : vectorOfSourceTimeSeries) {
			double[] destination = destIterator.next();
			boolean[] sourceValidity = sourceValidityIterator.next();
			boolean[] destValidity = destValidityIterator.next();
			int observationsAddedThisTime = 0;
			if (sourceValidity == null) {
				// Add the whole time-series
				observationsAddedThisTime = addObservationsWithGivenParams(
						condMiCalc_in_use, k_in_use, k_tau_in_use, l_in_use,
						l_tau_in_use, delay_in_use, source, destination);
			} else {
				observationsAddedThisTime = addObservationsWithGivenParams(
						condMiCalc_in_use, k_in_use, k_tau_in_use, l_in_use,
						l_tau_in_use, delay_in_use, source, destination,
						sourceValidity, destValidity);
			}
			separateNumObservationsArray[setNum++] = observationsAddedThisTime;
		}
		
		// TODO do we need to throw an exception if there are no observations to add?
		condMiCalc_in_use.finaliseAddObservations();
		
		return separateNumObservationsArray;
	}
	
	@Override
	public void setObservations(double[] source, double[] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		
		startAddObservations();
		// Add these observations and the indication of their validity
		//  for later analysis:
		vectorOfSourceTimeSeries.add(source);
		vectorOfDestinationTimeSeries.add(destination);
		vectorOfValidityOfSource.add(sourceValid);
		vectorOfValidityOfDestination.add(destValid);
		finaliseAddObservations();
	}

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid series of both source and destinations. (I.e. all points within the
	 *  embedding vectors must be valid, even if the invalid points won't be included
	 *  in any tuples)
	 * 
	 * <p>Made public so it can be used if one wants to compute the number of
	 *  observations prior to setting the observations.</p>
	 * 
	 * @param k_in_use k embedding dimension to use for target
	 * @param k_tau_in_use target tau embedding delay to use
	 * @param l_in_use l embedding dimension to use for source
	 * @param l_tau_in_use source tau embedding delay to use
	 * @param delay_in_use source-target delay to use
	 * @param sourceValid a time series (with indices the same as observations)
	 *  indicating whether the entry in observations at that index is valid for the source; 
	 * @param destValid as described for <code>sourceValid</code>
	 * @return a vector for start and end time pairs of valid series
	 *  of observations.
	 */
	public static Vector<int[]> computeStartAndEndTimePairs(
			int k_in_use, int k_tau_in_use, int l_in_use, int l_tau_in_use,
			int delay_in_use,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		
		if (sourceValid.length != destValid.length) {
			throw new Exception("Validity arrays must be of same length");
		}
		
		int lengthOfDestPastRequired = (k_in_use-1)*k_tau_in_use + 1;
		int lengthOfSourcePastRequired = (l_in_use-1)*l_tau_in_use + 1;
		// int numSourcePointsBeforeDestStart = delay_in_use - 1 + lengthOfSourcePastRequired
		//									- lengthOfDestPastRequired;

		// Scan along the data avoiding invalid values
		int startTime = 0;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();

		// Simple solution -- this takes more complexity in time, but is 
		//  much faster to code:
		boolean previousWasOk = false;
		int startTimeForFirstDestEmbedding_in_use =
				computeStartTimeForFirstDestEmbedding(k_in_use, k_tau_in_use, l_in_use,
						l_tau_in_use, delay_in_use);
		for (int t = startTimeForFirstDestEmbedding_in_use; t < destValid.length - 1; t++) {
			// Check the tuple with the history vector starting from
			//  t and running backwards
			if (previousWasOk) {
				// Just check the very next values of each:
				if (destValid[t + 1] && sourceValid[t + 1 - delay_in_use]) {
					// We can continue adding to this sequence
					continue;
				} else {
					// We need to shut down this sequence now
					previousWasOk = false;
					int[] timePair = new int[2];
					timePair[0] = startTime;
					timePair[1] = t; // Previous time step was last valid one
					startAndEndTimePairs.add(timePair);
					continue;
				}
			}
			// Otherwise we're trying to start a new sequence, so check all values
			if (!destValid[t + 1]) {
				continue;
			}
			boolean allOk = true;
			for (int tBack = 0; tBack < lengthOfDestPastRequired; tBack++) {
				if (!destValid[t - tBack]) {
					allOk = false;
					break;
				}
			}
			if (!allOk) {
				continue;
			}
			allOk = true;
			for (int tBack = delay_in_use - 1; tBack < delay_in_use - 1 + lengthOfSourcePastRequired; tBack++) {
				if (!sourceValid[t - tBack]) {
					allOk = false;
					break;
				}
			}
			if (!allOk) {
				continue;
			}
			// Postcondition: We've got a first valid tuple:
			startTime = t - startTimeForFirstDestEmbedding_in_use;
			previousWasOk = true;
		}
		// Now check if we were running a sequence and terminate it:
		if (previousWasOk) {
			// We need to shut down this sequence now
			previousWasOk = false;
			int[] timePair = new int[2];
			timePair[0] = startTime;
			timePair[1] = destValid.length - 1;
			startAndEndTimePairs.add(timePair);
		}

		return startAndEndTimePairs;
	}
	
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		return condMiCalc.computeAverageLocalOfObservations();
	}

	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] local = condMiCalc.computeLocalOfPreviousObservations();
		if (!condMiCalc.getAddedMoreThanOneObservationSet()) {
			double[] localsToReturn = new double[local.length + startTimeForFirstDestEmbedding + 1];
			System.arraycopy(local, 0, localsToReturn, startTimeForFirstDestEmbedding + 1, local.length);
			return localsToReturn;
		} else {
			return local;
		}
	}

	@Override
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations)
					throws Exception {
		if (newSourceObservations.length != newDestObservations.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					newSourceObservations.length, newDestObservations.length));
		}
		if (newDestObservations.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to compute for here
			return new double[newDestObservations.length];
		}
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
		double[] local = condMiCalc.computeLocalUsingPreviousObservations(
						newSourcePastVectors, newDestNextVectors, newDestPastVectors);
		// Pad the front of the array with zeros where local TE isn't defined:
		double[] localsToReturn = new double[local.length + startTimeForFirstDestEmbedding + 1];
		System.arraycopy(local, 0, localsToReturn, startTimeForFirstDestEmbedding + 1, local.length);
		return localsToReturn;

	}
	
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		// Reorder the source vectors in the surrogates, not the destination
		return condMiCalc.computeSignificance(1, numPermutationsToCheck);
	}

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		// Reorder the source vectors in the surrogates, not the destination
		return condMiCalc.computeSignificance(1, newOrderings);
	}

	@Override
	public double getLastAverage() {
		return condMiCalc.getLastAverage();
	}

	@Override
	public int getNumObservations() throws Exception {
		return condMiCalc.getNumObservations();
	}
	
	/**
	 * Retrieve an array of the number of observations that 
	 * were added by each call to {@link #addObservations(double[], double[])} etc.,
	 * in order of them being called.
	 * The actual number of observations for each call is computed <b>after</b>
	 * any auto-embedding is performed.
	 * 
	 * @return
	 */
	public int[] getSeparateNumObservations() {
		return separateNumObservations;
	}
	
	@Override
	public boolean getAddedMoreThanOneObservationSet() {
		return condMiCalc.getAddedMoreThanOneObservationSet();
	}

	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
		condMiCalc.setDebug(debug);
	}
}
