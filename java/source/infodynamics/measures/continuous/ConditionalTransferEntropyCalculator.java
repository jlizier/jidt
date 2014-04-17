package infodynamics.measures.continuous;

/**
 * <p>This specifies the interface for implementations of
 *  the conditional transfer entropy
 *  and local conditional transfer entropy
 *  (see Lizier et al. PRE, 2008, and Lizier et al., Chaos 2010).
 * </p>
 * 
 * <p>Specifically, this specifies the interface for computing
 * the transfer entropy for <i>continuous</i>-valued variables.</p>
 * 
 * @see "Schreiber, Physical Review Letters 85 (2) pp.461-464 (2000);
 *  <a href='http://dx.doi.org/10.1103/PhysRevLett.85.461'>download</a>
 *  (for definition of transfer entropy)"
 * @see "Lizier, Prokopenko and Zomaya, Physical Review E 77, 026110 (2008);
 * <a href='http://dx.doi.org/10.1103/PhysRevE.77.026110'>download</a>
 *  (for the extension to <i>conditional</i> transfer entropy 
 *  or <i>complete</i> where all other causal sources are conditioned on,
 *  and <i>local</i> transfer entropy)"
 * @see "Lizier, Prokopenko and Zomaya, Chaos 20, 3, 037109 (2010);
 * <a href='http://dx.doi.org/10.1063/1.3486801'>download</a>
 *  (for further clarification on <i>conditional</i> transfer entropy 
 *  or <i>complete</i> where all other causal sources are conditioned on)"
 *  
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 *
 */
public interface ConditionalTransferEntropyCalculator extends ChannelCalculatorCommon {

	/**
	 * Property name to specify the history length k
	 */
	public static final String K_PROP_NAME = "k_HISTORY";	
	/**
	 * Embedding delay for the destination past history vector
	 */
	public static final String K_TAU_PROP_NAME = "k_TAU";
	/**
	 * Embedding length for the source past history vector
	 */
	public static final String L_PROP_NAME = "l_HISTORY";
	/**
	 * Embedding delay for the source past history vector
	 */
	public static final String L_TAU_PROP_NAME = "l_TAU";
	/**
	 * Source-destination delay
	 */
	public static final String DELAY_PROP_NAME = "DELAY";
	/**
	 * Property name for embedding lengths of conditional variables
	 */
	public static final String COND_EMBED_LENGTHS_PROP_NAME = "COND_EMBED_LENGTHS";
	/**
	 * Property name for embedding delays of conditional variables
	 */
	public static final String COND_EMBED_DELAYS_PROP_NAME = "COND_TAUS";
	/**
	 * Property name for conditional-destination delays of conditional variables
	 */
	public static final String COND_DELAYS_PROP_NAME = "COND_DELAYS";

	/**
	 * Initialise the calculator for re-use with new observations.
	 * A new history length k can be supplied here; all other parameters
	 * remain unchanged.
	 * 
	 * @param k history length to be considered.
	 * @throws Exception
	 */
	public void initialise(int k) throws Exception;
	
	/**
	 * Initialise the calculator for a single conditional
	 *  variable, with the given destination,
	 *  source and conditional embedding length, setting all
	 *  embedding delays to 1, and the source-dest and
	 *  conditional-dest delays to 1.
	 * 
	 * @param k Length of destination past history to consider
	 * @param l length of source past history to consider
	 * @param condEmbedDim embedding length for one conditional variable.
	 *  Can be 0 if there are no conditional variables.
	 * @throws Exception
	 */
	public void initialise(int k, int l, int condEmbedDim) throws Exception;

	/**
	 * Initialise the calculator with all required parameters supplied,
	 *  for a single conditional variable.
	 * 
	 * @param k Length of destination past history to consider
	 * @param k_tau embedding delay for the destination variable
	 * @param l length of source past history to consider
	 * @param l_tau embedding delay for the source variable
	 * @param delay time lag between last element of source and destination next value
	 * @param condEmbedDim embedding lengths for one conditional variable.
	 *  Can be 0 if there are no conditional variables.
	 * @param cond_tau embedding delay for the conditional variable.
	 *  Ignored if condEmbedDim == 0.
	 * @param condDelay time lags between last element of the conditional variable
	 *  and destination next value.
	 *  Ignored if condEmbedDim == 0.
	 * @throws Exception for inconsistent arguments, e.g. if array lengths differ between 
	 *  condEmbedDims, cond_taus and condDelays.
	 */
	public void initialise(int k, int k_tau, int l, int l_tau, int delay,
			int condEmbedDim, int cond_tau, int condDelay) throws Exception;
	/**
	 * Initialise the calculator with all required parameters supplied.
	 * 
	 * @param k Length of destination past history to consider
	 * @param k_tau embedding delay for the destination variable
	 * @param l length of source past history to consider
	 * @param l_tau embedding delay for the source variable
	 * @param delay time lag between last element of source and destination next value
	 * @param condEmbedDims array of embedding lengths for each conditional variable.
	 *  Can be an empty array or null if there are no conditional variables.
	 * @param cond_taus array of embedding delays for the conditional variables.
	 *  Must be same length as condEmbedDims array.
	 * @param condDelays array of time lags between last element of each conditional variable
	 *  and destination next value.
	 *  Must be same length as condEmbedDims array.
	 * @throws Exception for inconsistent arguments, e.g. if array lengths differ between 
	 *  condEmbedDims, cond_taus and condDelays.
	 */
	public void initialise(int k, int k_tau, int l, int l_tau, int delay,
			int[] condEmbedDims, int[] cond_taus, int[] condDelays) throws Exception;
	
	/**
	 * <p>Set the given property to the given value.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * These can include:
	 * <ul>
	 * 		<li>{@link #K_PROP_NAME}</li>
	 * 		<li>{@link #K_TAU_PROP_NAME}</li>
	 * 		<li>{@link #L_PROP_NAME}</li>
	 * 		<li>{@link #L_TAU_PROP_NAME}</li>
	 * 		<li>{@link #DELAY_PROP_NAME}</li>
	 * 		<li>{@link #COND_EMBED_LENGTHS_PROP_NAME} -- as a comma separated integer list</li>
	 * 		<li>{@link #COND_EMBED_DELAYS_PROP_NAME} -- as a comma separated integer list</li>
	 * 		<li>{@link #COND_DELAYS_PROP_NAME} -- as a comma separated integer list</li>
	 * </ul>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value
	 */	
	public void setProperty(String propertyName, String propertyValue) throws Exception;
	
	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[], double[][])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param conditionals 2D time series array for the conditional variables
	 *        (first index is time, second index is variable number)
	 * @throws Exception
	 */
	public void setObservations(double[] source, double[] destination, double[][] conditionals) throws Exception;
	
	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[], double[][])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param conditionals 1D time series array for the conditional variables
	 *        -- valid only if the calculator was initialised for a single
	 *        conditional variable.
	 * @throws Exception for example if the calculator was not initialised for
	 *        a single conditional variable
	 */
	public void setObservations(double[] source, double[] destination, double[] conditionals) throws Exception;

	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * <p>Note that the arrays source, destination and conditionals must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param conditionals 2D time series array for the conditional variables
	 *        (first index is time, second index is variable number)
	 * @throws Exception
	 */
	public void addObservations(double[] source, double[] destination,  double[][] conditionals) throws Exception;

	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * <p>Note that the arrays source, destination and conditionals must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param conditionals 1D time series array for the conditional variables
	 *        -- valid only if the calculator was initialised for a single
	 *        conditional variable.
	 * @throws Exception for example if the calculator was not initialised for
	 *        a single conditional variable
	 */
	public void addObservations(double[] source, double[] destination,  double[] conditionals) throws Exception;

	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * <p>Note that the arrays source, destination and conditionals must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param conditionals 2D time series array for the conditional variables
	 *        (first index is time, second index is variable number)
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 * @throws Exception
	 */
	public void addObservations(double[] source, double[] destination,
			double[][] conditionals,
			int startTime, int numTimeSteps) throws Exception ;

	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * <p>Note that the arrays source, destination and conditionals must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param conditionals 1D time series array for the conditional variables
	 *        -- valid only if the calculator was initialised for a single
	 *        conditional variable.
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 * @throws Exception for example if the calculator was not initialised for
	 *        a single conditional variable
	 */
	public void addObservations(double[] source, double[] destination,
			double[] conditionals,
			int startTime, int numTimeSteps) throws Exception ;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * @param conditionals 2D time series array for the conditional variables
	 *        (first index is time, second index is variable number)
	 * @param sourceValid time series (with time indices the same as source)
	 *  indicating whether the source at that point is valid.
	 * @param destValid time series (with time indices the same as destination)
	 *  indicating whether the destination at that point is valid.
	 * @param conditionalsValid 2D time series (with time indices the same as conditionals)
	 *  indicating whether the conditional variables at that point are valid.
	 * @throws Exception
	 */
	public void setObservations(double[] source, double[] destination,
			double[][] conditionals,
			boolean[] sourceValid, boolean[] destValid, boolean[][] conditionalsValid) throws Exception;

	/**
	 * Compute local conditional transfer entropy values for the
	 *  observations in the given parameters,
	 *  using the PDFs computed from the previously supplied method calls.
	 * 
	 * @param newSourceObservations new observations for the source variable
	 * @param newDestObservations new observations for the destination variable
	 * @param newCondObservations new observations for the conditional variables
	 * @return
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations,
			double[][] newCondObservations) throws Exception;

	/**
	 * Compute local conditional transfer entropy values for the
	 *  observations in the given parameters,
	 *  using the PDFs computed from the previously supplied method calls.
	 * 
	 * @param newSourceObservations new observations for the source variable
	 * @param newDestObservations new observations for the destination variable
	 * @param newCondObservations 1D time series array for the conditional variables
	 *        -- valid only if the calculator was initialised for a single
	 *        conditional variable.
	 * @return
	 * @throws Exception for example if the calculator was not initialised for
	 *        a single conditional variable
	 */
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations,
			double[] newCondObservations) throws Exception;
}
