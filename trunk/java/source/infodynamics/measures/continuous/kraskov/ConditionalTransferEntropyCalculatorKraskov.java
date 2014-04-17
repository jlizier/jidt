package infodynamics.measures.continuous.kraskov;

import java.util.Hashtable;

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.ConditionalTransferEntropyCalculator;
import infodynamics.measures.continuous.ConditionalTransferEntropyCalculatorViaCondMutualInfo;

/**
 * 
 * <p>
 * Implements a conditional transfer entropy calculator using a conditional MI calculator
 * implementing the Kraskov-Grassberger estimator.
 * This is achieved by plugging in a {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}
 * as the calculator into {@link TransferEntropyCalculatorViaCondMutualInfo}.
 * </p>
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct: {@link #ConditionalTransferEntropyCalculatorKraskov()}</li>
 * 		<li>Set properties: {@link #setProperty(String, String)} for each relevant property, including those
 * 			of either {@link ConditionalTransferEntropyCalculatorViaCondMutualInfo#setProperty(String, String)}
 * 			or {@link ConditionalMutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)}.</li>
 *		<li>Initialise: by calling one of {@link #initialise()} etc.</li>
 * 		<li>Add observations to construct the PDFs: {@link #setObservations(double[], double[], double[][])},
 * 			or [{@link #startAddObservations()},
 * 			{@link #addObservations(double[])}*, {@link #finaliseAddObservations()}]
 *   		Note: If not using setObservations(), the results from computeLocal
 *   		will be concatenated directly, and getSignificance will mix up observations 
 *          from separate trials (added in separate {@link #addObservations(double[])} calls.</li> 
 * 		<li>Compute measures: e.g. {@link #computeAverageLocalOfObservations()} or
 * 			{@link #computeLocalOfPreviousObservations()} etc </li>
 * 	</ol>
 * </p>
 * 
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
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
 * @see "Kraskov, A., Stoegbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138;
 *  <a href='http://dx.doi.org/10.1103/PhysRevE.69.066138'>download</a>
 *  (for introduction of Kraskov-Grassberger method for MI)"
 * @see "G. Gomez-Herrero, W. Wu, K. Rutanen, M. C. Soriano, G. Pipa, and R. Vicente,
 *    arXiv:1008.0539, 2010;
 *  <a href='http://arxiv.org/abs/1008.0539'>download</a>
 *  (for introduction of Kraskov-Grassberger technique to transfer entropy)"
 * @see ConditionalTransferEntropyCalculator
 *
 */
public class ConditionalTransferEntropyCalculatorKraskov
	extends ConditionalTransferEntropyCalculatorViaCondMutualInfo {

	public static final String COND_MI_CALCULATOR_KRASKOV1 = ConditionalMutualInfoCalculatorMultiVariateKraskov1.class.getName();
	public static final String COND_MI_CALCULATOR_KRASKOV2 = ConditionalMutualInfoCalculatorMultiVariateKraskov2.class.getName();
	
	/**
	 * Property for setting which underlying Kraskov-Grassberger algorithm to use.
	 * Will only be applied at the next initialisation.
	 */
	public final static String PROP_KRASKOV_ALG_NUM = "ALG_NUM";
	
	protected int kraskovAlgorithmNumber = 1;
	protected boolean algChanged = false;
	/**
	 * Storage for the properties ready to pass onto the underlying conditional MI calculators should they change 
	 */
	protected Hashtable<String,String> props;

	/**
	 * Creates a new instance of the Kraskov-estimate style conditional transfer entropy calculator
	 * 
	 * Uses algorithm 1 by default, as per Gomez-Herro et al.
	 * 
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ConditionalTransferEntropyCalculatorKraskov() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(COND_MI_CALCULATOR_KRASKOV1);
		kraskovAlgorithmNumber = 1;
		props = new Hashtable<String,String>();
	}

	/**
	 * Creates a new instance of the Kraskov-Grassberger style conditional transfer entropy calculator,
	 *  with the supplied conditional MI calculator name
	 * 
	 * @param calculatorName fully qualified name of the underlying MI class.
	 *    Must be {@link #COND_MI_CALCULATOR_KRASKOV1} or {@link #COND_MI_CALCULATOR_KRASKOV2}
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ConditionalTransferEntropyCalculatorKraskov(String calculatorName) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(calculatorName);
		// Now check that it was one of our Kraskov-Grassberger calculators:
		if (calculatorName.equalsIgnoreCase(COND_MI_CALCULATOR_KRASKOV1)) {
			kraskovAlgorithmNumber = 1;
		} else if (calculatorName.equalsIgnoreCase(COND_MI_CALCULATOR_KRASKOV2)) {
			kraskovAlgorithmNumber = 2;
		} else {
			throw new ClassNotFoundException("Must be an underlying Kraskov-Grassberger conditional MI calculator");
		}
		props = new Hashtable<String,String>();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.TransferEntropyCalculatorViaCondMutualInfo#initialise(int, int, int, int, int)
	 */
	@Override
	public void initialise(int k, int k_tau, int l, int l_tau, int delay,
			int[] condEmbedDims, int[] cond_taus, int[] condDelays)
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
		
		super.initialise(k, k_tau, l, l_tau, delay, condEmbedDims, cond_taus, condDelays);
	}

	/**
	 * Sets properties for the calculator.
	 * Valid properties include:
	 * <ul>
	 * 	<li>{@link #PROP_KRASKOV_ALG_NUM} - which Kraskov algorithm number to use (1 or 2). Will only be applied at the next initialisation.</li>
	 *  <li>Any valid properties for {@link ConditionalTransferEntropyCalculatorViaCondMutualInfo#setProperty(String, String)}</li>
	 * 	<li>Any valid properties for {@link ConditionalMutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)}</li>
	 * </ul>
	 * One should set {@link ConditionalMutualInfoCalculatorMultiVariateKraskov#PROP_K} here, the number
	 *  of neighbouring points one should count up to in determining the joint kernel size. 
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property (as a string)
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
			props.put(propertyName, propertyValue); // This will keep properties for the super class as well as the cond MI calculator, but this is ok
		}
	}
}
