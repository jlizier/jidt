package infodynamics.measures.continuous.kernel;

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.measures.continuous.ActiveInfoStorageCalculatorViaMutualInfo;

/**
 * 
 * <p>
 * Implements an active information storage calculator using kernel estimation.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateKernel}
 * as the calculator into {@link ActiveInfoStorageCalculatorViaMutualInfo}.
 * </p> 
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct: {@link #ActiveInfoStorageCalculatorKernel()}</li>
 * 		<li>Set properties: {@link #setProperty(String, String)} for each relevant property, including those
 * 			of either {@link ActiveInfoStorageCalculatorViaMutualInfo#setProperty(String, String)}
 * 			or {@link MutualInfoCalculatorMultiVariateKernel#setProperty(String, String)}.</li>
 *		<li>Initialise: by calling one of {@link #initialise()} etc.</li>
 * 		<li>Add observations to construct the PDFs: {@link #setObservations(double[])}, or [{@link #startAddObservations()},
 * 			{@link #addObservations(double[])}*, {@link #finaliseAddObservations()}]
 *   		Note: If not using setObservations(), the results from computeLocal
 *   		will be concatenated directly, and getSignificance will mix up observations 
 *          from separate trials (added in separate {@link #addObservations(double[])} calls.</li> 
 * 		<li>Compute measures: e.g. {@link #computeAverageLocalOfObservations()} or
 * 			{@link #computeLocalOfPreviousObservations()} etc </li>
 * 	</ol>
 * </p>
 * 
 * @author Joseph Lizier
 * @see ActiveInfoStorageCalculator
 * @see ActiveInfoStorageCalculatorCorrelationIntegrals
 *
 */
public class ActiveInfoStorageCalculatorKernel
	extends ActiveInfoStorageCalculatorViaMutualInfo {
	
	public static final String MI_CALCULATOR_KERNEL = MutualInfoCalculatorMultiVariateKernel.class.getName();
		
	/**
	 * Creates a new instance of the kernel-estimate style active info storage calculator
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorKernel() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(MI_CALCULATOR_KERNEL);
	}

	/**
	 * Initialises the calculator with the existing values for k, tau and epsilon
	 * 
	 */
	public void initialise() throws Exception {
		initialise(k, tau, ((MutualInfoCalculatorMultiVariateKernel) miCalc).getKernelWidth());
	}

	/**
	 * Initialises the calculator using existing value for tau
	 * 
	 * @param k history length
	 * @param epsilon kernel width
	 */
	public void initialise(int k, double epsilon) throws Exception {
		initialise(k, tau, epsilon);
	}

	/**
	 * Initialises the calculator
	 * 
	 * @param k history length
	 * @param tau embedding delay
	 * @param epsilon kernel width
	 */
	public void initialise(int k, int tau, double epsilon) throws Exception {
		// Set the property before the calculator is initialised by the super class
		miCalc.setProperty(MutualInfoCalculatorMultiVariateKernel.KERNEL_WIDTH_PROP_NAME, Double.toString(epsilon));
		super.initialise(k, tau);
	}
}
