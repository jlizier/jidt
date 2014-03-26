package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.measures.continuous.ActiveInfoStorageCalculatorViaMutualInfo;

/**
 * 
 * <p>
 * Implements an active information storage calculator using model of 
 * Gaussian variables with linear interactions.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateGaussian}
 * as the calculator into {@link ActiveInfoStorageCalculatorViaMutualInfo}.
 * </p> 
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct: {@link #ActiveInfoStorageCalculatorGaussian()}</li>
 * 		<li>Set properties: {@link #setProperty(String, String)} for each relevant property, including those
 * 			of either {@link ActiveInfoStorageCalculatorViaMutualInfo#setProperty(String, String)}
 * 			or {@link MutualInfoCalculatorMultiVariateGaussian#setProperty(String, String)}.</li>
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
public class ActiveInfoStorageCalculatorGaussian
	extends ActiveInfoStorageCalculatorViaMutualInfo {
	
	public static final String MI_CALCULATOR_GAUSSIAN = MutualInfoCalculatorMultiVariateGaussian.class.getName();
	
	/**
	 * Creates a new instance of the Gaussian-estimate style active info storage calculator
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorGaussian() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(MI_CALCULATOR_GAUSSIAN);
	}
}
