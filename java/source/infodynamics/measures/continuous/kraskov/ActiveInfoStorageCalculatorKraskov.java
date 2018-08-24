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

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.measures.continuous.ActiveInfoStorageCalculatorViaMutualInfo;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;

/**
 * An Active Information Storage (AIS) calculator (implementing {@link ActiveInfoStorageCalculator})
 * which is affected using a 
 * Kraskov-Stoegbauer-Grassberger (KSG) Mutual Information (MI) calculator
 * ({@link MutualInfoCalculatorMultiVariateKraskov}) to make the calculations.
 * 
 * <p>
 * That is, this class implements an AIS calculator using the KSG nearest-neighbour approach.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateKraskov}
 * as the calculator into the parent class {@link ActiveInfoStorageCalculatorViaMutualInfo}.
 * </p> 
 * 
 * <p>Usage is as per the paradigm outlined for {@link ActiveInfoStorageCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step is either a simple call to {@link #ActiveInfoStorageCalculatorKraskov()},
 *      or else specifies which KSG algorithm to implement via {@link #ActiveInfoStorageCalculatorKraskov(int)}
 *      or {@link #ActiveInfoStorageCalculatorKraskov(String)};</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties for
 *      {@link MutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)}
 *      (except {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}
 *      and those of the super-class as outlined
 *      in {@link ActiveInfoStorageCalculatorViaMutualInfo#setProperty(String, String)}),
 *      including auto-embedding properties.</li>
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
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see ActiveInfoStorageCalculator
 * @see ActiveInfoStorageCalculatorViaMutualInfo
 * @see MutualInfoCalculatorMultiVariateKraskov
 *
 */
public class ActiveInfoStorageCalculatorKraskov
	extends ActiveInfoStorageCalculatorViaMutualInfo {
	
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
	 * Creates a new instance of the Kraskov-Stoegbauer-Grassberger style AIS calculator.
	 * 
	 * Uses algorithm 2 by default.
	 *
	 * 
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorKraskov() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
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
	public ActiveInfoStorageCalculatorKraskov(String calculatorName) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(calculatorName);
		// Now check that it was one of our Kraskov-Grassberger calculators:
		if (!calculatorName.equalsIgnoreCase(MI_CALCULATOR_KRASKOV1) &&
				!calculatorName.equalsIgnoreCase(MI_CALCULATOR_KRASKOV2)) {
			throw new ClassNotFoundException("Must be an underlying Kraskov-Grassberger calculator");
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
	public ActiveInfoStorageCalculatorKraskov(int algorithm) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(algorithm == 1 ? MI_CALCULATOR_KRASKOV1 : MI_CALCULATOR_KRASKOV2);
		if ((algorithm != 1) && (algorithm != 2)) {
			throw new ClassNotFoundException("Algorithm must be 1 or 2");
		}
		props = new Hashtable<String,String>();
	}

	@Override
	public void initialise(int k, int tau) throws Exception {
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
		super.initialise(k, tau);
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
	 * 		<li>Any properties accepted by {@link ActiveInfoStorageCalculatorViaMutualInfo#setProperty(String, String)}</li>
	 * 		<li>Or properties accepted by the underlying
	 * 		{@link MutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)} implementation.</li>
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
			props.put(propertyName, propertyValue); // This will keep properties for the super class as well as the MI calculator, but this is ok
		}
	}

	@Override
	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(PROP_KRASKOV_ALG_NUM)) {
			return Integer.toString(kraskovAlgorithmNumber);
		} else {
			// Assume it was a property for the parent class or underlying conditional MI calculator
			return super.getProperty(propertyName);
		}
	}

}
