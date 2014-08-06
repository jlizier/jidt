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

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.measures.continuous.ActiveInfoStorageCalculatorViaMutualInfo;

/**
 * 
 * <p>
 * Implements an active information storage calculator using Kraskov-Grassberger estimation.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateKraskov1}
 * or {@link MutualInfoCalculatorMultiVariateKraskov2}
 * as the calculator into {@link ActiveInfoStorageCalculatorViaMutualInfo}.
 * </p> 
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct: {@link #ActiveInfoStorageCalculatorKraskov()}</li>
 * 		<li>Set properties: {@link #setProperty(String, String)} for each relevant property, including those
 * 			of either {@link ActiveInfoStorageCalculatorViaMutualInfo#setProperty(String, String)}
 * 			or {@link MutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)}.</li>
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
public class ActiveInfoStorageCalculatorKraskov
	extends ActiveInfoStorageCalculatorViaMutualInfo {
	
	public static final String MI_CALCULATOR_KRASKOV1 = MutualInfoCalculatorMultiVariateKraskov1.class.getName();
	public static final String MI_CALCULATOR_KRASKOV2 = MutualInfoCalculatorMultiVariateKraskov2.class.getName();
		
	/**
	 * Creates a new instance of the Kraskov-Grassberger style active info storage calculator.
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
	}

	/**
	 * Creates a new instance of the Kraskov-Grassberger style active info storage calculator,
	 *  with the supplied MI calculator name
	 * 
	 * @param calculatorName fully qualified name of the underlying MI class.
	 *    Must be {@link #MI_CALCULATOR_KRASKOV1} or {@link #MI_CALCULATOR_KRASKOV2}
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
	}

	/**
	 * Creates a new instance of the Kraskov-Grassberger style active info storage calculator,
	 *  with the supplied Kraskov-Grassberger MI algorithm number
	 * 
	 * @param algorithm either 1 or 2
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
	}

}
