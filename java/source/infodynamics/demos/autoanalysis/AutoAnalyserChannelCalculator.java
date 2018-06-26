/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2015, Joseph T. Lizier
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

package infodynamics.demos.autoanalysis;

import infodynamics.measures.continuous.ChannelCalculator;
import infodynamics.measures.continuous.InfoMeasureCalculatorContinuous;
import infodynamics.measures.discrete.ChannelCalculatorDiscrete;
import infodynamics.measures.discrete.InfoMeasureCalculatorDiscrete;
import infodynamics.utils.MatrixUtils;

import java.util.Vector;

/**
 * This abstract class provides a GUI to build a simple calculation
 *  over a source-target channel,
 *  and supply the code to execute it.
 * Child classes fix this to a TE or MI calculation.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public abstract class AutoAnalyserChannelCalculator extends AutoAnalyser {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	// Property names for specific continuous calculators:
	protected String[] gaussianProperties;
	protected String[] gaussianPropertiesFieldNames;
	protected String[] gaussianPropertyDescriptions;
	protected String[][] gaussianPropertyValueChoices;
	protected String[] kernelProperties;
	protected String[] kernelPropertiesFieldNames;
	protected String[] kernelPropertyDescriptions;
	protected String[][] kernelPropertyValueChoices;
	protected String[] kraskovProperties;
	protected String[] kraskovPropertiesFieldNames;
	protected String[] kraskovPropertyDescriptions;
	protected String[][] kraskovPropertyValueChoices;

	
	public AutoAnalyserChannelCalculator() {
		super();
	}

	public AutoAnalyserChannelCalculator(String pathToAutoAnalyserDir) {
		super(pathToAutoAnalyserDir);
	}

	/**
	 * Constructor to initialise the GUI for a channel calculator
	 */
	protected void makeSpecificInitialisations() {
		numVariables = 2;
		variableColNumLabels = new String[] {"Source", "Destination"};
		useAllCombosCheckBox = true;
		useStatSigCheckBox = true;
		wordForCombinations = "pairs";
		variableRelationshipFormatString = "col_%d -> col_%d";
		
		disableVariableColTextFieldsForAllCombos = new boolean[] {true, true};
		indentsForAllCombos = 2;
	}

	@Override
	protected void fillOutAllCombinations(Vector<int[]> variableCombinations) {
		// All combinations here means all pairs
		for (int s = 0; s < dataColumns; s++) {
			for (int d = 0; d < dataColumns; d++) {
				variableCombinations.add(new int[] {s, d});
			}
		}
	}

	@Override
	protected String[] setUpLoopsForAllCombos(StringBuffer javaCode,
			StringBuffer pythonCode, StringBuffer matlabCode) {
		// Set up loops in the code:
		// 1. Java code
		javaCode.append("    \n");
		javaCode.append("    // Compute for all pairs:\n");
		javaCode.append("    for (int s = 0; s < " + dataColumns +
					"; s++) {\n");
		javaCode.append("        for (int d = 0; d < " + dataColumns +
				"; d++) {\n");
		String javaPrefix = "            ";
		javaCode.append(javaPrefix + "// For each source-dest pair:\n");
		javaCode.append(javaPrefix + "if (s == d) {\n");
		javaCode.append(javaPrefix + "    continue;\n");
		javaCode.append(javaPrefix + "}\n");
		// 2. Python code
		pythonCode.append("\n");
		pythonCode.append("# Compute for all pairs:\n");
		pythonCode.append("for s in range(" + dataColumns + "):\n");
		pythonCode.append("    for d in range(" + dataColumns + "):\n");
		String pythonPrefix = "        ";
		pythonCode.append(pythonPrefix+ "# For each source-dest pair:\n");
		pythonCode.append(pythonPrefix + "if (s == d):\n");
		pythonCode.append(pythonPrefix + "    continue\n");
		// 3. Matlab code
		matlabCode.append("\n");
		matlabCode.append("% Compute for all pairs:\n");
		matlabCode.append("for s = 1:" + dataColumns + "\n");
		matlabCode.append("\tfor d = 1:" + dataColumns + "\n");
		String matlabPrefix = "\t\t";
		matlabCode.append(matlabPrefix + "% For each source-dest pair:\n");
		matlabCode.append(matlabPrefix + "if (s == d)\n");
		matlabCode.append(matlabPrefix + "\tcontinue;\n");
		matlabCode.append(matlabPrefix + "end\n");
		
		// Return the variables to index each column:
		return new String[] {"s", "d"}; 
	}

	@Override
	protected void finaliseLoopsForAllCombos(StringBuffer javaCode,
			StringBuffer pythonCode, StringBuffer matlabCode) {
		
		// 1. Java code
		javaCode.append("        }\n");
		javaCode.append("    }\n");
		// 2. Python code
		// Nothing to do
		// 3. Matlab code
		matlabCode.append("\tend\n");
		matlabCode.append("end\n");
	}
	
	@Override
	protected String formatStringWithColumnNumbers(String formatStr, int[] columnNumbers) {
		// We format the source and target variable numbers into the
		//  return string here:
		return String.format(formatStr,
				columnNumbers[0], columnNumbers[1]);
	}
	
	@Override
	protected boolean skipColumnCombo(int[] columnCombo) {
		if (columnCombo[0] == columnCombo[1]) {
			// source and destination columns are the same here,
			//  so don't compute the channel measure
			return true;
		}
		return false;
	}
	
	@Override
	protected void setObservations(InfoMeasureCalculatorDiscrete calcDiscrete,
			InfoMeasureCalculatorContinuous calcContinuous,
			int[] columnCombo) throws Exception {
		
		String selectedCalcType = (String)
			calcTypeComboBox.getSelectedItem();
		
		int sourceColumn = columnCombo[0];
		int destColumn = columnCombo[1];
		
		// Set observations
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
			ChannelCalculatorDiscrete cc = (ChannelCalculatorDiscrete) calcDiscrete;
			cc.addObservations(
					MatrixUtils.selectColumn(dataDiscrete, sourceColumn),
					MatrixUtils.selectColumn(dataDiscrete, destColumn));
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_BINNED)) {
			ChannelCalculatorDiscrete cc = (ChannelCalculatorDiscrete) calcDiscrete;
			cc.addObservations(
					MatrixUtils.discretise(
						MatrixUtils.selectColumn(data, sourceColumn),
						Integer.parseInt(propertyValues.get(DISCRETE_PROPNAME_BASE))),
					MatrixUtils.discretise(
						MatrixUtils.selectColumn(data, destColumn),
						Integer.parseInt(propertyValues.get(DISCRETE_PROPNAME_BASE))));
		} else {
			ChannelCalculator cc = (ChannelCalculator) calcContinuous;
			cc.setObservations(
					MatrixUtils.selectColumn(data, sourceColumn),
					MatrixUtils.selectColumn(data, destColumn));
		}
	}
	
	protected CalcProperties assignCalcProperties(String selectedCalcType)
			throws Exception {
		// Let the super class handle discrete calculators
		CalcProperties calcProperties = super.assignCalcProperties(selectedCalcType);
		if (calcProperties == null) {
			// We need to assign properties for a continuous calculator
			calcProperties = new CalcProperties();
			calcProperties.calc = assignCalcObjectContinuous(selectedCalcType);
			calcProperties.calcClass = calcProperties.calc.getClass();
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
				calcProperties.classSpecificPropertyNames = gaussianProperties;
				calcProperties.classSpecificPropertiesFieldNames = gaussianPropertiesFieldNames;
				calcProperties.classSpecificPropertyDescriptions = gaussianPropertyDescriptions;
				calcProperties.classSpecificPropertyValueChoices = gaussianPropertyValueChoices;
			} else if (selectedCalcType.startsWith(CALC_TYPE_KRASKOV)) {
				// The if statement will work for both MI Kraskov calculators
				calcProperties.classSpecificPropertyNames = kraskovProperties;
				calcProperties.classSpecificPropertiesFieldNames = kraskovPropertiesFieldNames;
				calcProperties.classSpecificPropertyDescriptions = kraskovPropertyDescriptions;
				calcProperties.classSpecificPropertyValueChoices = kraskovPropertyValueChoices;
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
				calcProperties.classSpecificPropertyNames = kernelProperties;
				calcProperties.classSpecificPropertiesFieldNames = kernelPropertiesFieldNames;
				calcProperties.classSpecificPropertyDescriptions = kernelPropertyDescriptions;
				calcProperties.classSpecificPropertyValueChoices = kernelPropertyValueChoices;
			} else {
				calcProperties = null;
				throw new Exception("No recognised calculator selected: " +
						selectedCalcType);
			}
		}
		return calcProperties;
	}

}
