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

import infodynamics.measures.continuous.EntropyCalculator;
import infodynamics.measures.continuous.InfoMeasureCalculatorContinuous;
import infodynamics.measures.continuous.gaussian.EntropyCalculatorGaussian;
import infodynamics.measures.continuous.kernel.EntropyCalculatorKernel;
import infodynamics.measures.continuous.kozachenko.EntropyCalculatorMultiVariateKozachenko;
import infodynamics.measures.discrete.EntropyCalculatorDiscrete;
import infodynamics.measures.discrete.InfoMeasureCalculatorDiscrete;
import infodynamics.measures.discrete.MutualInformationCalculatorDiscrete;
import infodynamics.utils.MatrixUtils;

import java.util.Vector;

import javax.swing.JOptionPane;

/**
 * This class provides a GUI to build a simple calculation
 *  of entropy,
 *  and supply the code to execute it.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public class AutoAnalyserEntropy extends AutoAnalyser {

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
	protected String[] klProperties;
	protected String[] klPropertiesFieldNames;
	protected String[] klPropertyDescriptions;
	protected String[][] klPropertyValueChoices;

	public AutoAnalyserEntropy() {
		super();
	}

	public AutoAnalyserEntropy(String pathToAutoAnalyserDir) {
		super(pathToAutoAnalyserDir);
	}

	/**
	 * Constructor to initialise the GUI for a channel calculator
	 */
	protected void makeSpecificInitialisations() {
		numVariables = 1;
		variableColNumLabels = new String[] {"Variable"};
		useAllCombosCheckBox = true;
		useStatSigCheckBox = false;
		wordForCombinations = "variables";
		variableRelationshipFormatString = "col_%d";
		disableVariableColTextFieldsForAllCombos = new boolean[] {true};
		indentsForAllCombos = 1;
		
		// Set up the properties for Entropy:
		measureAcronym = "H";
		appletTitle = "JIDT Entropy Auto-Analyser";
		
		calcTypes = new String[] {
				CALC_TYPE_DISCRETE, CALC_TYPE_BINNED, CALC_TYPE_GAUSSIAN,
				CALC_TYPE_KOZ_LEO, CALC_TYPE_KERNEL};
		unitsForEachCalc = new String[] {"bits", "bits", "nats", "nats", "bits"};
		
		// Discrete:
		discreteClass = MutualInformationCalculatorDiscrete.class;
		discreteProperties = new String[] {
				DISCRETE_PROPNAME_BASE
		};
		discretePropertyDefaultValues = new String[] {
				"2"
		};
		discretePropertyDescriptions = new String[] {
				"Number of discrete states available for each variable (i.e. 2 for binary)"
		};
		discretePropertyValueChoices = new String[][] {
				null
		};
		
		// Continuous:
		abstractContinuousClass = EntropyCalculator.class;
		// Common properties for all continuous calcs:
		commonContPropertyNames = new String[] {
				// None
		};
		commonContPropertiesFieldNames = new String[] {
				// None
		};
		commonContPropertyDescriptions = new String[] {
				// None
		};
		commonContPropertyValueChoices = new String[][] {
				// None
		};
		// Gaussian properties:
		gaussianProperties = new String[] {
		};
		gaussianPropertiesFieldNames = new String[] {
		};
		gaussianPropertyDescriptions = new String[] {
		};
		gaussianPropertyValueChoices = new String[][] {	
		};
		// Kernel:
		kernelProperties = new String[] {
				EntropyCalculatorKernel.KERNEL_WIDTH_PROP_NAME,
				EntropyCalculatorKernel.NORMALISE_PROP_NAME,			
		};
		kernelPropertiesFieldNames = new String[] {
				"KERNEL_WIDTH_PROP_NAME",
				"NORMALISE_PROP_NAME"			
		};
		kernelPropertyDescriptions = new String[] {
				"Kernel width to be used in the calculation. <br/>If the property " +
						EntropyCalculatorKernel.NORMALISE_PROP_NAME +
						" is set, then this is a number of standard deviations; " +
						"otherwise it is an absolute value.",
				"(boolean) whether to normalise <br/>the incoming time-series to mean 0, standard deviation 1, or not (default true, recommended)",
		};
		kernelPropertyValueChoices = new String[][] {
				null,
				{"true", "false"}
		};
		// KSG (Kraskov):
		klProperties = new String[] {
				// There is a NUM_DIMENSIONS_PROP_NAME property, but this will
				//  only be relevant if we go to multivariate later.
		};
		klPropertiesFieldNames = new String[] {
		};
		klPropertyDescriptions = new String[] {
		};
		klPropertyValueChoices = new String[][] {
		};
	}

	@Override
	protected void fillOutAllCombinations(Vector<int[]> variableCombinations) {
		// All combinations here means all pairs
		for (int s = 0; s < dataColumns; s++) {
			variableCombinations.add(new int[] {s});
		}
	}

	@Override
	protected String[] setUpLoopsForAllCombos(StringBuffer javaCode,
			StringBuffer pythonCode, StringBuffer matlabCode) {
		// Set up loops in the code:
		// 1. Java code
		javaCode.append("    \n");
		javaCode.append("    // Compute for all variables:\n");
		javaCode.append("    for (int v = 0; v < " + dataColumns +
					"; v++) {\n");
		String javaPrefix = "        ";
		javaCode.append(javaPrefix + "// For each variable:\n");
		// 2. Python code
		pythonCode.append("\n");
		pythonCode.append("# Compute for all variables:\n");
		pythonCode.append("for v in range(" + dataColumns + "):\n");
		String pythonPrefix = "    ";
		pythonCode.append(pythonPrefix+ "# For each variable:\n");
		// 3. Matlab code
		matlabCode.append("\n");
		matlabCode.append("% Compute for all variables:\n");
		matlabCode.append("for v = 1:" + dataColumns + "\n");
		String matlabPrefix = "\t";
		matlabCode.append(matlabPrefix + "% For each variable:\n");
		
		// Return the variables to index each column:
		return new String[] {"v"}; 
	}

	@Override
	protected void finaliseLoopsForAllCombos(StringBuffer javaCode,
			StringBuffer pythonCode, StringBuffer matlabCode) {
		
		// 1. Java code
		javaCode.append("    }\n");
		// 2. Python code
		// Nothing to do
		// 3. Matlab code
		matlabCode.append("end\n");
	}
	
	@Override
	protected String formatStringWithColumnNumbers(String formatStr, int[] columnNumbers) {
		// We format the variable number into the
		//  return string here:
		return String.format(formatStr,
				columnNumbers[0]);
	}
	
	@Override
	protected boolean skipColumnCombo(int[] columnCombo) {
		// No reason to skip any columns here
		return false;
	}
	
	@Override
	protected void setObservations(InfoMeasureCalculatorDiscrete calcDiscrete,
			InfoMeasureCalculatorContinuous calcContinuous,
			int[] columnCombo) throws Exception {
		
		String selectedCalcType = (String)
			calcTypeComboBox.getSelectedItem();
		
		int variableColumn = columnCombo[0];
		
		// Set observations
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
			EntropyCalculatorDiscrete calc = (EntropyCalculatorDiscrete) calcDiscrete;
			calc.addObservations(
					MatrixUtils.selectColumn(dataDiscrete, variableColumn));
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_BINNED)) {
			EntropyCalculatorDiscrete calc = (EntropyCalculatorDiscrete) calcDiscrete;
			calc.addObservations(
					MatrixUtils.discretise(
							MatrixUtils.selectColumn(data, variableColumn),
							// Should be no parse error on the alphabet size by now
							Integer.parseInt(propertyValues.get(DISCRETE_PROPNAME_BASE))));
		} else {
			EntropyCalculator calc = (EntropyCalculator) calcContinuous;
			calc.setObservations(
					MatrixUtils.selectColumn(data, variableColumn));
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
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KOZ_LEO)) {
				calcProperties.classSpecificPropertyNames = klProperties;
				calcProperties.classSpecificPropertiesFieldNames = klPropertiesFieldNames;
				calcProperties.classSpecificPropertyDescriptions = klPropertyDescriptions;
				calcProperties.classSpecificPropertyValueChoices = klPropertyValueChoices;
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

	/**
	 * Method to assign and initialise our continuous calculator class
	 */
	@Override
	protected EntropyCalculator assignCalcObjectContinuous(String selectedCalcType) throws Exception {
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
			return new EntropyCalculatorGaussian();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KOZ_LEO)) {
			return new EntropyCalculatorMultiVariateKozachenko();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
			return new EntropyCalculatorKernel();
		} else {
			throw new Exception("No recognised continuous calculator selected: " +
					selectedCalcType);
		}

	}

	/**
	 * Method to assign and initialise our discrete calculator class
	 */
	protected DiscreteCalcAndArguments assignCalcObjectDiscrete() throws Exception {
		String basePropValueStr;
		try {
			basePropValueStr = propertyValues.get(DISCRETE_PROPNAME_BASE);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot find a value for property " + DISCRETE_PROPNAME_BASE);
			return null;
		}
		int base;
		try {
			base = Integer.parseInt(basePropValueStr);
		} catch (NumberFormatException nfe) {
			JOptionPane.showMessageDialog(this,
					"Cannot parse number for property " + DISCRETE_PROPNAME_BASE
					+ ": " + nfe.getMessage());
			resultsLabel.setText("Cannot parse number for property " +
					DISCRETE_PROPNAME_BASE + ": " + nfe.getMessage());
			return null;
		}
		
		return new DiscreteCalcAndArguments(
				new EntropyCalculatorDiscrete(base),
				base,
				Integer.toString(base));
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new AutoAnalyserEntropy();
	}
}
