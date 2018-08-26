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

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.measures.continuous.InfoMeasureCalculatorContinuous;
import infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov1;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov2;
import infodynamics.measures.discrete.ConditionalMutualInformationCalculatorDiscrete;
import infodynamics.measures.discrete.InfoMeasureCalculatorDiscrete;
import infodynamics.utils.MatrixUtils;

import javax.swing.JOptionPane;
import javax.swing.event.DocumentListener;

import java.awt.event.ActionListener;
import java.awt.event.MouseListener;
import java.util.Vector;

/**
 * This class provides a GUI to build a simple conditional mutual information calculation,
 *  and supply the code to execute it.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public class AutoAnalyserCMI extends AutoAnalyser
	implements ActionListener, DocumentListener, MouseListener {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	// Property names for specific continuous calculators:
	protected String[] gaussianProperties;
	protected String[] gaussianPropertiesFieldNames;
	protected String[] gaussianPropertyDescriptions;
	protected String[][] gaussianPropertyValueChoices;
	protected String[] kraskovProperties;
	protected String[] kraskovPropertiesFieldNames;
	protected String[] kraskovPropertyDescriptions;
	protected String[][] kraskovPropertyValueChoices;

	protected static final String CALC_TYPE_KRASKOV_ALG1 = CALC_TYPE_KRASKOV + " alg. 1";
	protected static final String CALC_TYPE_KRASKOV_ALG2 = CALC_TYPE_KRASKOV + " alg. 2";
	
	public AutoAnalyserCMI() {
		super();
	}

	public AutoAnalyserCMI(String pathToAutoAnalyserDir) {
		super(pathToAutoAnalyserDir);
	}

	/**
	 * Constructor to initialise the GUI for CMI
	 */
	protected void makeSpecificInitialisations() {
		
		numVariables = 3;
		variableColNumLabels = new String[] {"Source", "Destination", "Conditional"};
		useAllCombosCheckBox = true;
		useStatSigCheckBox = true;
		wordForCombinations = "pairs";
		variableRelationshipFormatString = "col_%d -> col_%d | col_%d";
		disableVariableColTextFieldsForAllCombos = new boolean[] {true, true, false};
		indentsForAllCombos = 2;
		
		// Set up the properties for CMI:
		measureAcronym = "CMI";
		appletTitle = "JIDT Conditional MI Auto-Analyser"; 
		
		calcTypes = new String[] {
				CALC_TYPE_DISCRETE, CALC_TYPE_BINNED, CALC_TYPE_GAUSSIAN,
				CALC_TYPE_KRASKOV_ALG1, CALC_TYPE_KRASKOV_ALG2};
				// No kernel calculator defined for CMI (yet, and unlikely to happen)
		unitsForEachCalc = new String[] {"bits", "bits", "nats", "nats", "nats"};
		
		// Discrete:
		discreteClass = ConditionalMutualInformationCalculatorDiscrete.class;
		discreteProperties = new String[] {
				DISCRETE_PROPNAME_BASE
		};
		discretePropertyDefaultValues = new String[] {
				"2"
		};
		discretePropertyDescriptions = new String[] {
				"Number of discrete states available for each variable (i.e. 2 for binary).<br/>" +
				"Can be set individually for each variable -- see code."
		};
		discretePropertyValueChoices = new String[][] {
				null
		};
		
		// Continuous:
		abstractContinuousClass = ConditionalMutualInfoCalculatorMultiVariate.class;
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
				ConditionalMutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION,
		};
		gaussianPropertiesFieldNames = new String[] {
				"PROP_BIAS_CORRECTION"
		};
		gaussianPropertyDescriptions = new String[] {
				"Whether the analytically determined bias (as the mean of the<br/>" +
						"surrogate distribution) will be subtracted from all" +
						"calculated values. Default is false."
		};
		gaussianPropertyValueChoices = new String[][] {
				{"true", "false"}
		};
		// KSG (Kraskov):
		kraskovProperties = new String[] {
				ConditionalMutualInfoMultiVariateCommon.PROP_NORMALISE,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_USE_GPU,
		};
		kraskovPropertiesFieldNames = new String[] {
				"ConditionalMutualInfoMultiVariateCommon.PROP_NORMALISE",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_USE_GPU"
		};
		kraskovPropertyDescriptions = new String[] {
				"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not (recommended)",
				"Number of k nearest neighbours to use <br/>in the full joint kernel space in the KSG algorithm",
				"Standard deviation for an amount <br/>of random Gaussian noise to add to each variable, " +
						"to avoid having neighbourhoods with artificially large counts. <br/>" +
						"(\"false\" may be used to indicate \"0\".). The amount is added in after any normalisation.",
				"Dynamic correlation exclusion time or <br/>Theiler window (see Kantz and Schreiber); " +
						"0 (default) means no dynamic exclusion window",
				"<br/>Norm type to use in KSG algorithm between the points in each marginal space. <br/>Options are: " +
						"\"MAX_NORM\" (default), otherwise \"EUCLIDEAN\" or \"EUCLIDEAN_SQUARED\" (both equivalent here)",
				"Number of parallel threads to use <br/>in computation: an integer > 0 or \"USE_ALL\" " +
						"(default, to indicate to use all available processors)",
				"Whether to enable the GPU module (number of threads then has no bearing); boolean, default false"
		};
		kraskovPropertyValueChoices = new String[][] {
				{"true", "false"},
				null,
				null,
				null,
				{"MAX_NORM", "EUCLIDEAN", "EUCLIDEAN_SQUARED"},
				null,
				{"true", "false"}
		};
	}

	@Override
	protected void fillOutAllCombinations(Vector<int[]> variableCombinations) throws Exception {
		// All combinations here means all pairs of sources and destinations,
		//  with the conditional fixed.
		int conditional = Integer.parseInt(variableColTextFields[2].getText());
		if (conditional >= dataColumns) {
			throw new Exception(String.format("%s column must be between 0 and %d for this data set",
					variableColNumLabels[2], dataColumns-1));
		}
		for (int s = 0; s < dataColumns; s++) {
			for (int d = 0; d < dataColumns; d++) {
				variableCombinations.add(new int[] {s, d, conditional});
			}
		}
	}

	@Override
	protected String[] setUpLoopsForAllCombos(StringBuffer javaCode,
			StringBuffer pythonCode, StringBuffer matlabCode) {
		// Set up loops in the code:
		int conditional = Integer.parseInt(variableColTextFields[2].getText());
		// 1. Java code
		javaCode.append("    \n");
		javaCode.append("    int c = " + conditional + ";\n");
		javaCode.append("    // Compute for all source-destination pairs:\n");
		javaCode.append("    for (int s = 0; s < " + dataColumns +
					"; s++) {\n");
		javaCode.append("        for (int d = 0; d < " + dataColumns +
				"; d++) {\n");
		String javaPrefix = "            ";
		javaCode.append(javaPrefix + "// For each source-dest pair (given conditional):\n");
		javaCode.append(javaPrefix + "if ((s == d) || (s == c) || (d == c)) {\n");
		javaCode.append(javaPrefix + "    continue;\n");
		javaCode.append(javaPrefix + "}\n");
		// 2. Python code
		pythonCode.append("\n");
		pythonCode.append("c = " + conditional + "\n");
		pythonCode.append("# Compute for all pairs:\n");
		pythonCode.append("for s in range(" + dataColumns + "):\n");
		pythonCode.append("    for d in range(" + dataColumns + "):\n");
		String pythonPrefix = "        ";
		pythonCode.append(pythonPrefix+ "# For each source-dest pair (given conditional):\n");
		pythonCode.append(pythonPrefix + "if ((s == d) or (s == c) or (d == c)):\n");
		pythonCode.append(pythonPrefix + "    continue\n");
		// 3. Matlab code
		matlabCode.append("\n");
		matlabCode.append("c = " + (conditional+1) + ";\n");
		matlabCode.append("% Compute for all pairs:\n");
		matlabCode.append("for s = 1:" + dataColumns + "\n");
		matlabCode.append("\tfor d = 1:" + dataColumns + "\n");
		String matlabPrefix = "\t\t";
		matlabCode.append(matlabPrefix + "% For each source-dest pair (given conditional):\n");
		matlabCode.append(matlabPrefix + "if ((s == d) || (s == c) || (d == c))\n");
		matlabCode.append(matlabPrefix + "\tcontinue;\n");
		matlabCode.append(matlabPrefix + "end\n");
		
		// Return the variables to index each column:
		return new String[] {"s", "d", "c"}; 
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
				columnNumbers[0], columnNumbers[1], columnNumbers[2]);
	}
	
	@Override
	protected boolean skipColumnCombo(int[] columnCombo) {
		if ((columnCombo[0] == columnCombo[1]) ||
				(columnCombo[0] == columnCombo[2]) ||
				(columnCombo[1] == columnCombo[2])) {
			// Two columns are the same here,
			//  so don't compute the conditional MI
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
		int condColumn = columnCombo[2];
		
		// Set observations
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
			ConditionalMutualInformationCalculatorDiscrete cmiCalc =
					(ConditionalMutualInformationCalculatorDiscrete) calcDiscrete;
			cmiCalc.addObservations(
					MatrixUtils.selectColumn(dataDiscrete, sourceColumn),
					MatrixUtils.selectColumn(dataDiscrete, destColumn),
					MatrixUtils.selectColumn(dataDiscrete, condColumn));
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_BINNED)) {
			ConditionalMutualInformationCalculatorDiscrete cmiCalc =
					(ConditionalMutualInformationCalculatorDiscrete) calcDiscrete;
			cmiCalc.addObservations(
					MatrixUtils.discretise(
						MatrixUtils.selectColumn(data, sourceColumn),
						Integer.parseInt(propertyValues.get(DISCRETE_PROPNAME_BASE))),
					MatrixUtils.discretise(
						MatrixUtils.selectColumn(data, destColumn),
						Integer.parseInt(propertyValues.get(DISCRETE_PROPNAME_BASE))),
					MatrixUtils.discretise(
							MatrixUtils.selectColumn(data, condColumn),
							Integer.parseInt(propertyValues.get(DISCRETE_PROPNAME_BASE))));
		} else {
			ConditionalMutualInfoCalculatorMultiVariate cmiCalc =
					(ConditionalMutualInfoCalculatorMultiVariate) calcContinuous;
			cmiCalc.setObservations(
					MatrixUtils.selectColumn(data, sourceColumn),
					MatrixUtils.selectColumn(data, destColumn),
					MatrixUtils.selectColumn(data, condColumn));
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
	protected ConditionalMutualInfoCalculatorMultiVariate assignCalcObjectContinuous(String selectedCalcType) throws Exception {
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
			return new ConditionalMutualInfoCalculatorMultiVariateGaussian();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV_ALG1)) {
			return new ConditionalMutualInfoCalculatorMultiVariateKraskov1();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV_ALG2)) {
			return new ConditionalMutualInfoCalculatorMultiVariateKraskov2();
		} else {
			throw new Exception("No recognised continuous calculator selected: " +
					selectedCalcType);
		}

	}

	/**
	 * Method to assign and initialise our discrete calculator class
	 */
	protected DiscreteCalcAndArguments assignCalcObjectDiscrete() throws Exception {
		int base;
		try {
			String basePropValueStr = propertyValues.get(DISCRETE_PROPNAME_BASE);
			base = Integer.parseInt(basePropValueStr);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot find a value for property " + DISCRETE_PROPNAME_BASE);
			return null;
		}
		
		return new DiscreteCalcAndArguments(
				new ConditionalMutualInformationCalculatorDiscrete(base, base, base),
				base,
				base + ", " + base + ", " + base);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new AutoAnalyserCMI();
	}
}
