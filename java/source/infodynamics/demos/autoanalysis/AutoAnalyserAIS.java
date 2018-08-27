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

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.measures.continuous.ActiveInfoStorageCalculatorViaMutualInfo;
import infodynamics.measures.continuous.InfoMeasureCalculatorContinuous;
import infodynamics.measures.continuous.gaussian.ActiveInfoStorageCalculatorGaussian;
import infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian;
import infodynamics.measures.continuous.kernel.ActiveInfoStorageCalculatorKernel;
import infodynamics.measures.continuous.kernel.ActiveInfoStorageCalculatorMultiVariateKernel;
import infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel;
import infodynamics.measures.continuous.kraskov.ActiveInfoStorageCalculatorKraskov;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov;
import infodynamics.measures.discrete.ActiveInformationCalculatorDiscrete;
import infodynamics.measures.discrete.InfoMeasureCalculatorDiscrete;
import infodynamics.utils.MatrixUtils;

import java.util.Vector;

import javax.swing.JOptionPane;

/**
 * This class provides a GUI to build a simple calculation
 *  of active information storage,
 *  and supply the code to execute it.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public class AutoAnalyserAIS extends AutoAnalyser {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	protected static final String DISCRETE_PROPNAME_K = "k_HISTORY";

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

	public AutoAnalyserAIS() {
		super();
	}

	public AutoAnalyserAIS(String pathToAutoAnalyserDir) {
		super(pathToAutoAnalyserDir);
	}

	/**
	 * Constructor to initialise the GUI for a channel calculator
	 */
	protected void makeSpecificInitialisations() {
		numVariables = 1;
		variableColNumLabels = new String[] {"Variable"};
		useAllCombosCheckBox = true;
		useStatSigCheckBox = true;
		wordForCombinations = "variables";
		variableRelationshipFormatString = "col_%d";
		disableVariableColTextFieldsForAllCombos = new boolean[] {true};
		indentsForAllCombos = 1;
		
		// Set up the properties for Entropy:
		measureAcronym = "AIS";
		appletTitle = "JIDT Active Information Storage Auto-Analyser";
		
		calcTypes = new String[] {
				CALC_TYPE_DISCRETE, CALC_TYPE_BINNED, CALC_TYPE_GAUSSIAN,
				CALC_TYPE_KRASKOV, CALC_TYPE_KERNEL};
		unitsForEachCalc = new String[] {"bits", "bits", "nats", "nats", "bits"};
		
		// Discrete:
		discreteClass = ActiveInformationCalculatorDiscrete.class;
		discreteProperties = new String[] {
				DISCRETE_PROPNAME_BASE,
				DISCRETE_PROPNAME_K
		};
		discretePropertyDefaultValues = new String[] {
				"2",
				"1"
		};
		discretePropertyDescriptions = new String[] {
				"Number of discrete states available for each variable (i.e. 2 for binary)",
				"History embedding length (k_HISTORY)"
		};
		discretePropertyValueChoices = new String[][] {
				null,
				null
		};
		
		// Continuous:
		abstractContinuousClass = ActiveInfoStorageCalculator.class;
		// Common properties for all continuous calcs:
		commonContPropertyNames = new String[] {
				ActiveInfoStorageCalculator.K_PROP_NAME,
				ActiveInfoStorageCalculator.TAU_PROP_NAME,
				ActiveInfoStorageCalculatorViaMutualInfo.PROP_AUTO_EMBED_METHOD,
				ActiveInfoStorageCalculatorViaMutualInfo.PROP_K_SEARCH_MAX,
				ActiveInfoStorageCalculatorViaMutualInfo.PROP_TAU_SEARCH_MAX,
		};
		commonContPropertiesFieldNames = new String[] {
				"K_PROP_NAME",
				"TAU_PROP_NAME",
				"ActiveInfoStorageCalculatorViaMutualInfo.PROP_AUTO_EMBED_METHOD",
				"ActiveInfoStorageCalculatorViaMutualInfo.PROP_K_SEARCH_MAX",
				"ActiveInfoStorageCalculatorViaMutualInfo.PROP_TAU_SEARCH_MAX",
		};
		commonContPropertyDescriptions = new String[] {
				"History embedding length (k_HISTORY)",
				"History embedding delay (k_TAU)",
				"Method to automatically determine embedding length (k_HISTORY)<br/> and delay (k_TAU) for " +
						"the samples. Default is \"" + ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_NONE +
						"\" meaning values are set manually; other values include: <br/>  -- \"" + ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ +
						"\" for use of the Ragwitz criteria for both source and destination (searching up to \"" + ActiveInfoStorageCalculatorKraskov.PROP_K_SEARCH_MAX +
						"\" and \"" + ActiveInfoStorageCalculatorKraskov.PROP_TAU_SEARCH_MAX + "\"); <br/>  -- \"" + ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_MAX_CORR_AIS +
						"\" for maximising the (bias corrected) Active Info Storage (searching up to \"" + ActiveInfoStorageCalculatorKraskov.PROP_K_SEARCH_MAX +
						"\" and \"" + ActiveInfoStorageCalculatorKraskov.PROP_TAU_SEARCH_MAX + "\"); <br/>Use of values other than \"" + ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_NONE +
						"\" leads to any previous settings for embedding lengths and delays to be overwritten after observations are supplied",
				"Max. embedding length to search to <br/>if auto embedding (as determined by " + ActiveInfoStorageCalculatorKraskov.PROP_AUTO_EMBED_METHOD + ")",
				"Max. embedding delay to search to <br/>if auto embedding (as determined by " + ActiveInfoStorageCalculatorKraskov.PROP_AUTO_EMBED_METHOD + ")",
		};
		commonContPropertyValueChoices = new String[][] {
				null,
				null,
				{ActiveInfoStorageCalculatorViaMutualInfo.AUTO_EMBED_METHOD_NONE,
					ActiveInfoStorageCalculatorViaMutualInfo.AUTO_EMBED_METHOD_RAGWITZ,
					ActiveInfoStorageCalculatorViaMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS},
				null,
				null,
		};
		// Gaussian properties:
		gaussianProperties = new String[] {
				MutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION,
				ActiveInfoStorageCalculatorGaussian.PROP_MAX_CORR_AIS_NUM_SURROGATES
		};
		gaussianPropertiesFieldNames = new String[] {
				"MutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION",
				"ActiveInfoStorageCalculatorGaussian.PROP_MAX_CORR_AIS_NUM_SURROGATES"
		};
		gaussianPropertyDescriptions = new String[] {
				"Whether the analytically determined bias (as the mean of the<br/>" +
						"surrogate distribution) will be subtracted from all" +
						"calculated values. Default is false.",
				"Number of surrogates to use in computing the bias correction<br/>if required for " +
						ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_MAX_CORR_AIS + " auto-embedding method.<br/>" +
						"(default is 0, meaning to use analytic method -- recommended)"
		};
		gaussianPropertyValueChoices = new String[][] {	
				{"true", "false"},
				null
		};
		// Kernel:
		kernelProperties = new String[] {
				MutualInfoCalculatorMultiVariateKernel.KERNEL_WIDTH_PROP_NAME,
				MutualInfoCalculatorMultiVariateKernel.DYN_CORR_EXCL_TIME_NAME,
				MutualInfoCalculatorMultiVariateKernel.NORMALISE_PROP_NAME,
				ActiveInfoStorageCalculatorMultiVariateKernel.PROP_MAX_CORR_AIS_NUM_SURROGATES
		};
		kernelPropertiesFieldNames = new String[] {
				"MutualInfoCalculatorMultiVariateKernel.KERNEL_WIDTH_PROP_NAME",
				"MutualInfoCalculatorMultiVariateKernel.DYN_CORR_EXCL_TIME_NAME",
				"MutualInfoCalculatorMultiVariateKernel.NORMALISE_PROP_NAME",
				"ActiveInfoStorageCalculatorMultiVariateKernel.PROP_MAX_CORR_AIS_NUM_SURROGATES"
		};
		kernelPropertyDescriptions = new String[] {
				"Kernel width to be used in the calculation. <br/>If the property " +
						MutualInfoCalculatorMultiVariateKernel.NORMALISE_PROP_NAME +
						" is set, then this is a number of standard deviations; " +
						"otherwise it is an absolute value.",
				"Dynamic correlation exclusion time or <br/>Theiler window (see Kantz and Schreiber); " +
						"0 (default) means no dynamic exclusion window",
				"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not  (default true, recommended)",
				"Number of surrogates to use in computing the bias correction<br/>if required for " +
						ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_MAX_CORR_AIS + " auto-embedding method.<br/>" +
						"(default is 20)"
		};
		kernelPropertyValueChoices = new String[][] {
				null,
				null,
				{"true", "false"},
				null
		};
		// KSG (Kraskov):
		kraskovProperties = new String[] {
				MutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE,
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE,
				MutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME,
				MutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
				MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
				MutualInfoCalculatorMultiVariateKraskov.PROP_USE_GPU,
				ActiveInfoStorageCalculatorKraskov.PROP_KRASKOV_ALG_NUM,
				ActiveInfoStorageCalculatorKraskov.PROP_RAGWITZ_NUM_NNS,
		};
		kraskovPropertiesFieldNames = new String[] {
				"MutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_K",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_USE_GPU",
				"PROP_KRASKOV_ALG_NUM",
				"PROP_RAGWITZ_NUM_NNS"
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
				"Whether to enable the GPU module (number of threads then has no bearing); boolean, default false",
				"Which KSG algorithm to use (1 or 2)",
				"Number of k nearest neighbours for <br/>Ragwitz auto embedding (if used; defaults to match property \"k\")"
		};
		kraskovPropertyValueChoices = new String[][] {
				{"true", "false"},
				null,
				null,
				null,
				{"MAX_NORM", "EUCLIDEAN", "EUCLIDEAN_SQUARED"},
				null,
				{"true", "false"},
				{"1", "2"},
				null,
		};
	}

	@Override
	protected void fillOutAllCombinations(Vector<int[]> variableCombinations) {
		// All combinations here means all variables
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
			ActiveInformationCalculatorDiscrete calc = (ActiveInformationCalculatorDiscrete) calcDiscrete;
			calc.addObservations(
					MatrixUtils.selectColumn(dataDiscrete, variableColumn));
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_BINNED)) {
			ActiveInformationCalculatorDiscrete calc = (ActiveInformationCalculatorDiscrete) calcDiscrete;
			calc.addObservations(
					MatrixUtils.discretise(
							MatrixUtils.selectColumn(data, variableColumn),
							// Should be no parse error on the alphabet size by now
							Integer.parseInt(propertyValues.get(DISCRETE_PROPNAME_BASE))));
		} else {
			ActiveInfoStorageCalculator calc = (ActiveInfoStorageCalculator) calcContinuous;
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
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV)) {
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

	/**
	 * Method to assign and initialise our continuous calculator class
	 */
	@Override
	protected ActiveInfoStorageCalculator assignCalcObjectContinuous(String selectedCalcType) throws Exception {
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
			return new ActiveInfoStorageCalculatorGaussian();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV)) {
			return new ActiveInfoStorageCalculatorKraskov();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
			return new ActiveInfoStorageCalculatorKernel();
		} else {
			throw new Exception("No recognised continuous calculator selected: " +
					selectedCalcType);
		}

	}

	/**
	 * Method to assign and initialise our discrete calculator class
	 */
	protected DiscreteCalcAndArguments assignCalcObjectDiscrete() throws Exception {
		int base, k;
		try {
			String basePropValueStr = propertyValues.get(DISCRETE_PROPNAME_BASE);
			base = Integer.parseInt(basePropValueStr);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot read a value for property " + DISCRETE_PROPNAME_BASE);
			return null;
		}
		try {
			String kPropValueStr = propertyValues.get(DISCRETE_PROPNAME_K);
			k = Integer.parseInt(kPropValueStr);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot read a value for property " + DISCRETE_PROPNAME_K);
			return null;
		}
		
		return new DiscreteCalcAndArguments(
				new ActiveInformationCalculatorDiscrete(base, k),
				base,
				base + ", " + k);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new AutoAnalyserAIS();
	}
}
