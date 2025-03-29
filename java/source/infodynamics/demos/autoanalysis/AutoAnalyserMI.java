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

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian;
import infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2;
import infodynamics.measures.discrete.MutualInformationCalculatorDiscrete;

import javax.swing.JOptionPane;
import javax.swing.event.DocumentListener;

import java.awt.event.ActionListener;
import java.awt.event.MouseListener;

/**
 * This class provides a GUI to build a simple mutual information calculation,
 *  and supply the code to execute it.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public class AutoAnalyserMI extends AutoAnalyserChannelCalculator
	implements ActionListener, DocumentListener, MouseListener {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	protected static final String DISCRETE_PROPNAME_TIME_DIFF = "time difference";
	
	protected static final String CALC_TYPE_KRASKOV_ALG1 = CALC_TYPE_KRASKOV + " alg. 1";
	protected static final String CALC_TYPE_KRASKOV_ALG2 = CALC_TYPE_KRASKOV + " alg. 2";
	
	public AutoAnalyserMI() {
		super();
	}

	public AutoAnalyserMI(String pathToAutoAnalyserDir) {
		super(pathToAutoAnalyserDir);
	}

	/**
	 * Constructor to initialise the GUI for MI
	 */
	protected void makeSpecificInitialisations() {
		
		super.makeSpecificInitialisations();
		
		// Set up the properties for MI:
		measureAcronym = "MI";
		appletTitle = "JIDT Mutual Information Auto-Analyser"; 
		
		calcTypes = new String[] {
				CALC_TYPE_DISCRETE, CALC_TYPE_BINNED, CALC_TYPE_GAUSSIAN,
				CALC_TYPE_KRASKOV_ALG1, CALC_TYPE_KRASKOV_ALG2,
				CALC_TYPE_KERNEL};
		unitsForEachCalc = new String[] {"bits", "bits", "nats", "nats", "nats", "bits"};
		
		// Discrete:
		discreteClass = MutualInformationCalculatorDiscrete.class;
		discreteProperties = new String[] {
				DISCRETE_PROPNAME_BASE,
				DISCRETE_PROPNAME_TIME_DIFF
		};
		discretePropertyDefaultValues = new String[] {
				"2",
				"0",
		};
		discretePropertyDescriptions = new String[] {
				"Number of discrete states available for each variable (i.e. 2 for binary)",
				"Time-lag from source to dest to consider MI across; must be >= 0 (0 for standard MI)",
		};
		discretePropertyValueChoices = new String[][] {
				null,
				null
		};
		
		// Continuous:
		abstractContinuousClass = MutualInfoCalculatorMultiVariate.class;
		// Common properties for all continuous calcs:
		commonContPropertyNames = new String[] {
				MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF,
				MutualInfoCalculatorMultiVariate.PROP_SURROGATE_TYPE,
				MutualInfoCalculatorMultiVariate.PROP_DYN_CORR_EXCL_TIME,
		};
		commonContPropertiesFieldNames = new String[] {
				"PROP_TIME_DIFF",
				"PROP_SURROGATE_TYPE",
				"PROP_DYN_CORR_EXCL_TIME"
		};
		commonContPropertyDescriptions = new String[] {
				"Time-lag from source to dest to consider MI across; must be >= 0 (0 for standard MI)",
				"Which strategy type to choose for selecting surrogates, default is " + MutualInfoCalculatorMultiVariate.PROP_SHUFFLE,
				"Dynamic correlation exclusion time or <br/>Theiler window (see Kantz and Schreiber); " +
					"0 (default) means no dynamic exclusion window. Only used for rotated surrogate selection for Gaussian Estimator",
		};
		commonContPropertyValueChoices = new String[][] {
				null,
				MutualInfoCalculatorMultiVariateKraskov.VALID_SURROGATE_TYPES,
				null
		};
		// Gaussian properties:
		gaussianProperties = new String[] {
				MutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION
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
		// Kernel:
		kernelProperties = new String[] {
				MutualInfoCalculatorMultiVariateKernel.KERNEL_WIDTH_PROP_NAME,
				MutualInfoCalculatorMultiVariateKernel.NORMALISE_PROP_NAME,			
		};
		kernelPropertiesFieldNames = new String[] {
				"KERNEL_WIDTH_PROP_NAME",
				"NORMALISE_PROP_NAME"			
		};
		kernelPropertyDescriptions = new String[] {
				"Kernel width to be used in the calculation. <br/>If the property " +
						MutualInfoCalculatorMultiVariateKernel.NORMALISE_PROP_NAME +
						" is set, then this is a number of standard deviations; " +
						"otherwise it is an absolute value.",
				"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not  (recommended)",
		};
		kernelPropertyValueChoices = new String[][] {
				null,
				null,
				{"true", "false"}
		};
		// KSG (Kraskov):
		kraskovProperties = new String[] {
				MutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE,
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE,
				MutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
				MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
				MutualInfoCalculatorMultiVariateKraskov.PROP_USE_GPU,
		};
		kraskovPropertiesFieldNames = new String[] {
				"MutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_K",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS",
				"MutualInfoCalculatorMultiVariateKraskov.PROP_USE_GPU"
		};
		kraskovPropertyDescriptions = new String[] {
				"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not (recommended)",
				"Number of k nearest neighbours to use <br/>in the full joint kernel space in the KSG algorithm",
				"Standard deviation for an amount <br/>of random Gaussian noise to add to each variable, " +
						"to avoid having neighbourhoods with artificially large counts. <br/>" +
						"(\"false\" may be used to indicate \"0\".). The amount is added in after any normalisation.",
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
				{"MAX_NORM", "EUCLIDEAN", "EUCLIDEAN_SQUARED"},
				null,
				{"true", "false"},
		};
	}

	/**
	 * Method to assign and initialise our continuous calculator class
	 */
	@Override
	protected MutualInfoCalculatorMultiVariate assignCalcObjectContinuous(String selectedCalcType) throws Exception {
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
			return new MutualInfoCalculatorMultiVariateGaussian();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV_ALG1)) {
			return new MutualInfoCalculatorMultiVariateKraskov1();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV_ALG2)) {
			return new MutualInfoCalculatorMultiVariateKraskov2();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
			return new MutualInfoCalculatorMultiVariateKernel();
		} else {
			throw new Exception("No recognised continuous calculator selected: " +
					selectedCalcType);
		}

	}

	/**
	 * Method to assign and initialise our discrete calculator class
	 */
	protected DiscreteCalcAndArguments assignCalcObjectDiscrete() throws Exception {
		int timeDiff, base;
		try {
			String timeDiffPropValueStr = propertyValues.get(DISCRETE_PROPNAME_TIME_DIFF);
			timeDiff = Integer.parseInt(timeDiffPropValueStr);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot find a value for property " + DISCRETE_PROPNAME_TIME_DIFF);
			return null;
		}
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
				new MutualInformationCalculatorDiscrete(base, base, timeDiff),
				base,
				base + ", " + base + ", " + timeDiff);
	}

	@Override
	protected String pythonSetObsSuffix() {
		String selectedCalcType = (String)
				calcTypeComboBox.getSelectedItem();
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE) ||
					selectedCalcType.equalsIgnoreCase(CALC_TYPE_BINNED)) {
			return "";
		} else {
			// For the moment we could direct all calls to the 1D arrays version, 
			//  but it is working fine with JPype 0.7; later
			//  when we have 2D inputs we should dynamically detect that and return "2D"
			return "";
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new AutoAnalyserMI();
	}
}
