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

import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.measures.continuous.TransferEntropyCalculator;
import infodynamics.measures.continuous.TransferEntropyCalculatorViaCondMutualInfo;
import infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian;
import infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian;
import infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov;
import infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete;

import javax.swing.JOptionPane;
import javax.swing.event.DocumentListener;

import java.awt.event.ActionListener;
import java.awt.event.MouseListener;

/**
 * This class provides a GUI to build a simple transfer entropy calculation,
 *  and supply the code to execute it.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public class AutoAnalyserTE extends AutoAnalyserChannelCalculator
	implements ActionListener, DocumentListener, MouseListener {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	protected static final String DISCRETE_PROPNAME_K = "k_HISTORY";
	protected static final String DISCRETE_PROPNAME_K_TAU = "k_TAU";
	protected static final String DISCRETE_PROPNAME_L = "l_HISTORY";
	protected static final String DISCRETE_PROPNAME_L_TAU = "l_TAU";
	protected static final String DISCRETE_PROPNAME_DELAY = "DELAY";
	
	public AutoAnalyserTE() {
		super();
	}

	public AutoAnalyserTE(String pathToAutoAnalyserDir) {
		super(pathToAutoAnalyserDir);
	}

	/**
	 * Constructor to initialise the GUI for TE
	 */
	protected void makeSpecificInitialisations() {
		
		super.makeSpecificInitialisations();
		
		// Set up the properties for TE:
		measureAcronym = "TE";
		appletTitle = "JIDT Transfer Entropy Auto-Analyser"; 
		
		calcTypes = new String[] {
				CALC_TYPE_DISCRETE, CALC_TYPE_BINNED, CALC_TYPE_GAUSSIAN,
				CALC_TYPE_KRASKOV, CALC_TYPE_KERNEL};
		unitsForEachCalc = new String[] {"bits", "bits", "nats", "nats", "bits"};
		
		// Discrete:
		discreteClass = TransferEntropyCalculatorDiscrete.class;
		discreteProperties = new String[] {
				DISCRETE_PROPNAME_BASE,
				DISCRETE_PROPNAME_K,
				DISCRETE_PROPNAME_K_TAU,
				DISCRETE_PROPNAME_L,
				DISCRETE_PROPNAME_L_TAU,
				DISCRETE_PROPNAME_DELAY
		};
		discretePropertyDefaultValues = new String[] {
				"2",
				"1",
				"1",
				"1",
				"1",
				"1"
		};
		discretePropertyDescriptions = new String[] {
				"Number of discrete states available for each variable (i.e. 2 for binary)",
				"Destination history embedding length (k_HISTORY)",
				"Destination history embedding delay (k_TAU)",
				"Source history embedding length (l_HISTORY)",
				"Source history embeding delay (l_TAU)",
				"Delay from source to destination (in time steps)",
		};
		discretePropertyValueChoices = new String[][] {
				null,
				null,
				null,
				null,
				null,
				null
		};

		// Continuous:
		abstractContinuousClass = TransferEntropyCalculator.class;
		// Common properties for all continuous calcs:
		commonContPropertyNames = new String[] {
				TransferEntropyCalculator.K_PROP_NAME
		};
		commonContPropertiesFieldNames = new String[] {
				"K_PROP_NAME"
		};
		commonContPropertyDescriptions = new String[] {
				"Destination history embedding length (k_HISTORY)"
		};
		commonContPropertyValueChoices = new String[][] {
				null
		};
		// Gaussian properties:
		gaussianProperties = new String[] {
				TransferEntropyCalculator.K_TAU_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.L_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.L_TAU_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.DELAY_PROP_NAME, // Not common to Kernel
				ConditionalMutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION,
				TransferEntropyCalculatorViaCondMutualInfo.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorViaCondMutualInfo.PROP_K_SEARCH_MAX,
				TransferEntropyCalculatorViaCondMutualInfo.PROP_TAU_SEARCH_MAX,
				TransferEntropyCalculatorViaCondMutualInfo.PROP_RAGWITZ_NUM_NNS
		};
		gaussianPropertiesFieldNames = new String[] {
				"TransferEntropyCalculator.K_TAU_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.L_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.L_TAU_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.DELAY_PROP_NAME", // Not common to Kernel
				"ConditionalMutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION",
				"TransferEntropyCalculatorViaCondMutualInfo.PROP_AUTO_EMBED_METHOD",
				"TransferEntropyCalculatorViaCondMutualInfo.PROP_K_SEARCH_MAX",
				"TransferEntropyCalculatorViaCondMutualInfo.PROP_TAU_SEARCH_MAX",
				"TransferEntropyCalculatorViaCondMutualInfo.PROP_RAGWITZ_NUM_NNS"
		};
		gaussianPropertyDescriptions = new String[] {
				"Destination history embedding delay (k_TAU)",
				"Source history embedding length (l_HISTORY)",
				"Source history embeding delay (l_TAU)",
				"Delay from source to destination (in time steps)",
				"Whether the analytically determined bias (as the mean of the<br/>" +
						"surrogate distribution) will be subtracted from all" +
						"calculated values. Default is false.",
				"Method to automatically determine embedding lengths (k_HISTORY,l_HISTORY)<br/> and delays (k_TAU, l_TAU) for " +
						"destination and potentially source time-series. Default is \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_NONE +
						"\" meaning values are set manually; other values include: <br/>  -- \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_RAGWITZ +
						"\" for use of the Ragwitz criteria for both source and destination (searching up to \"" + TransferEntropyCalculatorViaCondMutualInfo.PROP_K_SEARCH_MAX +
						"\" and \"" + TransferEntropyCalculatorViaCondMutualInfo.PROP_TAU_SEARCH_MAX + "\"); <br/>  -- \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY +
						"\" for use of the Ragwitz criteria for the destination only. <br/>  -- \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS +
						"\" for maximising the (bias corrected) Active Info Storage for both source and destination (searching up to \"" + TransferEntropyCalculatorViaCondMutualInfo.PROP_K_SEARCH_MAX +
						"\" and \"" + TransferEntropyCalculatorViaCondMutualInfo.PROP_TAU_SEARCH_MAX + "\"); <br/>  -- \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS_DEST_ONLY +
						"\" for maximising the (bias corrected) Active Info Storage for the destination only; <br/> -- \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE +
						"\" for maximising (bias corrected) Active Info Storage for the destination first, and then selecting source embedding to maximise (bias-corrected) transfer entropy." +
						"<br/>Use of values other than \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_NONE +
						"\" leads to any previous settings for embedding lengths and delays for the destination and perhaps source to be overwritten after observations are supplied",
				"Max. embedding length to search to <br/>if auto embedding (as determined by " + TransferEntropyCalculatorViaCondMutualInfo.PROP_AUTO_EMBED_METHOD + ")",
				"Max. embedding delay to search to <br/>if auto embedding (as determined by " + TransferEntropyCalculatorViaCondMutualInfo.PROP_AUTO_EMBED_METHOD + ")",
				"Number of k nearest neighbours for <br/>Ragwitz auto embedding (if used; defaults to match property \"k\")"
		};
		gaussianPropertyValueChoices = new String[][] {	
				null,
				null,
				null,
				null,
				{"true", "false"},
				{TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_NONE,
					TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_RAGWITZ,
					TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY,
					TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS,
					TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS_DEST_ONLY,
					TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE},
				null,
				null,
				null,
		};
		// Kernel:
		kernelProperties = new String[] {
				TransferEntropyCalculatorKernel.KERNEL_WIDTH_PROP_NAME,
				TransferEntropyCalculatorKernel.DYN_CORR_EXCL_TIME_NAME,
				TransferEntropyCalculatorKernel.NORMALISE_PROP_NAME,			
		};
		kernelPropertiesFieldNames = new String[] {
				"KERNEL_WIDTH_PROP_NAME",
				"DYN_CORR_EXCL_TIME_NAME",
				"NORMALISE_PROP_NAME"			
		};
		kernelPropertyDescriptions = new String[] {
				"Kernel width to be used in the calculation. <br/>If the property " +
						TransferEntropyCalculatorKernel.NORMALISE_PROP_NAME +
						" is set, then this is a number of standard deviations; " +
						"otherwise it is an absolute value.",
				"Dynamic correlation exclusion time or <br/>Theiler window (see Kantz and Schreiber); " +
						"0 (default) means no dynamic exclusion window",
				"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not  (default true, recommended)",
		};
		kernelPropertyValueChoices = new String[][] {
				null,
				null,
				{"true", "false"}
		};
		// KSG (Kraskov):
		kraskovProperties = new String[] {
				TransferEntropyCalculator.K_TAU_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.L_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.L_TAU_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.DELAY_PROP_NAME, // Not common to Kernel
				ConditionalMutualInfoMultiVariateCommon.PROP_NORMALISE,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_USE_GPU,
				TransferEntropyCalculatorKraskov.PROP_KRASKOV_ALG_NUM,
				TransferEntropyCalculatorViaCondMutualInfo.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorViaCondMutualInfo.PROP_K_SEARCH_MAX,
				TransferEntropyCalculatorViaCondMutualInfo.PROP_TAU_SEARCH_MAX,
				TransferEntropyCalculatorViaCondMutualInfo.PROP_RAGWITZ_NUM_NNS
		};
		kraskovPropertiesFieldNames = new String[] {
				"TransferEntropyCalculator.K_TAU_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.L_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.L_TAU_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.DELAY_PROP_NAME", // Not common to Kernel
				"ConditionalMutualInfoMultiVariateCommon.PROP_NORMALISE",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_USE_GPU",
				"PROP_KRASKOV_ALG_NUM",
				"TransferEntropyCalculatorViaCondMutualInfo.PROP_AUTO_EMBED_METHOD",
				"TransferEntropyCalculatorViaCondMutualInfo.PROP_K_SEARCH_MAX",
				"TransferEntropyCalculatorViaCondMutualInfo.PROP_TAU_SEARCH_MAX",
				"TransferEntropyCalculatorViaCondMutualInfo.PROP_RAGWITZ_NUM_NNS"
		};
		kraskovPropertyDescriptions = new String[] {
				"Destination history embedding delay (k_TAU)",
				"Source history embedding length (l)",
				"Source history embeding delay (l_TAU)",
				"Delay from source to destination (in time steps)",
				"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not (default true, recommended)",
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
				"Method to automatically determine embedding lengths (k_HISTORY,l_HISTORY)<br/> and delays (k_TAU, l_TAU) for " +
						"destination and potentially source time-series. Default is \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_NONE +
						"\" meaning values are set manually; other values include: <br/>  -- \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_RAGWITZ +
						"\" for use of the Ragwitz criteria for both source and destination (searching up to \"" + TransferEntropyCalculatorViaCondMutualInfo.PROP_K_SEARCH_MAX +
						"\" and \"" + TransferEntropyCalculatorViaCondMutualInfo.PROP_TAU_SEARCH_MAX + "\"); <br/>  -- \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY +
						"\" for use of the Ragwitz criteria for the destination only. <br/>  -- \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS +
						"\" for maximising the (bias corrected) Active Info Storage for both source and destination (searching up to \"" + TransferEntropyCalculatorViaCondMutualInfo.PROP_K_SEARCH_MAX +
						"\" and \"" + TransferEntropyCalculatorViaCondMutualInfo.PROP_TAU_SEARCH_MAX + "\"); <br/>  -- \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS_DEST_ONLY +
						"\" for maximising the (bias corrected) Active Info Storage for the destination only; <br/> -- \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE +
						"\" for maximising (bias corrected) Active Info Storage for the destination first, and then selecting source embedding to maximise (bias-corrected) transfer entropy." +
						"<br/>Use of values other than \"" + TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_NONE +
						"\" leads to any previous settings for embedding lengths and delays for the destination and perhaps source to be overwritten after observations are supplied",
				"Max. embedding length to search to <br/>if auto embedding (as determined by " + TransferEntropyCalculatorViaCondMutualInfo.PROP_AUTO_EMBED_METHOD + ")",
				"Max. embedding delay to search to <br/>if auto embedding (as determined by " + TransferEntropyCalculatorViaCondMutualInfo.PROP_AUTO_EMBED_METHOD + ")",
				"Number of k nearest neighbours for <br/>Ragwitz auto embedding (if used; defaults to match property \"k\")"
		};
		kraskovPropertyValueChoices = new String[][] {
				null,
				null,
				null,
				null,
				{"true", "false"},
				null,
				null,
				null,
				{"MAX_NORM", "EUCLIDEAN", "EUCLIDEAN_SQUARED"},
				null,
				{"true", "false"},
				{"1", "2"},
				{TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_NONE,
					TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_RAGWITZ,
					TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY,
					TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS,
					TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS_DEST_ONLY,
					TransferEntropyCalculatorViaCondMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE},
				null,
				null,
				null,
		};
	}

	/**
	 * Method to assign and initialise our continuous calculator class
	 */
	protected TransferEntropyCalculator assignCalcObjectContinuous(String selectedCalcType) throws Exception {
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
			return new TransferEntropyCalculatorGaussian();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV)) {
			return new TransferEntropyCalculatorKraskov();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
			return new TransferEntropyCalculatorKernel();
		} else {
			throw new Exception("No recognised continuous calculator selected: " +
					selectedCalcType);
		}

	}

	/**
	 * Method to assign and initialise our discrete calculator class
	 */
	protected DiscreteCalcAndArguments assignCalcObjectDiscrete() throws Exception {
		int base, k, k_tau, l, l_tau, delay;
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
		try {
			String kTauPropValueStr = propertyValues.get(DISCRETE_PROPNAME_K_TAU);
			k_tau = Integer.parseInt(kTauPropValueStr);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot read a value for property " + DISCRETE_PROPNAME_K_TAU);
			return null;
		}
		try {
			String lPropValueStr = propertyValues.get(DISCRETE_PROPNAME_L);
			l = Integer.parseInt(lPropValueStr);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot read a value for property " + DISCRETE_PROPNAME_L);
			return null;
		}
		try {
			String lTauPropValueStr = propertyValues.get(DISCRETE_PROPNAME_L_TAU);
			l_tau = Integer.parseInt(lTauPropValueStr);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot read a value for property " + DISCRETE_PROPNAME_L_TAU);
			return null;
		}
		try {
			String delayPropValueStr = propertyValues.get(DISCRETE_PROPNAME_DELAY);
			delay = Integer.parseInt(delayPropValueStr);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot read a value for property " + DISCRETE_PROPNAME_DELAY);
			return null;
		}
		
		return new DiscreteCalcAndArguments(
				new TransferEntropyCalculatorDiscrete(base, k, k_tau, l, l_tau, delay),
				base,
				base + ", " + k + ", " + k_tau + ", " + l + ", " + l_tau + ", " + delay);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new AutoAnalyserTE();
	}
}
