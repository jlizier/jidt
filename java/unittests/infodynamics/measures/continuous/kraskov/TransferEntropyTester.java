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

import java.util.Arrays;
import java.util.Random;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

public class TransferEntropyTester
	extends infodynamics.measures.continuous.TransferEntropyAbstractTester {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testLocalsAverageCorrectly() throws Exception {
		
		TransferEntropyCalculatorKraskov teCalc =
				new TransferEntropyCalculatorKraskov();
		
		String kraskov_K = "4";
		
		teCalc.setProperty(
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testLocalsAverageCorrectly(teCalc, 100, 1);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		TransferEntropyCalculatorKraskov teCalc =
				new TransferEntropyCalculatorKraskov();
		
		String kraskov_K = "4";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKraskov.PROP_KRASKOV_ALG_NUM,
				"2");
		teCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testComputeSignificanceDoesntAlterAverage(teCalc, 100, 1);
	}

	/**
	 * Test the computed univariate TE
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 1, 1, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=1,
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testUnivariateTEforCoupledVariablesFromFile() throws Exception {
		// Test set 1:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		double[] col0 = MatrixUtils.selectColumn(data, 0);
		double[] col1 = MatrixUtils.selectColumn(data, 1);
		// Need to normalise these ourselves rather than letting the calculator do it -
		//  this ensures the extra values in the time series (e.g. last value in source)
		//  are taken into account, in line with TRENTOOL
		col0 = MatrixUtils.normaliseIntoNewArray(col0);
		col1 = MatrixUtils.normaliseIntoNewArray(col1);
		
		// Use various Kraskov k nearest neighbours parameter
		int kNNs = 4;
		// Expected values from TRENTOOL:
		double expectedFromTRENTOOL0to1 = 0.3058006;
		double expectedFromTRENTOOL1to0 = -0.0029744;
		
		System.out.println("Kraskov TE comparison 1 to TRENTOOL - univariate coupled data 1");
		TransferEntropyCalculatorKraskov teCalc =
				new TransferEntropyCalculatorKraskov();
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K, Integer.toString(kNNs));
		// We already normalised above, and this will do a different
		//  normalisation without taking the extra values in to account if we did it
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE, "false");
		// Need consistency for unit tests:
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0");
		
		teCalc.initialise(1);
		teCalc.setObservations(col0, col1);
		double result = teCalc.computeAverageLocalOfObservations();
		System.out.printf("From 2coupledRandomCols 0->1 expecting %.6f, got %.6f\n",
				expectedFromTRENTOOL0to1, result);
		assertEquals(expectedFromTRENTOOL0to1, result, 0.000001);
		
		teCalc.initialise(1);
		teCalc.setObservations(col1, col0);
		result = teCalc.computeAverageLocalOfObservations();
		assertEquals(expectedFromTRENTOOL1to0, result, 0.000001);
		System.out.printf("From 2coupledRandomCols 1->0 expecting %.6f, got %.6f\n",
				expectedFromTRENTOOL1to0, result);

		assertEquals(99, teCalc.getNumObservations());
	}
	
	/**
	 * Test the computed univariate TE
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data, adding dynamic correlation exclusion
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 1, 1, 1, kraskovK, dynCorrExcl)
	 * with these values ensuring source-dest lag 1, history k=1,
	 * embedding lag 1, and dynamic correlation exclusion window dynCorrExcl
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testUnivariateTEforCoupledVariablesFromFileDynCorrExcl() throws Exception {
		// Test set 1:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		double[] col0 = MatrixUtils.selectColumn(data, 0);
		double[] col1 = MatrixUtils.selectColumn(data, 1);
		// Need to normalise these ourselves rather than letting the calculator do it -
		//  this ensures the extra values in the time series (e.g. last value in source)
		//  are taken into account, in line with TRENTOOL
		col0 = MatrixUtils.normaliseIntoNewArray(col0);
		col1 = MatrixUtils.normaliseIntoNewArray(col1);
		
		// Use various Kraskov k nearest neighbours parameter
		int kNNs = 4;
		// Expected values from TRENTOOL for correlation exclusion window 10:
		int exclWindow = 10;
		double expectedFromTRENTOOL0to1 = 0.2930714;
		double expectedFromTRENTOOL1to0 = -0.0387031;
		
		System.out.println("Kraskov TE comparison 1b to TRENTOOL - univariate coupled data with dynamic correlation exclusion");
		TransferEntropyCalculatorKraskov teCalc =
				new TransferEntropyCalculatorKraskov();
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K, Integer.toString(kNNs));
		// We already normalised above, and this will do a different
		//  normalisation without taking the extra values in to account if we did it
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE, "false");
		// Set dynamic correlation exclusion window:
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME,
				Integer.toString(exclWindow));
		
		teCalc.initialise(1);
		teCalc.setObservations(col0, col1);
		double result = teCalc.computeAverageLocalOfObservations();
		System.out.printf("From 2coupledRandomCols 0->1, Theiler window 10, expecting %.6f, got %.6f\n",
				expectedFromTRENTOOL0to1, result);
		assertEquals(expectedFromTRENTOOL0to1, result, 0.000001);
		
		teCalc.initialise(1);
		teCalc.setObservations(col1, col0);
		result = teCalc.computeAverageLocalOfObservations();
		assertEquals(expectedFromTRENTOOL1to0, result, 0.000001);
		System.out.printf("From 2coupledRandomCols 1->0, Theiler window 10, expecting %.6f, got %.6f\n",
				expectedFromTRENTOOL1to0, result);

		assertEquals(99, teCalc.getNumObservations());
		
		// Change dynamic correlation exclusion window:
		exclWindow = 20;
		expectedFromTRENTOOL0to1 = 0.2995997;
		expectedFromTRENTOOL1to0 = -0.0381608;
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME,
				Integer.toString(exclWindow));
		
		teCalc.initialise(1);
		teCalc.setObservations(col0, col1);
		result = teCalc.computeAverageLocalOfObservations();
		System.out.printf("From 2coupledRandomCols 0->1, Theiler window 20, expecting %.6f, got %.6f\n",
				expectedFromTRENTOOL0to1, result);
		assertEquals(expectedFromTRENTOOL0to1, result, 0.000001);
		
		teCalc.initialise(1);
		teCalc.setObservations(col1, col0);
		result = teCalc.computeAverageLocalOfObservations();
		assertEquals(expectedFromTRENTOOL1to0, result, 0.000001);
		System.out.printf("From 2coupledRandomCols 1->0, Theiler window 20, expecting %.6f, got %.6f\n",
				expectedFromTRENTOOL1to0, result);

		assertEquals(99, teCalc.getNumObservations());
	}
	
	/**
	 * Test the computed univariate TE
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 1, 1, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=1,
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testUnivariateTEforRandomDataFromFile() throws Exception {
		
		// Test set 2:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		double[] col0 = MatrixUtils.selectColumn(data, 0);
		double[] col1 = MatrixUtils.selectColumn(data, 1);
		double[] col2 = MatrixUtils.selectColumn(data, 2);
		// Need to normalise these ourselves rather than letting the calculator do it -
		//  this ensures the extra values in the time series (e.g. last value in source)
		//  are taken into account, in line with TRENTOOL
		col0 = MatrixUtils.normaliseIntoNewArray(col0);
		col1 = MatrixUtils.normaliseIntoNewArray(col1);
		col2 = MatrixUtils.normaliseIntoNewArray(col2);

		// Use various Kraskov k nearest neighbours parameter
		int kNNs = 4;
		// Expected values from TRENTOOL:
		double expectedFromTRENTOOL0to1 = -0.0096556;
		double expectedFromTRENTOOL1to2 = 0.0175389;
		double expectedFromTRENTOOL1to0 = 0.0026367;
		double expectedFromTRENTOOL0to2 = -0.00012474;
		double expectedFromTRENTOOL2to0 = -5.4437e-03;
		
		TransferEntropyCalculatorKraskov teCalc =
				new TransferEntropyCalculatorKraskov();
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K, Integer.toString(kNNs));
		// We already normalised above, and this will do a different
		//  normalisation without taking the extra values in to account if we did it
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE, "false");
		// Need consistency for unit tests:
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0");

		System.out.printf("Kraskov TE comparison 2 to TRENTOOL - univariate random data 1 (col 0->1)");
		teCalc.initialise(1);
		teCalc.setObservations(col0, col1);
		double result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL0to1, result, 0.000001);
				
		System.out.printf("  (col 1->2):");
		teCalc.initialise(1);
		teCalc.setObservations(col1, col2);
		result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL1to2, result, 0.000001);

		System.out.printf("  (col 1->0):");
		teCalc.initialise(1);
		teCalc.setObservations(col1, col0);
		result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL1to0, result, 0.000001);

		System.out.printf("  (col 0->2):");
		teCalc.initialise(1);
		teCalc.setObservations(col0, col2);
		result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL0to2, result, 0.000001);

		System.out.printf("  (col 2->0):");
		teCalc.initialise(1);
		teCalc.setObservations(col2, col0);
		result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL2to0, result, 0.000001);
	}

	/**
	 * Test the computed multivariate TE
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * It's multivariate because we use embedding dimension 2 on both source
	 *  and destination.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 2, 1, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=2,
	 * history embedding dimension l=2 on source as well.
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMultivariateTEforCoupledDataFromFile() throws Exception {
		
		// Test set 3:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/4ColsPairedOneStepNoisyDependence-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		double[] col0 = MatrixUtils.selectColumn(data, 0);
		double[] col1 = MatrixUtils.selectColumn(data, 1);
		double[] col2 = MatrixUtils.selectColumn(data, 2);
		double[] col3 = MatrixUtils.selectColumn(data, 3);
		// Need to normalise these ourselves rather than letting the calculator do it -
		//  this ensures the extra values in the time series (e.g. last value in source)
		//  are taken into account, in line with TRENTOOL
		col0 = MatrixUtils.normaliseIntoNewArray(col0);
		col1 = MatrixUtils.normaliseIntoNewArray(col1);
		col2 = MatrixUtils.normaliseIntoNewArray(col2);
		col3 = MatrixUtils.normaliseIntoNewArray(col3);

		// Use various Kraskov k nearest neighbours parameter
		int kNNs = 4;
		// Expected values from TRENTOOL:
		double expectedFromTRENTOOL0to2 = 0.1400645;
		double expectedFromTRENTOOL2to0 = -0.0181459;
		double expectedFromTRENTOOL1to3 = 0.1639186;
		double expectedFromTRENTOOL3to1 = 0.0036976;
		
		TransferEntropyCalculatorKraskov teCalc =
				new TransferEntropyCalculatorKraskov();
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K, Integer.toString(kNNs));
		// We already normalised above, and this will do a different
		//  normalisation without taking the extra values in to account if we did it
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE, "false");
		// Need consistency for unit tests:
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0");

		System.out.println("Kraskov Cond MI as TE - multivariate coupled data 1, k=2,l=2");
		System.out.println("  (0->2)");
		teCalc.initialise(2, 1, 2, 1, 1);
		teCalc.setObservations(col0, col2);
		double result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL0to2, result, 0.000001);

		System.out.println("  (2->0):");
		teCalc.initialise(2, 1, 2, 1, 1);
		teCalc.setObservations(col2, col0);
		result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL2to0, result, 0.000001);

		System.out.println("  (1->3):");
		teCalc.initialise(2, 1, 2, 1, 1);
		teCalc.setObservations(col1, col3);
		result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL1to3, result, 0.000001);

		System.out.println("  (3->1):");
		teCalc.initialise(2, 1, 2, 1, 1);
		teCalc.setObservations(col3, col1);
		result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL3to1, result, 0.000001);
		
		// -------------
		// And finally, confirm that we get different results for k=1,l=1,
		//  which match TRENTOOL
		double expectedFromTRENTOOL0to1_k1l1 = 0.0072169;
		double expectedFromTRENTOOL1to2_k1l1 = 0.0011738;
		
		System.out.println("  (0->1) but with k=1,l=1:");
		teCalc.initialise(1, 1, 1, 1, 1);
		teCalc.setObservations(col0, col1);
		result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL0to1_k1l1, result, 0.000001);

		// And in reverse
		System.out.println("  (1->2) but with k=1,l=1:");
		teCalc.initialise(1, 1, 1, 1, 1);
		teCalc.setObservations(col1, col2);
		result = teCalc.computeAverageLocalOfObservations();
		System.out.printf(" %.5f\n", result);
		assertEquals(expectedFromTRENTOOL1to2_k1l1, result, 0.000001);
	}
	
	public void testAutoEmbeddingRagwitz() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/SFI-heartRate_breathVol_bloodOx.txt");
		double[][] data = afr.getDouble2DMatrix();
		// Select data points 2350:3550
		data = MatrixUtils.selectRows(data, 2349, 3550-2350+1);

		TransferEntropyCalculatorKraskov teCalc = new TransferEntropyCalculatorKraskov();
		// teCalc.setDebug(true);
		
		// Use one thread to test first:
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS, "1");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ); // Embed both source and target
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_K_SEARCH_MAX, "6");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_TAU_SEARCH_MAX, "4");
		// Need consistency for unit tests:
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0");
		teCalc.initialise();
		teCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1));
		int optimisedK = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME));
		int optimisedKTau = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_TAU_PROP_NAME));
		int optimisedL = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME));
		int optimisedLTau = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_TAU_PROP_NAME));
		double teOptimisedSingleThread = teCalc.computeAverageLocalOfObservations();
		System.out.println("TE was " + teOptimisedSingleThread + " for k=" + optimisedK +
				",k_tau=" + optimisedKTau + ",l=" + optimisedL + ",l_tau=" + optimisedLTau +
				" optimised over kNNs=" +
				teCalc.getProperty(TransferEntropyCalculatorKraskov.PROP_RAGWITZ_NUM_NNS));
		
		// Test that the answer was k=5, k_tau=1, l=2, l_tau=1 for this data set
		//  (I've not checked this anywhere else, just making sure our result stays stable)
		assertEquals(5, optimisedK);
		assertEquals(1, optimisedKTau);
		assertEquals(2, optimisedL);
		assertEquals(1, optimisedLTau);
		// Test that kNNs are equal to that used by the MI calculator when we have not set this
		assertEquals(teCalc.getProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K),
				teCalc.getProperty(TransferEntropyCalculatorKraskov.PROP_RAGWITZ_NUM_NNS));
		
		// Test that we get the same answer by a multi-threaded approach
		// Use one thread to test first:
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS, 
				ConditionalMutualInfoCalculatorMultiVariateKraskov.USE_ALL_THREADS);
		teCalc.initialise();
		teCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1));
		double teOptimisedMultiThread = teCalc.computeAverageLocalOfObservations();
		assertEquals(teOptimisedSingleThread, teOptimisedMultiThread, 0.00000001);
		System.out.println("Answer unchanged by multi-threading");

		// Test that optimisation looks the same if source and target swapped
		teCalc.initialise();
		teCalc.setObservations(MatrixUtils.selectColumn(data, 1), MatrixUtils.selectColumn(data, 0));
		assertEquals(teOptimisedSingleThread, teOptimisedMultiThread, 0.00000001);
		int optimisedSwappedK = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME));
		int optimisedSwappedKTau = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_TAU_PROP_NAME));
		int optimisedSwappedL = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME));
		int optimisedSwappedLTau = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_TAU_PROP_NAME));
		assertEquals(optimisedK, optimisedSwappedL);
		assertEquals(optimisedKTau, optimisedSwappedLTau);
		assertEquals(optimisedL, optimisedSwappedK);
		assertEquals(optimisedLTau, optimisedSwappedKTau);
		System.out.println("Ragwitz auto-embedding for source and destination swaps around when we swap source and target");

		// Test that we can turn optimisation off now and we can hard code the parameters:
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_NONE); // No auto embedding
		teCalc.initialise(2, 2, 2, 2, 2);
		teCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1));
		double teDifferentParams = teCalc.computeAverageLocalOfObservations();
		assertFalse(teDifferentParams == teOptimisedSingleThread);
		assertEquals(2, Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME)));
		assertEquals(2, Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_TAU_PROP_NAME)));
		assertEquals(2, Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME)));
		assertEquals(2, Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_TAU_PROP_NAME)));
		System.out.printf("Auto-embedding goes away when we request no auto embedding (result now  TE(k=%d,k_tau=%d,l=%d,l_tau=%d,u=%d)=%.3f)\n",
				2, 2, 2, 2, 2, teDifferentParams);

		// Test that we get the same answer by setting these parameters
		teCalc = new TransferEntropyCalculatorKraskov();
		// Need consistency for unit tests:
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0");
		teCalc.initialise(optimisedK, optimisedKTau, optimisedL, optimisedLTau, 1);
		teCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1));
		double teManual = teCalc.computeAverageLocalOfObservations();
		assertEquals(teOptimisedSingleThread, teManual, 0.00000001);
		System.out.println("Result stable to hard-coding auto-embedded parameters");
		
		// Test optimising destination only
		teCalc.initialise();
		// auto embed destination only
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY);
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_K_SEARCH_MAX, "5");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_TAU_SEARCH_MAX, "5");
		// Explicitly set the source embedding params
		teCalc.setProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME, "1");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.L_TAU_PROP_NAME, "1");
		teCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1));
		double teOptimisedDestOnly = teCalc.computeAverageLocalOfObservations();
		assertFalse(teOptimisedDestOnly == teOptimisedSingleThread);
		assertEquals(optimisedK, Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME)));
		assertEquals(optimisedKTau, Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_TAU_PROP_NAME)));
		assertEquals(1, Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME)));
		assertEquals(1, Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_TAU_PROP_NAME)));
		System.out.printf("Auto-embedding for dest only does not embed the source (result now  TE(k=%d,k_tau=%d,l=%d,l_tau=%d,u=%d)=%.3f)\n",
				optimisedK, optimisedKTau, 1, 1, 1, teOptimisedDestOnly);

		// Finally, test that we can use a different number of kNNs to the MI calculator
		teCalc.initialise();
		// auto embed source and destination:
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ);
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_RAGWITZ_NUM_NNS, "8");
		teCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1));

		// TODO what happens with method to select start and end times if we 
		// have not selected embedding parameters yet???
		teCalc = new TransferEntropyCalculatorKraskov();
		// Need consistency for unit tests:
		teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ);
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_K_SEARCH_MAX, "5");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_TAU_SEARCH_MAX, "5");
		teCalc.initialise();
		boolean[] validity = new boolean[data.length];
		Arrays.fill(validity, true);
		teCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1),
				validity, validity);
		double teOptimisedWithValidity = teCalc.computeAverageLocalOfObservations();
		assertEquals(teOptimisedSingleThread, teOptimisedWithValidity, 0.00000001);
		System.out.println("Answer unchanged by setting validity");
	}
	
	public void testGetSeparateNumObservations() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/SFI-heartRate_breathVol_bloodOx.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		TransferEntropyCalculatorKraskov teCalc = new TransferEntropyCalculatorKraskov();
		
		teCalc.initialise();
		teCalc.startAddObservations();
		int timeStepsPerCall = 100;
		int calls = 10;
		for (int i = 0; i < calls; i++) {
			// Add more samples
			teCalc.addObservations(MatrixUtils.selectColumn(data, 0, i*timeStepsPerCall, timeStepsPerCall),
					MatrixUtils.selectColumn(data, 1, i*timeStepsPerCall, timeStepsPerCall));
		}
		teCalc.finaliseAddObservations();
		@SuppressWarnings("unused")
		double result = teCalc.computeAverageLocalOfObservations();
		
		// Now we want to check how many observations were added at each call:
		int[] samplesPerCall = teCalc.getSeparateNumObservations();
		assertEquals(calls, samplesPerCall.length);
		for (int i = 0; i < calls; i++) {
			// For k = l = 1, we should have timeStepsPerCall - 1 samples per addObservations() call: 
			assertEquals(timeStepsPerCall - 1, samplesPerCall[i]);
		}
		
		// =====================
		// Now run it again with different k and l and embedding lags, etc:
		teCalc.initialise();
		teCalc.startAddObservations();
		// auto embed destination only
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY);
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_K_SEARCH_MAX, "5");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_TAU_SEARCH_MAX, "5");
		// Explicitly set the source embedding params
		teCalc.setProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME, "1");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.L_TAU_PROP_NAME, "1");
		int[] timeStepsPerCallArray = new int[] {100, 200, 150, 300, 99, 54};
		int startTime = 0;
		for (int i = 0; i < timeStepsPerCallArray.length; i++) {
			// Add more samples
			teCalc.addObservations(MatrixUtils.selectColumn(data, 0, startTime, timeStepsPerCallArray[i]),
					MatrixUtils.selectColumn(data, 1, startTime, timeStepsPerCallArray[i]));
			startTime += timeStepsPerCallArray[i];
		}
		teCalc.finaliseAddObservations();
		result = teCalc.computeAverageLocalOfObservations();
		int optimisedK = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME));
		int optimisedKTau = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_TAU_PROP_NAME));
		int timeOfFirstObservationPerSet = (optimisedK - 1)*optimisedKTau + 1;
		System.out.printf("In testing tracking of observations per addObservations() call" +
				" we have auto-embedding dimension %d and lag %d, timeOfFirstObservationPerSet %d\n",
				optimisedK, optimisedKTau, timeOfFirstObservationPerSet);
		
		// Now we want to check how many observations were added at each call:
		samplesPerCall = teCalc.getSeparateNumObservations();
		assertEquals(timeStepsPerCallArray.length, samplesPerCall.length);
		for (int i = 0; i < timeStepsPerCallArray.length; i++) {
			// For timeOfFirstObservationPerSet, we should have timeStepsPerCall - timeOfFirstObservationPerSet
			//  samples per addObservations() call: 
			assertEquals(timeStepsPerCallArray[i] - timeOfFirstObservationPerSet, samplesPerCall[i]);
		}
	}

  public void testAutoEmbeddingAIS() throws Exception {
    System.out.println("Start AIS autoembedding test.");
    // Generate multivariate data
		RandomGenerator rg = new RandomGenerator();
		double[] source = rg.generateNormalData(5000, 0, 1);
		double[] target = rg.generateNormalData(5000, 0, 1);

    for (int i=3; i < source.length; i++) {
      source[i] = 0.2*source[i-3] + 0.2*source[i-2] + 0.2*source[i-1] + source[i];
      target[i] = 0.2*target[i-3] + 0.2*target[i-2] + 0.2*target[i-1] + target[i];
    }

    int correctK = 3;
    int correctL = 3;

    // Instantiate calculator and set search bounds
		TransferEntropyCalculatorKraskov teCalc = 
				new TransferEntropyCalculatorKraskov();
		teCalc.setProperty("k", "4");
		teCalc.setProperty(teCalc.PROP_K_SEARCH_MAX, "3");
		teCalc.setProperty(teCalc.PROP_TAU_SEARCH_MAX, "1");
		teCalc.setProperty(teCalc.PROP_AUTO_EMBED_METHOD, teCalc.AUTO_EMBED_METHOD_MAX_CORR_AIS);
    teCalc.setDebug(true);

    // Run optimisation
		teCalc.initialise();
		teCalc.setObservations(source, target);
		int optimisedK = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME));
		int optimisedL = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME));

    // Test that answer was correct
		assertEquals(correctK, optimisedK);
		assertEquals(correctL, optimisedL);
  }

  public void testAutoEmbeddingTE() throws Exception {
    System.out.println("Start AIS+TE autoembedding test.");

    // Generate multivariate data (note that source time series has no memory)
		RandomGenerator rg = new RandomGenerator();
		double[] source = rg.generateNormalData(10000, 0, 1);
		double[] target = rg.generateNormalData(10000, 0, 1);

    for (int i=3; i < source.length; i++) {
      target[i] = 0.2*source[i-2] + 0.2*source[i-1] + 
                     0.2*target[i-2] + 0.2*target[i-1] + target[i];
    }

    int correctK = 2;
    int correctL = 2;

    // Instantiate calculator and set search bounds
		TransferEntropyCalculatorKraskov teCalc = 
				new TransferEntropyCalculatorKraskov();
		teCalc.setProperty("k", "4");
		teCalc.setProperty(teCalc.PROP_K_SEARCH_MAX, "2");
		teCalc.setProperty(teCalc.PROP_TAU_SEARCH_MAX, "1");
		teCalc.setProperty(teCalc.PROP_AUTO_EMBED_METHOD, teCalc.AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE);
    teCalc.setDebug(true);

    // Run optimisation
		teCalc.initialise();
		teCalc.setObservations(source, target);
		int optimisedK = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME));
		int optimisedL = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME));

    // Test that answer was correct
		assertEquals(correctK, optimisedK);
		assertEquals(correctL, optimisedL);
  }

  public void testAutoEmbeddingTEWithValidity() throws Exception {
	    System.out.println("Start AIS+TE autoembedding test.");

	    // Generate multivariate data (note that source time series has no memory)
			RandomGenerator rg = new RandomGenerator();
			double[] source = rg.generateNormalData(20000, 0, 1);
			double[] target = rg.generateNormalData(20000, 0, 1);
		boolean[] validity = new boolean[target.length];
		Random random = new Random();
		for (int t = 0; t < target.length; t++) {
			if (random.nextDouble() < 0.05) {
				validity[t] = false;
			} else {
				validity[t] = true;
			}
		}
			

	    for (int i=3; i < source.length; i++) {
	      target[i] = 0.2*source[i-2] + 0.2*source[i-1] + 
	                     0.2*target[i-2] + 0.2*target[i-1] + target[i];
	    }

	    int correctK = 2;
	    int correctL = 2;

	    // Instantiate calculator and set search bounds
			TransferEntropyCalculatorKraskov teCalc = 
					new TransferEntropyCalculatorKraskov();
			teCalc.setProperty("k", "4");
			teCalc.setProperty(teCalc.PROP_K_SEARCH_MAX, "2");
			teCalc.setProperty(teCalc.PROP_TAU_SEARCH_MAX, "1");
			teCalc.setProperty(teCalc.PROP_AUTO_EMBED_METHOD, teCalc.AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE);
	    teCalc.setDebug(true);

	    // Run optimisation
			teCalc.initialise();
			teCalc.setObservations(source, target, validity, validity);
			int optimisedK = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME));
			int optimisedL = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME));

	    // Test that answer was correct
			assertEquals(correctK, optimisedK);
			assertEquals(correctL, optimisedL);
	  }

}
