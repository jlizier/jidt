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

package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateAbstractTester;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

public class ConditionalMutualInfoMultiVariateTester extends
		ConditionalMutualInfoMultiVariateAbstractTester {

	public void testLocalsAverageCorrectly() throws Exception {
		ConditionalMutualInfoCalculatorMultiVariateGaussian condMiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateGaussian();
		super.testLocalsAverageCorrectly(condMiCalc, 2, 100);
	}

	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		ConditionalMutualInfoCalculatorMultiVariateGaussian condMiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateGaussian();
		super.testComputeSignificanceDoesntAlterAverage(condMiCalc, 2, 100);
	}
	
	/**
	 * Test the construction of the joint covariance matrix against a known
	 *  test case verified using covariance calculations in matlab/octave
	 */
	public void testJointCovariance() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/4ColsPairedOneStepNoisyDependence-1.txt");
		double[][] data = afr.getDouble2DMatrix();

		//============================
		// Case 1: Autocovariance conditioned on other variables:
		
		double[][] source = MatrixUtils.selectRowsAndColumns(data,
				MatrixUtils.range(0, data.length-2), new int[] {0});
		double[][] dest = MatrixUtils.selectRowsAndColumns(data,
				MatrixUtils.range(1, data.length-1), new int[] {0});
		double[][] others = MatrixUtils.selectRowsAndColumns(data,
				MatrixUtils.range(0, data.length-2), new int[] {1,2,3});
		
		ConditionalMutualInfoCalculatorMultiVariateGaussian condMiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateGaussian();
		condMiCalc.initialise(1, 1, 3);
		condMiCalc.setObservations(source, dest, others);
		
		// Now check that the Cholesky decomposition matches that for the
		//  expected covariance matrix:
		// (Note that this was computed from all available observations; here
		//  we're cutting off some of the first and last observations where there is
		//  no matching pair in the source/dest, so results will differ slightly)
		double[][] expectedCov = new double[][]
			{{0.9647348336238838, 3.206553219847798E-5, -0.0013932612411635703, 0.04178350449818639, -0.01494202491454874},
				  {3.206553219847798E-5, 0.9647348336238838, -0.055547119949140286, -0.0020067804899770256, 0.02693742557840663},
				  {-0.0013932612411635703, -0.055547119949140286, 1.0800072991165575, -0.009974731537464664, -2.1485745647111378E-4},
				  {0.04178350449818639, -0.0020067804899770256, -0.009974731537464664, 0.48319024794457854, -0.011333013565018278},
				  {-0.01494202491454874, 0.02693742557840663, -2.1485745647111378E-4, -0.011333013565018278, 0.5018806693655076}};
		double[][] expectedCholesky = MatrixUtils.CholeskyDecomposition(expectedCov);
		
		for (int r = 0; r < expectedCholesky.length; r++) {
			for (int c = 0; c < expectedCholesky[r].length; c++) {
				// As above, results will differ slightly., so allow larger than
				//  usual margin for error (plus amplification then occurs in
				//  computing the Cholesky decomposition):
				assertEquals(expectedCholesky[r][c], condMiCalc.L[r][c], 0.001);
			}
		}
		
		//============================
		// Case 2: Covariance conditioned on other variables:
		// The joint covariance matrix here is just what it would be if we were
		//  measuring the complete transfer entropy 
		
		source = MatrixUtils.selectRowsAndColumns(data,
				MatrixUtils.range(0, data.length-2), new int[] {1});
		dest = MatrixUtils.selectRowsAndColumns(data,
				MatrixUtils.range(1, data.length-1), new int[] {0});
		others = MatrixUtils.selectRowsAndColumns(data,
				MatrixUtils.range(0, data.length-2), new int[] {0,2,3});
		
		condMiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateGaussian();
		condMiCalc.initialise(1, 1, 3);
		condMiCalc.setObservations(source, dest, others);
		
		// Now check that the Cholesky decomposition matches that for the
		//  expected covariance matrix:
		// (Note that this was computed from all available observations; here
		//  we're cutting off some of the first and last observations where there is
		//  no matching pair in the source/dest, so results will differ slightly)
		expectedCov = new double[][]
			{{1.0800072991165575, -0.055547119949140286, -0.0013932612411635703, -0.009974731537464664, -2.1485745647111378E-4},
				  {-0.055547119949140286, 0.9647348336238838, 3.206553219847798E-5, -0.0020067804899770256, 0.02693742557840663},
				  {-0.0013932612411635703, 3.206553219847798E-5, 0.9647348336238838, 0.04178350449818639, -0.01494202491454874},
				  {-0.009974731537464664, -0.0020067804899770256, 0.04178350449818639, 0.48319024794457854, -0.011333013565018278},
				  {-2.1485745647111378E-4, 0.02693742557840663, -0.01494202491454874, -0.011333013565018278, 0.5018806693655076}};
		expectedCholesky = MatrixUtils.CholeskyDecomposition(expectedCov);
		
		for (int r = 0; r < expectedCholesky.length; r++) {
			for (int c = 0; c < expectedCholesky[r].length; c++) {
				// As above, results will differ slightly., so allow larger than
				//  usual margin for error (plus amplification then occurs in
				//  computing the Cholesky decomposition):
				assertEquals(expectedCholesky[r][c], condMiCalc.L[r][c], 0.001);
			}
		}
		
		// For future reference, the covariance matrix for TE with k=2 on this
		//  example (same source and dest) should be (verified with octave):
		/* expectedCov = new double[][]
				{{1.0800072991165575, -0.055547119949140286, -0.0013932612411635703, -0.020520351423877373, -0.009974731537464664, -2.1485745647111378E-4},
				  {-0.055547119949140286, 0.9647348336238838, 3.206553219847798E-5, 0.05415562847372847, -0.0020067804899770256, 0.02693742557840663},
				  {-0.0013932612411635703, 3.206553219847798E-5, 0.9647348336238838, 3.206553219847798E-5, 0.04178350449818639, -0.01494202491454874},
				  {-0.020520351423877373, 0.05415562847372847, 3.206553219847798E-5, 0.9647348336238838, 0.350905073828977, -0.013825394184539444},
				  {-0.009974731537464664, -0.0020067804899770256, 0.04178350449818639, 0.350905073828977, 0.48319024794457854, -0.011333013565018278},
				  {-2.1485745647111378E-4, 0.02693742557840663, -0.01494202491454874, -0.013825394184539444, -0.011333013565018278, 0.5018806693655076}};
		*/

	}
	
	/**
	 * Test whether, if the conditional variable has zero covariance,
	 *  that the method just returns the MI
	 *  
	 * @throws Exception
	 */
	public void testZeroCovarianceConditional() throws Exception {
		ConditionalMutualInfoCalculatorMultiVariateGaussian condMiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateGaussian();
		condMiCalc.initialise(1, 1, 1);
		double covar1 = 1;
		double covar2 = 0.8;
		double crossCovar = 0.6;
		double[][] covariance = new double[][] {
				{covar1, crossCovar, 0.0},
				{crossCovar, covar2, 0.0},
				{0.0, 0.0, 0.0}};
		condMiCalc.setCovariance(covariance, false);
		double condMi = condMiCalc.computeAverageLocalOfObservations();
		assertEquals(0.5 * Math.log(covar1 * covar2 / (covar1 * covar2 - crossCovar*crossCovar)),
				condMi, 0.0000000001);
	}
	
	public void testBiasCorrectionDoesNotChangeAnalyticPValue() throws Exception {
		ConditionalMutualInfoCalculatorMultiVariateGaussian cmiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateGaussian();
		
		int dimensions = 2;
		int timeSteps = 100;
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] condData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		cmiCalc.setProperty(MutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION, "false");
		cmiCalc.initialise(dimensions, dimensions, dimensions);
		cmiCalc.setObservations(sourceData, destData, condData);
		double avNotBiasCorrected = cmiCalc.computeAverageLocalOfObservations();
		ChiSquareMeasurementDistribution distroNotBiasCorrected = cmiCalc.computeSignificance();
		assertEquals(avNotBiasCorrected, distroNotBiasCorrected.actualValue, 0.0000001);
		
		// Now run again with bias correction:
		cmiCalc.setProperty(MutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION, "true");
		cmiCalc.initialise(dimensions, dimensions, dimensions);
		cmiCalc.setObservations(sourceData, destData, condData);
		double avBiasCorrected = cmiCalc.computeAverageLocalOfObservations();
		ChiSquareMeasurementDistribution distroBiasCorrected = cmiCalc.computeSignificance();
		assertEquals(avBiasCorrected, distroBiasCorrected.actualValue, 0.0000001);
		// And now check that the pValues are unchanged whether we bias correct or not:
		assertEquals(distroNotBiasCorrected.pValue, distroBiasCorrected.pValue);
	}

	protected int timeStepsDepCheck = 100;

	/**
	 * Check that the tests of linear dependencies that we also use for the MI calculator
	 * apply for conditional MI when the conditional is empty in some fashion
	 * (i.e. no dimensions, null or all constants)
	 * 
	 */
	public void testHandlingLinearDependenciesNoConditional() throws Exception {
		// 2D array with 0 variables (columns):
		checkHandlingLinearDependenciesUnrelatedOrNoConditional(new double[timeStepsDepCheck][0], 0);
		
		// Array of constants for each variable
		for (int c = 1; c < 5; c++) {
			double[][] conditional = new double[timeStepsDepCheck][c];
			for (int j = 0; j < c; j++) {
				MatrixUtils.copyIntoColumn(conditional, j,
						MatrixUtils.constantArray(timeStepsDepCheck, j)); // Just use a constant, not necessraily zeros
			}
			checkHandlingLinearDependenciesUnrelatedOrNoConditional(conditional, c);
		}
		// Null conditional
		checkHandlingLinearDependenciesUnrelatedOrNoConditional(null, 0);
	}
	
	/**
	 * Check that the tests of linear dependencies that we also use for the MI calculator
	 * apply for conditional MI when the conditional is irrelevant to source and dest
	 * 
	 */
	public void testHandlingLinearDependenciesIrrelevantConditional() throws Exception {
		RandomGenerator rg = new RandomGenerator();
		for (int c = 0; c < 5; c++) {
			// Try with just noisy conditionals:
			double[][] condData = rg.generateNormalData(timeStepsDepCheck, c,
				0, 1);
			checkHandlingLinearDependenciesUnrelatedOrNoConditional(condData, c);
			// Try also if one of the variables is redundant with others
			if (c > 2) {
				condData = rg.generateNormalData(timeStepsDepCheck, c, 0, 1);
				MatrixUtils.copyIntoColumn(condData, 2,
						MatrixUtils.add(
								MatrixUtils.selectColumn(condData, 0),
								MatrixUtils.selectColumn(condData, 1)));
				checkHandlingLinearDependenciesUnrelatedOrNoConditional(condData, c);
			}
		}
	}

	/**
	 * Check that the tests of linear dependencies that we also use for the MI calculator
	 * apply for conditional MI when the conditional is empty in some fashion or irrelevant
	 * 
	 * @param emptyConditional
	 * @param conditionalDims
	 * @throws Exception
	 */
	protected void checkHandlingLinearDependenciesUnrelatedOrNoConditional(
			double[][] emptyConditional, int conditionalDims) throws Exception {
		
		ConditionalMutualInfoCalculatorMultiVariateGaussian condMiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateGaussian();
				
		int dimensions = 2; // Assumed to be >-= 2
		assertTrue(dimensions == 2);
		RandomGenerator rg = new RandomGenerator();

		// Generate some random data and do an MI on copied version
		//  - both dimensions are copied
		double[][] sourceData = rg.generateNormalData(timeStepsDepCheck, dimensions,
				0, 1);
		double[][] destData = MatrixUtils.arrayCopy(sourceData);
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, emptyConditional);
		double miCopied = condMiCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isInfinite(miCopied));
		
		// Now overwrite one of the columns, and check result is still infinite
		//  with two columns across the variables the same (only one dimension copied):
		MatrixUtils.copyIntoColumn(sourceData, 1,
				rg.generateNormalData(timeStepsDepCheck, 0, 1));
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, emptyConditional);
		double miOneCopied = condMiCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isInfinite(miOneCopied));
		
		// Now make the source columns a copy of themselves, and check that it can ignore this
		double[] indpColumn = MatrixUtils.selectColumn(sourceData, 1);
		MatrixUtils.copyIntoColumn(sourceData, 0,
				indpColumn);
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, emptyConditional);
		double mi1SourceFrom2 = condMiCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isFinite(mi1SourceFrom2));
		condMiCalc.initialise(1,dimensions, conditionalDims);
		// Should be the same whether we compute this from both source columns or only 1
		double[][] source1Column = new double[timeStepsDepCheck][1];
		MatrixUtils.copyIntoColumn(source1Column, 0, indpColumn);
		condMiCalc.setObservations(source1Column, destData, emptyConditional);
		double mi1DSourceFrom1 = condMiCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isFinite(mi1DSourceFrom1));
		assertEquals(mi1DSourceFrom1, mi1SourceFrom2, 0.00000001);
		// And check it works if we flip source and dest:
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(destData, sourceData, emptyConditional);
		double mi1SourceFrom2Flipped = condMiCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isFinite(mi1SourceFrom2Flipped));
		assertEquals(mi1DSourceFrom1, mi1SourceFrom2Flipped, 0.00000001);
		
		// Refresh data, with a third dependent source variable in and check that this doesn't change things: 
		sourceData = rg.generateNormalData(timeStepsDepCheck, dimensions,
				0, 1);
		destData = rg.generateNormalData(timeStepsDepCheck, dimensions + 1,
				0, 1);
		MatrixUtils.copyIntoColumn(destData, dimensions,
				MatrixUtils.add(MatrixUtils.selectColumn(destData, 0), MatrixUtils.selectColumn(destData, 1)));
		condMiCalc.initialise(dimensions, dimensions + 1, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, emptyConditional);
		double mi1RedundantDest = condMiCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isFinite(mi1RedundantDest));
		// check that it works flipped
		condMiCalc.initialise(dimensions+1, dimensions, conditionalDims);
		condMiCalc.setObservations(destData, sourceData, emptyConditional);
		double mi1Redundantsource = condMiCalc.computeAverageLocalOfObservations();
		assertEquals(mi1RedundantDest, mi1Redundantsource, 1e-7);
		// Check that it's the same if we only use the independent variables:
		destData = MatrixUtils.selectColumns(destData, new int[] {0, 1});
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, emptyConditional);
		double miDestWithoutRedundant = condMiCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isFinite(miDestWithoutRedundant));
		assertEquals(mi1RedundantDest, miDestWithoutRedundant, 0.00000001);
		
		// First check that the MI is not precisely zero if only one of the sub-variables is a zero:
		MatrixUtils.copyIntoColumn(destData, 0,
				MatrixUtils.constantArray(timeStepsDepCheck, 0));
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, emptyConditional);
		double miOneColZero = condMiCalc.computeAverageLocalOfObservations();
		// Check that the MI here is not precisely zero, it shouldn't be for a finite sample length
		// System.out.println(miOneColZero);
		assertTrue(Math.abs(miOneColZero) > 1e-12);
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(destData, sourceData, emptyConditional); // check if we swap the variables all is the same
		double miOneColZeroSwapped = condMiCalc.computeAverageLocalOfObservations();
		assertEquals(miOneColZero, miOneColZeroSwapped, 1e-7);
		// Now do the same to a sub-variable of the source:
		MatrixUtils.copyIntoColumn(sourceData, 1,
				MatrixUtils.constantArray(timeStepsDepCheck, 0));
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, emptyConditional);
		double miOneColZeroSourceAlso = condMiCalc.computeAverageLocalOfObservations();
		assertTrue(Math.abs(miOneColZeroSourceAlso) > 1e-12);
		// Check that the result is the same if we only supplied the non-zero columns:
		condMiCalc.initialise(1, 1, conditionalDims);
		condMiCalc.setObservations(MatrixUtils.selectColumns(sourceData, new int[] {0}),
				MatrixUtils.selectColumns(destData, new int[] {1}),
				emptyConditional);
		double miOneColZeroSourceAlsoUnivariateCalcs = condMiCalc.computeAverageLocalOfObservations();
		assertEquals(miOneColZeroSourceAlsoUnivariateCalcs, miOneColZeroSourceAlso, 1e-7);
		
		// Test that we get a zero if no variance within any source or target variables
		MatrixUtils.copyIntoColumn(destData, 1,
				MatrixUtils.constantArray(timeStepsDepCheck, 0));
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, emptyConditional);
		double miDestAllZeros = condMiCalc.computeAverageLocalOfObservations();
		assertEquals(miDestAllZeros, 0, 1e-13);
		// And the other way around:
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(destData, sourceData, emptyConditional);
		double miSourceAllZeros = condMiCalc.computeAverageLocalOfObservations();
		assertEquals(miSourceAllZeros, 0, 1e-13);
		
		// And if only a univariate is all zeros:
		condMiCalc.initialise(1, dimensions, conditionalDims);
		condMiCalc.setObservations(MatrixUtils.selectColumns(destData, 1, 1),
				rg.generateNormalData(timeStepsDepCheck, dimensions, 0, 1), emptyConditional);
		double miSourceZero = condMiCalc.computeAverageLocalOfObservations();
		assertEquals(miSourceZero, 0, 1e-13);
		// and the other way around:
		condMiCalc.initialise(dimensions, 1, conditionalDims);
		condMiCalc.setObservations(rg.generateNormalData(timeStepsDepCheck, dimensions, 0, 1),
				MatrixUtils.selectColumns(destData, 1, 1), emptyConditional);
		double miDestZero = condMiCalc.computeAverageLocalOfObservations();
		assertEquals(miDestZero, 0, 1e-13);
	}

	public void testHandlingLinearDependenciesWithConditionalRelated() throws Exception {
		
		ConditionalMutualInfoCalculatorMultiVariateGaussian condMiCalc =
				new ConditionalMutualInfoCalculatorMultiVariateGaussian();
				
		int dimensions = 2; // Assumed to be 2
		assertTrue(dimensions == 2);
		int conditionalDims = 2; // Assumed to be 2
		assertTrue(conditionalDims == 2);
		RandomGenerator rg = new RandomGenerator();

		// Generate some random data and do an MI on copied version
		//  - both dimensions are copied
		double[][] sourceData = rg.generateNormalData(timeStepsDepCheck, dimensions,
				0, 1);
		double[][] condData = rg.generateNormalData(timeStepsDepCheck, dimensions,
				0, 1);
		double[][] destData = MatrixUtils.arrayCopy(condData);

		// Conditional renders dest fully redundant (straight copy of conditional) -> 0
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, condData);
		double miDestRedundantWithCond = condMiCalc.computeAverageLocalOfObservations();
		assertEquals(0, miDestRedundantWithCond, 1e-7);
		// Or is sum of two variables
		MatrixUtils.copyIntoColumn(destData, 0,
				MatrixUtils.add(MatrixUtils.selectColumn(condData, 0),
							MatrixUtils.selectColumn(condData, 1)));
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, condData);
		double miDestRedundantWithCondSum = condMiCalc.computeAverageLocalOfObservations();
		assertEquals(0, miDestRedundantWithCondSum, 1e-7);
		// And the other way around:
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(destData, sourceData, condData);
		assertEquals(0, condMiCalc.computeAverageLocalOfObservations(), 1e-7);
		// And what if we only use a univariate dest:
		condMiCalc.initialise(dimensions, 1, conditionalDims);
		condMiCalc.setObservations(sourceData, MatrixUtils.selectColumns(destData, 0, 1), condData);
		assertEquals(0, condMiCalc.computeAverageLocalOfObservations(), 1e-7);		
		
		// Conditional rendering a dest column (0) fully dependent doesn't change other result (for its column 1)
		MatrixUtils.copyIntoColumn(destData, 1, rg.generateNormalData(timeStepsDepCheck, 0, 1));
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, condData);
		double miFromOneIndpCol = condMiCalc.computeAverageLocalOfObservations();
		condMiCalc.initialise(dimensions, 1, conditionalDims);
		condMiCalc.setObservations(sourceData, MatrixUtils.selectColumns(destData, 1, 1), condData);
		assertEquals(miFromOneIndpCol, condMiCalc.computeAverageLocalOfObservations(), 1e-7);
		
		// Conditional renders source and dest dependent (e.g. by summing the two) - inf
		// Whether it's only one variable or both
		destData = MatrixUtils.add(sourceData, condData);
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(sourceData, destData, condData);
		assertTrue(Double.isInfinite(condMiCalc.computeAverageLocalOfObservations()));
		condMiCalc.initialise(dimensions, 1, conditionalDims);
		condMiCalc.setObservations(sourceData, MatrixUtils.selectColumns(destData, 0, 1), condData);
		assertTrue(Double.isInfinite(condMiCalc.computeAverageLocalOfObservations()));
		// And now flip the source and dest
		condMiCalc.initialise(dimensions, dimensions, conditionalDims);
		condMiCalc.setObservations(destData, sourceData, condData);
		assertTrue(Double.isInfinite(condMiCalc.computeAverageLocalOfObservations()));
		condMiCalc.initialise(1, dimensions, conditionalDims);
		condMiCalc.setObservations(MatrixUtils.selectColumns(destData, 0, 1), sourceData, condData);
		assertTrue(Double.isInfinite(condMiCalc.computeAverageLocalOfObservations()));
	}
}
