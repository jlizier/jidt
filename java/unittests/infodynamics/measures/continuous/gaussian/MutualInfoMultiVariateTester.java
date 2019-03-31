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

import junit.framework.TestCase;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.MatrixUtilsTest;
import infodynamics.utils.RandomGenerator;

public class MutualInfoMultiVariateTester extends TestCase {

	public void testCovarianceDoesntMatchDimensions() throws Exception {
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		miCalc.initialise(1, 1);
		boolean caughtException = false;
		// Check that we catch a covariance matrix which doesn't have
		//  the right number of rows
		try {
			miCalc.setCovariance(new double[][] {{2,1,0.5}, {1,2,0.5}, {0.5,0.5,2}}, 1);
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		// Check that we catch a covariance matrix which isn't square
		caughtException = false;
		try {
			miCalc.setCovariance(new double[][] {{2,1}, {1,2,0.5}}, 1);
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		// Check that we catch a covariance matrix which isn't symmetric
		caughtException = false;
		try {
			miCalc.setCovariance(new double[][] {{2,1}, {1.000001,2}}, 1);
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		// Check that no exception is thrown for an ok covariance matrix
		caughtException = false;
		double[][] goodCovariance = new double[][] {{2,1}, {1,2}};
		try {
			miCalc.setCovariance(goodCovariance, 93);
		} catch (Exception e) {
			caughtException = true;
		}
		assertFalse(caughtException);
		// and that this covariance has been set (by verifying the Cholesky
		// decomposition of it is stored):
		MatrixUtilsTest.checkMatrix(MatrixUtils.CholeskyDecomposition(goodCovariance),
				miCalc.L, 0.00001);
	}
	
	public void testAnalyticMatchesCouplingValue() throws Exception {
		double[][] covarianceMatrix = {{1.0, 0}, {0, 1.0}};
		
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		
		for (double covar = 0.0; covar < 1.0; covar += 0.01) {
			// Insert the new covariance into the matrix:
			covarianceMatrix[0][1] = covar;
			covarianceMatrix[1][0] = covar;

			miCalc.initialise(1, 1);
			miCalc.setCovariance(covarianceMatrix, 1);
			assertEquals(-0.5 * Math.log(1.0 - covar*covar),
					miCalc.computeAverageLocalOfObservations(), 0.00000000001);
		}
	}
	
	public void testMIfromSuppliedCovariance() throws Exception {
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();

		double[][] covarianceMatrix = {{5, 3}, {3, 4}}; // det is 11
		miCalc.initialise(1, 1);
		miCalc.setCovariance(covarianceMatrix, 100);
		assertEquals(0.5 * Math.log(20.0 / 11.0),
				miCalc.computeAverageLocalOfObservations(), 0.00000000001);
		
		double[][] covarianceMatrix2 = {{5, 3, 1}, {3, 4, 1.5}, {1, 1.5, 2}}; // det is 15.75
		miCalc.initialise(2, 1);
		miCalc.setCovariance(covarianceMatrix2, 100);
		assertEquals(0.5 * Math.log(11.0 * 2.0 / 15.75), // marginal dets are 11 and 2
				miCalc.computeAverageLocalOfObservations(), 0.00000000001);

	}
	
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		
		int dimensions = 2;
		int timeSteps = 100;
		
		miCalc.initialise(dimensions, dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		miCalc.setObservations(sourceData, destData);
		
		//miCalc.setDebug(true);
		double mi = miCalc.computeAverageLocalOfObservations();
		//miCalc.setDebug(false);
		//double[] miLocal = miCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", mi);
		
		// Now look at statistical significance tests
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(
				timeSteps, 100);

		EmpiricalMeasurementDistribution measDist =
				miCalc.computeSignificance(newOrderings);
		
		System.out.printf("pValue of sig test was %.3f\n", measDist.pValue);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double averageCheck1 = miCalc.computeAverageLocalOfObservations();
			assertEquals(mi, averageCheck1);
		}
	}
	
	public void testBiasCorrectionDoesNotChangeAnalyticPValue() throws Exception {
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		
		int dimensions = 2;
		int timeSteps = 100;
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		miCalc.setProperty(MutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION, "false");
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(sourceData, destData);
		double avNotBiasCorrected = miCalc.computeAverageLocalOfObservations();
		ChiSquareMeasurementDistribution distroNotBiasCorrected = miCalc.computeSignificance();
		assertEquals(avNotBiasCorrected, distroNotBiasCorrected.actualValue, 0.0000001);
		
		// Now run again with bias correction:
		miCalc.setProperty(MutualInfoCalculatorMultiVariateGaussian.PROP_BIAS_CORRECTION, "true");
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(sourceData, destData);
		double avBiasCorrected = miCalc.computeAverageLocalOfObservations();
		ChiSquareMeasurementDistribution distroBiasCorrected = miCalc.computeSignificance();
		assertEquals(avBiasCorrected, distroBiasCorrected.actualValue, 0.0000001);
		// And now check that the pValues are unchanged whether we bias correct or not:
		assertEquals(distroNotBiasCorrected.pValue, distroBiasCorrected.pValue);
	}

	public void testHandlingLinearDependencies() throws Exception {
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		
		int dimensions = 2; // Assumed to be >-= 2
		assertTrue(dimensions == 2);
		int timeSteps = 100;
		RandomGenerator rg = new RandomGenerator();
		
		// Generate some random data and do an MI on copied version
		//  - both dimensions are copied
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = MatrixUtils.arrayCopy(sourceData);
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(sourceData, destData);
		double miCopied = miCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isInfinite(miCopied));
		
		// Now overwrite one of the columns, and check result is still infinite
		//  with two columns across the variables the same (only one dimension copied):
		MatrixUtils.copyIntoColumn(sourceData, 1,
				rg.generateNormalData(timeSteps, 0, 1));
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(sourceData, destData);
		double miOneCopied = miCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isInfinite(miOneCopied));
		
		// Now make the source columns a copy of themselves, and check that it can ignore this
		double[] indpColumn = MatrixUtils.selectColumn(sourceData, 1);
		MatrixUtils.copyIntoColumn(sourceData, 0,
				indpColumn);
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(sourceData, destData);
		double mi1SourceFrom2 = miCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isFinite(mi1SourceFrom2));
		miCalc.initialise(1,dimensions);
		// Should be the same whether we compute this from both source columns or only 1
		double[][] source1Column = new double[timeSteps][1];
		MatrixUtils.copyIntoColumn(source1Column, 0, indpColumn);
		miCalc.setObservations(source1Column, destData);
		double mi1DSourceFrom1 = miCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isFinite(mi1DSourceFrom1));
		assertEquals(mi1DSourceFrom1, mi1SourceFrom2, 0.00000001);
		// And check it works if we flip source and dest:
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(destData, sourceData);
		double mi1SourceFrom2Flipped = miCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isFinite(mi1SourceFrom2Flipped));
		assertEquals(mi1DSourceFrom1, mi1SourceFrom2Flipped, 0.00000001);
		
		// Refresh data, with a third dependent source variable in and check that this doesn't change things: 
		sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		destData = rg.generateNormalData(timeSteps, dimensions + 1,
				0, 1);
		MatrixUtils.copyIntoColumn(destData, dimensions,
				MatrixUtils.add(MatrixUtils.selectColumn(destData, 0), MatrixUtils.selectColumn(destData, 1)));
		miCalc.initialise(dimensions, dimensions + 1);
		miCalc.setObservations(sourceData, destData);
		double mi1RedundantDest = miCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isFinite(mi1RedundantDest));
		// check that it works flipped
		miCalc.initialise(dimensions+1, dimensions);
		miCalc.setObservations(destData, sourceData);
		double mi1Redundantsource = miCalc.computeAverageLocalOfObservations();
		assertEquals(mi1RedundantDest, mi1Redundantsource, 1e-7);
		// Check that it's the same if we only use the independent variables:
		destData = MatrixUtils.selectColumns(destData, new int[] {0, 1});
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(sourceData, destData);
		double miDestWithoutRedundant = miCalc.computeAverageLocalOfObservations();
		assertTrue(Double.isFinite(miDestWithoutRedundant));
		assertEquals(mi1RedundantDest, miDestWithoutRedundant, 0.00000001);
		
		// First check that the MI is not precisely zero if only one of the sub-variables is a zero:
		MatrixUtils.copyIntoColumn(destData, 0,
				MatrixUtils.constantArray(timeSteps, 0));
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(sourceData, destData);
		double miOneColZero = miCalc.computeAverageLocalOfObservations();
		// Check that the MI here is not precisely zero, it shouldn't be for a finite sample length
		// System.out.println(miOneColZero);
		assertTrue(Math.abs(miOneColZero) > 1e-12);
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(destData, sourceData); // check if we swap the variables all is the same
		double miOneColZeroSwapped = miCalc.computeAverageLocalOfObservations();
		assertEquals(miOneColZero, miOneColZeroSwapped, 1e-7);
		// Now do the same to a sub-variable of the source:
		MatrixUtils.copyIntoColumn(sourceData, 1,
				MatrixUtils.constantArray(timeSteps, 0));
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(sourceData, destData);
		double miOneColZeroSourceAlso = miCalc.computeAverageLocalOfObservations();
		assertTrue(Math.abs(miOneColZeroSourceAlso) > 1e-12);
		// Check that the result is the same if we only supplied the non-zero columns:
		miCalc.initialise(1, 1);
		miCalc.setObservations(MatrixUtils.selectColumn(sourceData, 0),
				MatrixUtils.selectColumn(destData, 1));
		double miOneColZeroSourceAlsoUnivariateCalcs = miCalc.computeAverageLocalOfObservations();
		assertEquals(miOneColZeroSourceAlsoUnivariateCalcs, miOneColZeroSourceAlso, 1e-7);
		
		// Test that we get a zero if no variance within any source or target variables
		MatrixUtils.copyIntoColumn(destData, 1,
				MatrixUtils.constantArray(timeSteps, 0));
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(sourceData, destData);
		double miDestAllZeros = miCalc.computeAverageLocalOfObservations();
		assertEquals(miDestAllZeros, 0, 1e-13);
		// And the other way around:
		miCalc.initialise(dimensions, dimensions);
		miCalc.setObservations(destData, sourceData);
		double miSourceAllZeros = miCalc.computeAverageLocalOfObservations();
		assertEquals(miSourceAllZeros, 0, 1e-13);
		
		// And if only a univariate is all zeros:
		miCalc.initialise(1, dimensions);
		miCalc.setObservations(MatrixUtils.selectColumns(destData, 1, 1),
				rg.generateNormalData(timeSteps, dimensions, 0, 1));
		double miSourceZero = miCalc.computeAverageLocalOfObservations();
		assertEquals(miSourceZero, 0, 1e-13);
		// and the other way around:
		miCalc.initialise(dimensions, 1);
		miCalc.setObservations(rg.generateNormalData(timeSteps, dimensions, 0, 1),
				MatrixUtils.selectColumns(destData, 1, 1));
		double miDestZero = miCalc.computeAverageLocalOfObservations();
		assertEquals(miDestZero, 0, 1e-13);
	}
}
