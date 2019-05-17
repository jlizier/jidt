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

package infodynamics.utils;

import java.util.List;

import junit.framework.TestCase;

/**
 * Test functionality of MatrixUtils methods.
 * 
 * @author Joseph Lizier joseph.lizier_at_gmail.com
 *
 */
public class MatrixUtilsTest extends TestCase {

	private static double OCTAVE_RESOLUTION = 0.00001;

	public void testIdentityMatrix() {
		// Test identity matrix generation, including for size 0
		for (int n = 0; n < 10; n++) {
			double[][] I = MatrixUtils.identityMatrix(n);
			assertNotNull(I);
			assertEquals(n, I.length);
			for (int i = 0; i < n; i++) {
				assertEquals(n, I[i].length);
				for (int j = 0; j < n; j++) {
					assertEquals(i == j ? 1 : 0, I[i][j], 0.000000001);
				}
			}
		}
	}
	
	public void testCovariance() throws Exception {
		// Load a test data file which contains a time series and
		//  covariance matrix as computed by matlab:
		// We'll just take the first two columns from this data set
		OctaveFileReader ofr = new OctaveFileReader(
			"demos/data/Network-GaussianLinear-N100-T100-p0.04-b0.50-c0.50-dir-disc-repeat1.txt");
		double[][] data = ofr.getDouble2DMatrix("timeseries");
		double[][] expectedCovariance = ofr.getDouble2DMatrix("empiricalCovariance");
		double[][] computedCovariance = MatrixUtils.covarianceMatrix(data);
		checkMatrix(expectedCovariance, computedCovariance, 0.00001);
		
		// And test that it's still correct if we supply the data in 2 separate
		//  parts:
		double[][] part1 = MatrixUtils.selectColumns(data,
				MatrixUtils.range(0, 49));
		double[][] part2 = MatrixUtils.selectColumns(data,
				MatrixUtils.range(50, 99));
		double[][] split2ComputedCovariance = MatrixUtils.covarianceMatrix(
				part1, part2);
		checkMatrix(expectedCovariance, split2ComputedCovariance, 0.00001);

		// And test that it's still correct if we supply the data in 3 separate
		//  parts:
		double[][] part2a = MatrixUtils.selectColumns(data,
				MatrixUtils.range(50, 74));
		double[][] part2b = MatrixUtils.selectColumns(data,
				MatrixUtils.range(75, 99));
		double[][] split3ComputedCovariance = MatrixUtils.covarianceMatrix(
				part1, part2a, part2b);
		checkMatrix(expectedCovariance, split3ComputedCovariance, 0.00001);
	}
	
	/**
	 * Test our Cholesky decomposition implementation
	 * 
	 * @throws Exception
	 */
	public void testCholesky() throws Exception {
		
		// Check some ordinary Cholesky decompositions:
		
		double[][] A = {{6, 2, 3}, {2, 5, 1}, {3, 1, 4}};
		// Expected result from Octave:
		double[][] expectedL = {{2.44949, 0, 0}, {0.81650, 2.08167, 0},
				{1.22474, 0, 1.58114}};
		double[][] L = MatrixUtils.CholeskyDecomposition(A);
		checkMatrix(expectedL, L, OCTAVE_RESOLUTION);

		double[][] A2 = {{6, 2, 3, 1}, {2, 5, 1, 0.5}, {3, 1, 4, 2}, {1, 0.5, 2, 3}};
		// Expected result from Octave:
		double[][] expectedL2 = {{2.44949, 0, 0, 0}, {0.81650, 2.08167, 0, 0},
				{1.22474, 0, 1.58114, 0}, {0.40825, 0.08006, 0.94868, 1.38814}};
		double[][] L2 = MatrixUtils.CholeskyDecomposition(A2);
		checkMatrix(expectedL2, L2, OCTAVE_RESOLUTION);

		// Now check that it picks up asymmetric A:
		double[][] asymmetricA = {{6, 2, 3}, {2, 5, 1}, {3, 1.0001, 4}};
		boolean flaggedException = false;
		try {
			MatrixUtils.CholeskyDecomposition(asymmetricA);
		} catch (Exception e) {
			flaggedException = true;
		}
		assertTrue(flaggedException);
	}
	
	public void testCholeskyNotPositiveDefinite() throws Exception {

		// Now check that it picks up if A is not positive definite:
		//  Variable 2 is a scaled version of variable 1 here
		double[][] notpositiveDefiniteA = {{1, 2}, {2, 4}};
		boolean flaggedException = false;
		int problematicRow = -1;
		try {
			MatrixUtils.CholeskyDecomposition(notpositiveDefiniteA);
		} catch (NonPositiveDefiniteMatrixException npdme) {
			flaggedException = true;
			problematicRow = npdme.problematicRow;
		}
		assertTrue(flaggedException);
		assertEquals(1, problematicRow);
				
		// Test that we pick up a variable which duplicates another
		double[][] copiedVariableA = {{1, 1.9, 1}, {1.9, 4, 1.9}, {1, 1.9, 1}};
		problematicRow = -1;
		flaggedException = false;
		try {
			MatrixUtils.CholeskyDecomposition(copiedVariableA);
		} catch (NonPositiveDefiniteMatrixException npdme) {
			flaggedException = true;
			problematicRow = npdme.problematicRow;
		}
		assertTrue(flaggedException);
		assertEquals(2, problematicRow);
		double[][] copiedVariableRow1A = {{1, 1, 1.9}, {1, 1, 1.9}, {1.9, 1.9, 4}};
		problematicRow = -1;
		flaggedException = false;
		try {
			MatrixUtils.CholeskyDecomposition(copiedVariableRow1A);
		} catch (NonPositiveDefiniteMatrixException npdme) {
			flaggedException = true;
			problematicRow = npdme.problematicRow;
		}
		assertTrue(flaggedException);
		assertEquals(1, problematicRow);
		
		// Now do this with some real data:
		RandomGenerator rg = new RandomGenerator();
		int N = 1000;
		double[][] data = new double[N][3];
		MatrixUtils.copyIntoColumn(data, 0, rg.generateNormalData(N, 0, 1));
		MatrixUtils.copyIntoColumn(data, 1, rg.generateNormalData(N, 0, 1));
		MatrixUtils.copyIntoColumn(data, 2, rg.generateNormalData(N, 0, 1));
		double[][] covarianceMatrix = MatrixUtils.covarianceMatrix(data);
		// This one should be fine and not throw an exception:
		MatrixUtils.CholeskyDecomposition(covarianceMatrix);
		// But if we substitute in a linearly redundant variable:
		MatrixUtils.copyIntoColumn(data, 2,
				MatrixUtils.add(MatrixUtils.selectColumn(data, 0),
								MatrixUtils.selectColumn(data, 1)));
		covarianceMatrix = MatrixUtils.covarianceMatrix(data);
		problematicRow = -1;
		flaggedException = false;
		double[][] L;
		try {
			// This should throw an exception:
			L = MatrixUtils.CholeskyDecomposition(covarianceMatrix);
			// Debug prints in case it does not (because tolerance in the CholeskyDecomposition was too large):
			System.out.println(L[2][2]);
			System.out.println("Det " + MatrixUtils.determinantViaCholeskyResult(L) + " > 0, unexpectedly, but it happens - precision is approx: " + 
					java.lang.Math.ulp(0));
			MatrixUtils.printMatrix(System.out, covarianceMatrix);
			System.out.printf("%.6f, %.6f\n", covarianceMatrix[0][2] - covarianceMatrix[0][1] - covarianceMatrix[0][0],
					covarianceMatrix[1][2] - covarianceMatrix[1][1] - covarianceMatrix[1][0]);
		} catch (NonPositiveDefiniteMatrixException npdme) {
			flaggedException = true;
			problematicRow = npdme.problematicRow;
			L = null;
		}
		assertTrue(flaggedException);
		assertEquals(2, problematicRow);
		// Now if we swap the second and third variables in the data, it should still complain about the last one
		//  since that's the one it saw last:
		double[] temp = MatrixUtils.selectColumn(data, 2);
		MatrixUtils.copyIntoColumn(data, 2, MatrixUtils.selectColumn(data, 1));
		MatrixUtils.copyIntoColumn(data, 1, temp);
		covarianceMatrix = MatrixUtils.covarianceMatrix(data);
		problematicRow = -1;
		flaggedException = false;
		try {
			// This should throw an exception:
			L = MatrixUtils.CholeskyDecomposition(covarianceMatrix);
		} catch (NonPositiveDefiniteMatrixException npdme) {
			flaggedException = true;
			problematicRow = npdme.problematicRow;
			L = null;
		}
		assertTrue(flaggedException);
		assertEquals(2, problematicRow);
	}
	
	/**
	 * Test the inversion of symmetric positive definite matrices
	 * 
	 * @throws Exception
	 */
	public void testInverseOfSymmPosDefMatrices() throws Exception {
		// Check some ordinary matrices:
		
		double[][] A = {{6, 2, 3}, {2, 5, 1}, {3, 1, 4}};
		// Expected result from Octave:
		double[][] expectedInv = {{0.29231, -0.07692, -0.2},
				  {-0.07692, 0.23077, 0}, {-0.2, 0, 0.4}};
		double[][] inv = MatrixUtils.invertSymmPosDefMatrix(A);
		checkMatrix(expectedInv, inv, OCTAVE_RESOLUTION);

		double[][] A2 = {{6, 2, 3, 1}, {2, 5, 1, 0.5}, {3, 1, 4, 2}, {1, 0.5, 2, 3}};
		// Expected result from Octave:
		double[][] expectedInv2 = {{0.303393, -0.079840, -0.245509, 0.075848},
				  {-0.079840, 0.231537, 0.011976, -0.019960},
				  {-0.245509, 0.011976, 0.586826, -0.311377},
				  {0.075848, -0.019960, -0.311377, 0.518962}};
		double[][] inv2 = MatrixUtils.invertSymmPosDefMatrix(A2);
		checkMatrix(expectedInv2, inv2, OCTAVE_RESOLUTION);

		// Now check that it picks up asymmetric A:
		double[][] asymmetricA = {{6, 2, 3}, {2, 5, 1}, {3, 1.0001, 4}};
		boolean flaggedException = false;
		try {
			MatrixUtils.invertSymmPosDefMatrix(asymmetricA);
		} catch (Exception e) {
			flaggedException = true;
		}
		assertTrue(flaggedException);

		// Now check that it picks up if A is not positive definite:
		double[][] notpositiveDefiniteA = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
		flaggedException = false;
		try {
			MatrixUtils.invertSymmPosDefMatrix(notpositiveDefiniteA);
		} catch (Exception e) {
			flaggedException = true;
		}
		assertTrue(flaggedException);
	}
	
	/**
	 * Test the solving of matrix equations via Cholesky decomposition.
	 * Solving A*X = B
	 * 
	 * @throws Exception
	 */
	public void testSolvingMatrixEquationsOfSymmPosDefMatrices() throws Exception {
		// Check some ordinary matrices:
		
		double[][] A = {{6, 2, 3}, {2, 5, 1}, {3, 1, 4}};
		double[][] B = {{5}, {4}, {3}};
		// Expected result from Octave:
		double[][] expectedX = {{0.55385}, {0.53846}, {0.20000}};
		double[][] X = MatrixUtils.solveViaCholeskyResult(
				MatrixUtils.CholeskyDecomposition(A), B);
		checkMatrix(expectedX, X, OCTAVE_RESOLUTION);

		//  Check more complicated example
		double[][] A2 = {{6, 2, 3, 1}, {2, 5, 1, 0.5}, {3, 1, 4, 2}, {1, 0.5, 2, 3}};
		double[][] B2 = {{10, 5, 4, 12}, {4, 6, -1, 4.3}, {20, 1, 0, -5}, {6, 3, 2, 1}};
		double[][] expectedX2 = {{-1.740519, 1.019960, 1.445110, 4.600798},
			   {0.247505, 0.942116, -0.590818, -0.042315},
			   {7.461078, -1.502994, -1.616766, -6.140120},
			   {-2.435130, 1.504990, 1.361277, 2.900200}};
		double[][] X2 = MatrixUtils.solveViaCholeskyResult(
				MatrixUtils.CholeskyDecomposition(A2), B2);
		checkMatrix(expectedX2, X2, OCTAVE_RESOLUTION);
		
		// TODO Check error conditions
	}
	
	public void testDeterminant() throws Exception {

		// test some error conditions:
		double[][] AnonSquare = {{6, 2}, {2, 5, 1}, {3, 1, 4}};
		boolean flaggedException = false;
		try {
			MatrixUtils.determinant(AnonSquare);
		} catch (Exception e) {
			flaggedException = true;
		}
		assertTrue(flaggedException);
		
		// Test some simple examples:
		double[][] A1 = {{3.445454}};
		assertEquals(3.445454, MatrixUtils.determinant(A1), OCTAVE_RESOLUTION);
		
		double[][] A2 = {{6, 2}, {2, 5}};
		assertEquals(26, MatrixUtils.determinant(A2), OCTAVE_RESOLUTION);

		// Check against value computed by Octave:
		double[][] A = {{6, 2, 3}, {2, 5, 1}, {3, 1, 4}};
		assertEquals(65, MatrixUtils.determinant(A), OCTAVE_RESOLUTION);
		
		// Check zero determinant case
		double[][] AzeroDet = {{6, 2, 3}, {2, 5, 1}, {10, -14, 5}};
		assertEquals(0, MatrixUtils.determinant(AzeroDet), OCTAVE_RESOLUTION);		
	}
	
	public void testDeterminantSymmPosDef() throws Exception {
		// Test some simple examples:
		double[][] A1 = {{3.445454}};
		assertEquals(3.445454,
				MatrixUtils.determinantSymmPosDefMatrix(A1), OCTAVE_RESOLUTION);
		
		double[][] A2 = {{6, 2}, {2, 5}};
		assertEquals(26, MatrixUtils.determinantSymmPosDefMatrix(A2), OCTAVE_RESOLUTION);

		// Check against value computed by Octave:
		double[][] A = {{6, 2, 3}, {2, 5, 1}, {3, 1, 4}};
		assertEquals(65, MatrixUtils.determinantSymmPosDefMatrix(A), OCTAVE_RESOLUTION);
		
		// Now check that it picks up asymmetric A:
		double[][] asymmetricA = {{6, 2, 3}, {2, 5, 1}, {3, 1.0001, 4}};
		boolean flaggedException = false;
		try {
			MatrixUtils.determinantSymmPosDefMatrix(asymmetricA);
		} catch (Exception e) {
			flaggedException = true;
		}
		assertTrue(flaggedException);

		// Now check that it picks up if A is not positive definite:
		double[][] notpositiveDefiniteA = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
		flaggedException = false;
		try {
			MatrixUtils.determinantSymmPosDefMatrix(notpositiveDefiniteA);
		} catch (Exception e) {
			flaggedException = true;
		}
		assertTrue(flaggedException);
	}
	
	public void testSortIndices() {
		double[] array1 = {0.1, 0.2, 0.3, 0.4, 0.5};
		checkArray(new int[] {0, 1, 2, 3, 4}, MatrixUtils.sortIndices(array1));

		double[] array2 = {0.5, 0.4, 0.3, 0.2, 0.1};
		checkArray(new int[] {4, 3, 2, 1, 0}, MatrixUtils.sortIndices(array2));

		double[] array3 = {0.3, 0.1, 0.5, 0.4, 0.2};
		checkArray(new int[] {1, 4, 0, 3, 2}, MatrixUtils.sortIndices(array3));
	}
	
	public void testDelayEmbeddings() throws Exception {
		double[] array1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
		
		// Do a standard delay embedding with tau 1
		checkMatrix(new double[][] { {4, 3, 2, 1, 0}, {5, 4, 3, 2, 1}, 
				{6, 5, 4, 3, 2}, {7, 6, 5, 4, 3}, 
				{8, 7, 6, 5, 4}, {9, 8, 7, 6, 5} },
				MatrixUtils.makeDelayEmbeddingVector(array1, 5, 4, 6),
				0.00001);
		// Now specify tau explicitly
		checkMatrix(new double[][] { {4, 3, 2, 1, 0}, {5, 4, 3, 2, 1}, 
				{6, 5, 4, 3, 2}, {7, 6, 5, 4, 3}, 
				{8, 7, 6, 5, 4}, {9, 8, 7, 6, 5} },
				MatrixUtils.makeDelayEmbeddingVector(array1, 5, 1, 4, 6),
				0.00001);
		
		// Do same standard delay embedding but starting at an offset
		checkMatrix(new double[][] { {8, 7, 6, 5, 4}, {9, 8, 7, 6, 5} },
				MatrixUtils.makeDelayEmbeddingVector(array1, 5, 8, 2),
				0.00001);
		// Now specify tau explicitly
		checkMatrix(new double[][] { {8, 7, 6, 5, 4}, {9, 8, 7, 6, 5} },
				MatrixUtils.makeDelayEmbeddingVector(array1, 5, 1, 8, 2),
				0.00001);

		// Try with tau 2
		checkMatrix(new double[][] { {8, 6, 4, 2, 0}, {9, 7, 5, 3, 1} },
				MatrixUtils.makeDelayEmbeddingVector(array1, 5, 2, 8, 2),
				0.00001);
		
		// Try with tau 3
		checkMatrix(new double[][] { {6, 3, 0}, {7, 4, 1}, {8, 5, 2}, {9, 6, 3} },
				MatrixUtils.makeDelayEmbeddingVector(array1, 3, 3, 6, 4),
				0.00001);
	}
	
	/**
	 * Check that all entries in the given matrix match those of the expected
	 *  matrix
	 * 
	 * @param expected
	 * @param actual
	 * @param resolution
	 */
	public static void checkMatrix(double[][] expected, double[][] actual, double resolution) {
		for (int r = 0; r < expected.length; r++) {
			for (int c = 0; c < expected[r].length; c++) {
				assertEquals(expected[r][c], actual[r][c], resolution);
			}
		}		
	}

	/**
	 * Check that all entries in the given array match those of the expected
	 *  array
	 * 
	 * @param expected
	 * @param actual
	 */
	public static void checkArray(int[] expected, int[] actual) {
		for (int r = 0; r < expected.length; r++) {
			assertEquals(expected[r], actual[r]);
		}		
	}
	
	public void test2DArrayCopy() {
		double[][] temp = {{1,2,3,4}, {4,5,6,7}, {7,8,9,10}, {10,11,12,13}};
		double[][] newMatrix = new double[10][10];
		
		MatrixUtils.arrayCopy(temp, 1, 1, newMatrix, 3, 3, 3, 3);
		for (int r = 0; r < 3; r++) {
			for (int c = 0; c < 3; c++) {
				assertEquals(temp[r+1][c+1], newMatrix[r+3][c+3]);
			}
		}
	}
	
	public void testAddAndRemoveFromMean() {
		RandomGenerator rg = new RandomGenerator();
		int N = 1000;
		int cols = 3;
		int rowsToRemove = 10;
		double[][] data = rg.generateNormalData(N, cols, 0, 1);
		
		for (int c = 0; c < cols; c++) {
			double mean = MatrixUtils.mean(data, c);
			
			// Remove mean from row c as well
			double meanAfterRemoval = MatrixUtils.removeFromColumnMean(
					data, c, c, rowsToRemove, mean, N);
			
			// And now try manually:
			double sum = mean * (double) N;
			for (int r = c; r < c + rowsToRemove; r++) {
				sum -= data[r][c];
			}
			double meanAfterRemovalManual = sum / (double) (N - rowsToRemove);
			assertEquals(meanAfterRemovalManual, meanAfterRemoval, 1e-5);
			
			// Then add back in:
			double meanAfterAddingBackIn = MatrixUtils.addToColumnMean(
					data, c, c, rowsToRemove, meanAfterRemoval, N - rowsToRemove);
			assertEquals(mean, meanAfterAddingBackIn, 1e-5);
			
			// Take one row out then try swapping it with another:
			double meanAfterOneRemoval = MatrixUtils.removeFromColumnMean(
					data, c, c, 1, mean, N);
			double meanAfterSwappedRemoval = MatrixUtils.swapIntoColumnMean(
					data, c, c, c+1, meanAfterOneRemoval, N-1);
			double meanAfterSwapManual = (mean * (double) N - data[c+1][c]) / ((double) (N - 1));
			assertEquals(meanAfterSwapManual, meanAfterSwappedRemoval, 1e-5);
		}
	}

	public void testAddAndRemoveFromCovariance() {
		RandomGenerator rg = new RandomGenerator();
		
		// Run a short test:
		double[][] shortData = new double[10][1];
		shortData[0][0] = 1; // Only item that contributes to the variance (otherwise it will be zero)
		double[][] shortCovariance = MatrixUtils.covarianceMatrix(shortData);
		double initialCovariance = shortCovariance[0][0];
		System.out.printf("Initial covariance = %.5f\n", initialCovariance);
		MatrixUtils.removeFromCovarianceMatrix(shortCovariance,
				shortData, null, null,
				new int[] {0}, new int[] {}, new int[] {},
				new double[] {0.1}, new double[] {}, new double[] {},
				0, 1, 10);
		// and now check what happens if we remove this row manually:
		double[][] cutShortData = MatrixUtils.selectRows(shortData, 1, 9);
		double[][] covShortAfterManualCut = MatrixUtils.covarianceMatrix(cutShortData);
		double covMatRemoval = MatrixUtils.computeCovarianceMatrixRemoval(shortData, shortData,
				0, 0, 0.1, 0.1, 10, 0, 1);
		System.out.printf("computeCovarianceMatrixRemoval returns %.5f\n", covMatRemoval);
		assertEquals(covShortAfterManualCut[0][0], shortCovariance[0][0], 1e-6);
		assertEquals(0, shortCovariance[0][0], 1e-6);
		// or indeed if and now check what happens if we add this row back in:
		MatrixUtils.addToCovarianceMatrix(shortCovariance,
				shortData, null, null,
				new int[] {0}, new int[] {}, new int[] {},
				new double[] {0}, new double[] {}, new double[] {},
				0, 1, 9);
		assertEquals(initialCovariance, shortCovariance[0][0], 1e-6);
		
		int N = 100;
		int cols = 3;
		int rowsToRemove = 10;
		double[][] data = rg.generateNormalData(N, cols, 0, 1);
		
		double[][] covariances = MatrixUtils.covarianceMatrix(data);
		double[][] covariancesOriginal = MatrixUtils.arrayCopy(covariances);
		MatrixUtils.removeFromCovarianceMatrix(covariances,
				MatrixUtils.selectColumns(data, new int[] {0}), MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {2}),
				new int[] {0}, new int[] {0}, new int[] {0},
				new double[] {MatrixUtils.mean(data, 0)}, new double[] {MatrixUtils.mean(data, 1)}, new double[] {MatrixUtils.mean(data, 2)},
				0, rowsToRemove, N);
		
		// and now check what happens if we remove these rows manually:
		double[][] cutData = MatrixUtils.selectRows(data, rowsToRemove, N-rowsToRemove);
		double[][] covAfterManualCut = MatrixUtils.covarianceMatrix(cutData);
		
		// Need fairly large tolerance here, because the two are being computed in different
		///  ways and so have slightly different numerical errors
		checkMatrix(covAfterManualCut, covariances, 1e-6);
		
		// And check that all is resolved once we insert them again
		double runningMean0 = MatrixUtils.mean(data, 0, rowsToRemove, N-rowsToRemove);
		double runningMean1 = MatrixUtils.mean(data, 1, rowsToRemove, N-rowsToRemove);
		double runningMean2 = MatrixUtils.mean(data, 2, rowsToRemove, N-rowsToRemove);
		MatrixUtils.addToCovarianceMatrix(covariances,
				MatrixUtils.selectColumns(data, new int[] {0}), MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {2}),
				new int[] {0}, new int[] {0}, new int[] {0},
				new double[] {runningMean0},
				new double[] {runningMean1},
				new double[] {runningMean2},
				0, rowsToRemove, N - rowsToRemove);
		// Maybe need larger tolerance here, because the two are being computed in different
		///  ways and so have slightly different numerical errors
		checkMatrix(covariancesOriginal, covariances, 1e-6);
		
		// Now need to test that we can swap rows in and out ok with new code
		List<Integer> rows = MatrixUtils.createArrayList(MatrixUtils.range(rowsToRemove, N-1));
		for (int r = rowsToRemove; r < rowsToRemove + 10; r++) {
			// Swap row r - 1 in, and row r out:
			rows.remove(0); // row r will be the first index for this test
			rows.add(r-1);
			// First from covariance, since it uses the old means:
			MatrixUtils.swapIntoCovarianceMatrix(
				covAfterManualCut,
				MatrixUtils.selectColumns(data, new int[] {0}), MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {2}),
				new int[] {0}, new int[] {0}, new int[] {0},
				new double[] {runningMean0}, new double[] {runningMean1}, new double[] {runningMean2},
				r - 1, r, N-rowsToRemove);
			// Then from the means:
			runningMean0 = MatrixUtils.swapIntoColumnMean(data, 0, r-1, r, runningMean0, N-rowsToRemove);
			runningMean1 = MatrixUtils.swapIntoColumnMean(data, 1, r-1, r, runningMean1, N-rowsToRemove);
			runningMean2 = MatrixUtils.swapIntoColumnMean(data, 2, r-1, r, runningMean2, N-rowsToRemove);
			// Now check them:
			double[][] extractedData = MatrixUtils.selectRows(data, rows);
			double[] means = MatrixUtils.means(extractedData);
			assertEquals(means[0], runningMean0, 1e-6);
			assertEquals(means[1], runningMean1, 1e-6);
			assertEquals(means[2], runningMean2, 1e-6);
			double[][] covarianceMatrix = MatrixUtils.covarianceMatrix(extractedData, means);
			checkMatrix(covarianceMatrix, covAfterManualCut, 1e-6);
		}
	}
}
