package infodynamics.utils;

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
		
		// Now check that it picks up if A is not positive definite:
		double[][] notpositiveDefiniteA = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
		flaggedException = false;
		try {
			MatrixUtils.CholeskyDecomposition(notpositiveDefiniteA);
		} catch (Exception e) {
			flaggedException = true;
		}
		assertTrue(flaggedException);
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
		
		// Check error conditions
	}
	
	/**
	 * Check that all entries in the given matrix match those of the expected
	 *  matrix
	 * 
	 * @param expected
	 * @param actual
	 * @param resolution
	 */
	protected void checkMatrix(double[][] expected, double[][] actual, double resolution) {
		for (int r = 0; r < expected.length; r++) {
			for (int c = 0; c < expected[r].length; c++) {
				assertEquals(expected[r][c], actual[r][c], resolution);
			}
		}		
	}

}
