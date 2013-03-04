package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateAbstractTester;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;

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
	
}
