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

import java.util.Calendar;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

public class ConditionalMutualInfoMultiVariateTester
	extends infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateAbstractTester {

	protected String NUM_THREADS_TO_USE_DEFAULT = ConditionalMutualInfoCalculatorMultiVariateKraskov.USE_ALL_THREADS;
	protected String NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;
	
	/**
	 * Utility function to create a calculator for the given algorithm number
	 * 
	 * @param algNumber
	 * @return
	 */
	public ConditionalMutualInfoCalculatorMultiVariateKraskov getNewCalc(int algNumber) {
		ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc = null;
		if (algNumber == 1) {
			condMiCalc = new ConditionalMutualInfoCalculatorMultiVariateKraskov1();
		} else if (algNumber == 2) {
			condMiCalc = new ConditionalMutualInfoCalculatorMultiVariateKraskov2();
		}
		return condMiCalc;
	}
	
	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void checkLocalsAverageCorrectly(int algNumber) throws Exception {
		
		ConditionalMutualInfoCalculatorMultiVariateKraskov miCalc = getNewCalc(algNumber);
		
		String kraskov_K = "4";
		
		miCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testLocalsAverageCorrectly(miCalc, 2, 100);
	}
	public void testLocalsAverageCorrectly() throws Exception {
		checkLocalsAverageCorrectly(1);
		checkLocalsAverageCorrectly(2);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @throws Exception
	 */
	public void checkComputeSignificanceDoesntAlterAverage(int algNumber) throws Exception {
		
		ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc = getNewCalc(algNumber);
		
		String kraskov_K = "4";
		
		condMiCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testComputeSignificanceDoesntAlterAverage(condMiCalc, 2, 100);
	}
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		checkComputeSignificanceDoesntAlterAverage(1);
		checkComputeSignificanceDoesntAlterAverage(2);
	}
	
	/**
	 * Utility function to run Kraskov conditional MI algorithm 1
	 *  as transfer entropy for data with known results
	 *  from TRENTOOL. (with default parameter settings k=1, l=1)
	 * 
	 * @param var1 source multivariate data set
	 * @param var2 dest multivariate data set
	 * @param kNNs array of Kraskov k nearest neighbours parameter to check
	 * @param expectedResults array of expected results for each k
	 */
	protected void checkTEForGivenData(double[][] var1, double[][] var2,
			int[] kNNs, double[] expectedResults) throws Exception {
		checkTEForGivenData(var1, var2, 1, 1, kNNs, expectedResults);
	}
	
	/**
	 * Utility function to run Kraskov conditional MI algorithm 1
	 *  as transfer entropy for data with known results
	 *  from TRENTOOL.
	 * 
	 * @param var1 source multivariate data set
	 * @param var2 dest multivariate data set
	 * @param historyK history length k of destination
	 * @param historyL history length l of source
	 * @param kNNs array of Kraskov k nearest neighbours parameter to check
	 * @param expectedResults array of expected results for each k
	 */
	protected void checkTEForGivenData(double[][] var1, double[][] var2,
			int historyK, int historyL, int[] kNNs, double[] expectedResults) throws Exception {
				
		ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc = getNewCalc(1);
		
		// Which is the first time index for the dest next state?
		//  It depends on the values of k and l for embedding the past state
		//  of destination and source.
		int firstDestTimeIndex = Math.max(historyK, historyL);
		
		// Normalise the data ourselves rather than letting the calculator do it -
		//  this ensures the extra values in the time series (e.g. last value in source)
		//  are taken into account, in line with TRENTOOL
		var1 = MatrixUtils.normaliseIntoNewArray(var1);
		var2 = MatrixUtils.normaliseIntoNewArray(var2);
		
		for (int kIndex = 0; kIndex < kNNs.length; kIndex++) {
			int k = kNNs[kIndex];
			condMiCalc.setProperty(
					ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
					Integer.toString(k));
			// We already normalised above, and this will do a different
			//  normalisation without taking the extra values in to account if we did it
			condMiCalc.setProperty(
					ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE,
					Boolean.toString(false));
			// No longer need to set this property as it's set by default:
			//condMiCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
			//		EuclideanUtils.NORM_MAX_NORM_STRING);
			condMiCalc.setProperty(
					ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
					NUM_THREADS_TO_USE);
			condMiCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0"); // Need consistency of results for unit test
			condMiCalc.initialise(var1[0].length * historyL,
					var2[0].length, var2[0].length * historyK);
			// Construct the joint vectors of the source states
			double[][] sources = null;
			if (historyL == 1) {
				sources = MatrixUtils.selectRows(var1, firstDestTimeIndex - historyL,
						var1.length - firstDestTimeIndex);
			} else {
				// Build the storage for the source states
				int sourceVars = var1[0].length;
				sources = new double[var1.length - firstDestTimeIndex][sourceVars * historyL];
				for (int t = 0; t < historyL; t++) {
					MatrixUtils.arrayCopy(
							var1, firstDestTimeIndex - historyL + t, 0,
							sources, 0, t*sourceVars,
							var1.length - firstDestTimeIndex, sourceVars);
				}
			}
			// Construct the joint vectors of the conditionals
			double[][] conditionals = null;
			if (historyK == 1) {
				conditionals = MatrixUtils.selectRows(var2, firstDestTimeIndex - historyK,
						var2.length - firstDestTimeIndex);
			} else {
				// Build the storage for the conditional observations
				int destVars = var2[0].length;
				conditionals = new double[var2.length - firstDestTimeIndex][destVars * historyK];
				for (int t = 0; t < historyK; t++) {
					MatrixUtils.arrayCopy(
							var2, firstDestTimeIndex - historyK + t, 0,
							conditionals, 0, t*destVars,
							var2.length - firstDestTimeIndex, destVars);
				}
			}
			// And set the observations using these
			condMiCalc.setObservations(sources,
					MatrixUtils.selectRows(var2, firstDestTimeIndex, var2.length - firstDestTimeIndex),
					conditionals);
			double condMi = condMiCalc.computeAverageLocalOfObservations();
			//miCalc.setDebug(false);
			
			System.out.printf("k=%d: Average MI %.8f (expected %.8f)\n",
					k, condMi, expectedResults[kIndex]);
			// 6 decimal places is Matlab accuracy
			assertEquals(expectedResults[kIndex], condMi, 0.000001);			
		}
	}

	/**
	 * Test the computed univariate TE as a conditional MI
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
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {0.3058006};
		
		System.out.println("Kraskov Cond MI as TE comparison 1 - univariate coupled data 1");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromTRENTOOL);
		
		// And now in the reverse direction:
		expectedFromTRENTOOL = new double[] {-0.0029744};
		
		System.out.println("  reverse direction:");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {0}),
				kNNs, expectedFromTRENTOOL);

	}

	/**
	 * Test the computed univariate TE as a conditional MI
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
	public void testUnivariateTEforCoupledLogisticMapFromFile() throws Exception {
		
		// Test set 1:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/coupledLogisticMapXY.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {0.508417};
		
		System.out.println("Kraskov Cond MI as TE comparison 1 - univariate coupled logistic map data 1");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromTRENTOOL);
		
		// And now in the reverse direction:
		expectedFromTRENTOOL = new double[] {0.016257};
		
		System.out.println("  reverse direction:");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {0}),
				kNNs, expectedFromTRENTOOL);

	}

	/**
	 * Test the computed univariate TE as a conditional MI
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
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {-0.0096556};
		
		System.out.println("Kraskov Cond MI as TE comparison 2 - univariate random data 1 (col 0->1)");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromTRENTOOL);
		
		// And now for other columns
		expectedFromTRENTOOL = new double[] {0.0175389};
		
		System.out.println("  (col 1->2):");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {2}),
				kNNs, expectedFromTRENTOOL);

		// And now for other columns
		expectedFromTRENTOOL = new double[] {0.0026367};
		
		System.out.println("  (col 1->0):");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {0}),
				kNNs, expectedFromTRENTOOL);

		// And now for other columns
		expectedFromTRENTOOL = new double[] {-0.00012474};
		
		System.out.println("  (col 0->2):");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {2}),
				kNNs, expectedFromTRENTOOL);

		// And now for other columns
		expectedFromTRENTOOL = new double[] {-5.4437e-03};
		
		System.out.println("  (col 2->0):");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {2}),
				MatrixUtils.selectColumns(data, new int[] {0}),
				kNNs, expectedFromTRENTOOL);
	}

	/**
	 * Test the computed multivariate conditional MI
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * It's multivariate because we use embedding dimension 2 on both source
	 *  and destination.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 2, 2, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=2,
	 * history embedding dimension l=2 on source as well.
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMultivariateCondMIforCoupledDataFromFile() throws Exception {
		
		// Test set 3:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/4ColsPairedOneStepNoisyDependence-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {0.1400645};
		
		System.out.println("Kraskov Cond MI as TE - multivariate coupled data 1, k=2,l=2");
		System.out.println("  (0->2)");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {2}),
				2, 2,
				kNNs, expectedFromTRENTOOL);
		
		// And now for reverse direction:
		expectedFromTRENTOOL = new double[] {-0.0181459};
		
		System.out.println("  (2->0):");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {2}),
				MatrixUtils.selectColumns(data, new int[] {0}),
				2, 2,
				kNNs, expectedFromTRENTOOL);		

		// And now for other columns:
		expectedFromTRENTOOL = new double[] {0.1639186};
		System.out.println("  (1->3):");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {3}),
				2, 2,
				kNNs, expectedFromTRENTOOL);		
		// And in reverse:
		expectedFromTRENTOOL = new double[] {0.0036976};
		System.out.println("  (3->1):");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {3}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				2, 2,
				kNNs, expectedFromTRENTOOL);
		
		// -------------
		// And finally, confirm that we get different results for k=1,l=1,
		//  which match TRENTOOL
		expectedFromTRENTOOL = new double[] {0.0072169};
		System.out.println("  (0->1) but with k=1,l=1:");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				1, 1,
				kNNs, expectedFromTRENTOOL);		
		// And in reverse
		expectedFromTRENTOOL = new double[] {0.0011738};
		System.out.println("  (1->2) but with k=1,l=1:");
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {1}),
				MatrixUtils.selectColumns(data, new int[] {2}),
				1, 1,
				kNNs, expectedFromTRENTOOL);		
	}

	/**
	 * Unit test for conditional MI on new observations.
	 * We can test this against the calculator itself. If we send in the original data set as new observations,
	 *  we can recreate the neighbour counts (plus one) by setting K to 1 larger (to account for the data point itself), and
	 *  account for the change in bias. 
	 * 
	 * @throws Exception
	 */
	public void testMultivariateCondMIForNewObservations() throws Exception {
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/4ColsPairedOneStepNoisyDependence-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		// RandomGenerator rg = new RandomGenerator();
		//double[][] data = rg.generateNormalData(50, 4, 0, 1);

		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		
		System.out.println("Kraskov Cond MI testing new Observations:");
		
		for (int ki = 0; ki < kNNs.length; ki++) {
			ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc = getNewCalc(1);
			ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalcForNew = getNewCalc(1);
			// Let it normalise by default
			// And no noise addition to protect the integrity of our neighbour counts under both techniques here:
			condMiCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0");
			condMiCalcForNew.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0");
			double[][] var1 = MatrixUtils.selectColumns(data, new int[] {0});
			double[][] var2 = MatrixUtils.selectColumns(data, new int[] {1});
			double[][] condVar = MatrixUtils.selectColumns(data, new int[] {2, 3});
			
			// Compute MI(0;1|2,3) :
			condMiCalc.setProperty(
					ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
					Integer.toString(kNNs[ki]));
			System.out.println("Main calc normalisation is " + condMiCalc.getProperty(condMiCalc.PROP_NORMALISE));
			condMiCalc.initialise(var1[0].length, var2[0].length, condVar[0].length);
			condMiCalc.setObservations(var1, var2, condVar);
			double miAverage = condMiCalc.computeAverageLocalOfObservations();
			// Now compute as new observations:
			condMiCalcForNew.setProperty(
					ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
					Integer.toString(kNNs[ki] + 1)); // Using K = K + 1
			// condMiCalcForNew.setProperty(
			//		ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
			//		"1");
			System.out.println("New obs calc normalisation is " + condMiCalcForNew.getProperty(condMiCalc.PROP_NORMALISE));
			condMiCalcForNew.initialise(var1[0].length, var2[0].length, condVar[0].length);
			condMiCalcForNew.setObservations(var1, var2, condVar);
			// condMiCalc.setDebug(true);
			//condMiCalcForNew.setDebug(true);
			double[] newLocals = condMiCalcForNew.computeLocalUsingPreviousObservations(var1, var2, condVar);
			//condMiCalcForNew.setDebug(false);
			double averageFromNewObservations = MatrixUtils.mean(newLocals);
			// We can't check this directly, so test each point individually:
			for (int t = 0; t < data.length; t++) {
				double[] originalNeighbourCounts = condMiCalc.partialComputeFromObservations(t, 1, false);
				// Need to normalise the data before passing it in here -- this is
				//  what is happening inside computeLocalUsingPreviousObservations above
				double[] newObsNeighbourCounts = condMiCalcForNew.partialComputeFromNewObservations(
						t, 1,
						MatrixUtils.normaliseIntoNewArray(var1),
						MatrixUtils.normaliseIntoNewArray(var2),
						MatrixUtils.normaliseIntoNewArray(condVar), false);
				// Now check each return count in the array:
				if (originalNeighbourCounts[1] != newObsNeighbourCounts[1] - 1) {
					System.out.println("Assertion failure for t=" + t + ": expected " + originalNeighbourCounts[1] +
							" from original, plus 1, but got " + newObsNeighbourCounts[1]);
					System.out.print("Actual raw data was: ");
					MatrixUtils.printArray(System.out, data[0]);
				}
				assertEquals(originalNeighbourCounts[1], newObsNeighbourCounts[1] - 1); // Nxz should be 1 higher
				assertEquals(originalNeighbourCounts[2], newObsNeighbourCounts[2] - 1); // Nyz should be 1 higher
				assertEquals(originalNeighbourCounts[3], newObsNeighbourCounts[3] - 1); // Nz should be 1 higher
				// Now check the local value at each point using these verified counts:
				double newLocalValue = condMiCalcForNew.digammaK -
						MathsUtils.digamma((int) newObsNeighbourCounts[1] + 1) -
						MathsUtils.digamma((int) newObsNeighbourCounts[2] + 1) +
						MathsUtils.digamma((int) newObsNeighbourCounts[3] + 1);
				/* System.out.printf("n=%d, n_xz=%d, n_yz=%d, n_z=%d, local=%.4f," +
						" digamma(n_xz+1)=%.5f, digamma(n_yz+1)=%.5f, digamma(n_z+1)=%.5f, \n",
						t, (int) newObsNeighbourCounts[1], (int) newObsNeighbourCounts[2],
						(int) newObsNeighbourCounts[3], newLocalValue,
						MathsUtils.digamma((int) newObsNeighbourCounts[1] + 1), 
						MathsUtils.digamma((int) newObsNeighbourCounts[2] + 1),
						MathsUtils.digamma((int) newObsNeighbourCounts[3] + 1)); */
				if (Math.abs(newLocalValue - newLocals[t]) > 0.00000001) {
					System.out.printf("t=%d: Assertion failed: computed local was %.5f, local from nn counts was %.5f\n",
							t, newLocals[t], newLocalValue);
				}
				assertEquals(newLocalValue, newLocals[t], 0.00000001);
			}
			
			// And check that the two methods are equivalent, after correcting for the different K biases:
			// CAN'T EXPECT THIS TO WORK though, because different neighbour counts (by one in each case) puts
			//  different inputs to the digamma functions, and we have no way to check that for the average here;
			//  Our above checks on the counts for each point suffice though.
			// assertEquals(miAverage, averageFromNewObservations - MathsUtils.digamma(kNNs[ki]+1) + MathsUtils.digamma(kNNs[ki]),
			//		0.00000001);
		}
	}

	/**
	 * Test the computed univariate TE as a conditional MI
	 * using various numbers of threads.
	 * 
	 * Test the computed univariate TE as a conditional MI
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
	public void testUnivariateTEVariousNumberThreads() throws Exception {
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {-0.0096556};
		
		System.out.println("Kraskov Cond MI as TE - multivariate coupled data 1, k=2,l=2 (0->2)");
		System.out.println(" with various numbers of threads:");
		System.out.println(" -- 1 thread:");
		NUM_THREADS_TO_USE = "1";
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromTRENTOOL);
		System.out.println(" -- 2 threads:");
		NUM_THREADS_TO_USE = "2";
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromTRENTOOL);
		System.out.println(" -- 3 threads:");
		NUM_THREADS_TO_USE = "3";
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromTRENTOOL);
		System.out.println(" -- all threads:");
		NUM_THREADS_TO_USE = ConditionalMutualInfoCalculatorMultiVariateKraskov.USE_ALL_THREADS;
		checkTEForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromTRENTOOL);
		
		// And finally test that multithreading is still ok if we have
		//  an imbalanced number of data between each thread.
		// Expected value is only generated from our own code; we're not so 
		//  interested in checking for this precise value, as we are in
		//  checking that the value is stable when we change the number of
		//  threads and have an uneven amount of data in each thread.
		double[] expectedValue = new double[] {0.026517704};
		NUM_THREADS_TO_USE = "2";
		checkTEForGivenData(
				MatrixUtils.selectRows(
						MatrixUtils.selectColumns(data, new int[] {0}),
						0, 501),
				MatrixUtils.selectRows(
						MatrixUtils.selectColumns(data, new int[] {1}),
						0, 501),
				kNNs, expectedValue);
		NUM_THREADS_TO_USE = "3";
		checkTEForGivenData(
				MatrixUtils.selectRows(
						MatrixUtils.selectColumns(data, new int[] {0}),
						0, 501),
				MatrixUtils.selectRows(
						MatrixUtils.selectColumns(data, new int[] {1}),
						0, 501),
				kNNs, expectedValue);
		NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;
	}
	
	/**
	 * Test the conditional MI without any conditionals
	 * against the MI calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIxnyn <dataFile> 1 1 3000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testNoConditionalsIsMIforRandomVariablesFromFile() throws Exception {
		
		// Test set 1:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/2randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {1, 2, 3, 4, 5, 6, 10, 15};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA = {-0.05294175, -0.03944338, -0.02190217,
				0.00120807, -0.00924771, -0.00316402, -0.00778205, -0.00565778};
		
		System.out.println("Kraskov comparison 1 - univariate random data 1 --");
		// Check for algorithm 1 first (no expected results from MILCA included):
		System.out.println("Algorithm 1:");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, null);
		// then check for algorithm 2 with expected results from MILCA
		System.out.println("Algorithm 2:");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromMILCA);
		
		//------------------
		// Test set 2:
		
		// We'll just take the first two columns from this data set
		afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		data = afr.getDouble2DMatrix();
		
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA_2 = {-0.04614525, -0.00861460, -0.00164540,
				-0.01130354, -0.01339670, -0.00964035, -0.00237072, -0.00096891};
		
		System.out.println("Kraskov comparison 2 - univariate random data 2 --");
		// Check for algorithm 1 first (no expected results from MILCA included):
		System.out.println("Algorithm 1:");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, null);
		// then check for algorithm 2 with expected results from MILCA
		System.out.println("Algorithm 2:");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromMILCA_2);

	}

	/**
	 * Utility function to run Kraskov MI for data with known results
	 * 
	 * @param var1
	 * @param var2
	 * @param kNNs array of Kraskov k nearest neighbours parameter to check
	 * @param expectedResults array of expected results for each k
	 */
	protected void checkMIForGivenData(double[][] var1, double[][] var2,
			int[] kNNs, double[] expectedResults) throws Exception {
				
		// The Kraskov MILCA toolkit MIhigherdim executable 
		//  uses algorithm 2 by default (this is what it means by rectangular);
		// These are the only ones we have expected results for
		ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc = getNewCalc(expectedResults == null ? 1 : 2);
		
		for (int kIndex = 0; kIndex < kNNs.length; kIndex++) {
			int k = kNNs[kIndex];
			condMiCalc.setProperty(
					MutualInfoCalculatorMultiVariateKraskov.PROP_K,
					Integer.toString(k));
			condMiCalc.setProperty(
					MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
					NUM_THREADS_TO_USE);
			// No longer need to set this property as it's set by default:
			//miCalc.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
			//		EuclideanUtils.NORM_MAX_NORM_STRING);
			condMiCalc.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0"); // Need consistency for unit tests
			condMiCalc.initialise(var1[0].length, var2[0].length, 0);
			condMiCalc.setObservations(var1, var2, null);
			// condMiCalc.setDebug(true);
			double mi = condMiCalc.computeAverageLocalOfObservations();
			// condMiCalc.setDebug(false);
			
			// And compute also passing in an array of vectors of length 0
			condMiCalc.initialise(var1[0].length, var2[0].length, 0);
			condMiCalc.setObservations(var1, var2, new double[var1.length][]);
			// condMiCalc.setDebug(true);
			double miFromEmptyVectors = condMiCalc.computeAverageLocalOfObservations();
			
			double expectedMi = 0.0;
			if (expectedResults == null) {
				// We don't have expected results from MILCA toolkit for algorithm 1.
				// In this case, we'll generate our own expectations:
				MutualInfoCalculatorMultiVariateKraskov1 miCalc = new MutualInfoCalculatorMultiVariateKraskov1();
				miCalc.setProperty(
						MutualInfoCalculatorMultiVariateKraskov.PROP_K,
						Integer.toString(k));
				miCalc.setProperty(
						MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
						NUM_THREADS_TO_USE);
				miCalc.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0"); // Need consistency for unit tests
				miCalc.initialise(var1[0].length, var2[0].length);
				miCalc.setObservations(var1, var2);
				expectedMi = miCalc.computeAverageLocalOfObservations();
			} else {
				expectedMi = expectedResults[kIndex];
			}
			System.out.printf("k=%d: Average MI %.8f (expected %.8f)\n",
					k, mi, expectedMi);
			// Dropping required accuracy by one order of magnitude, due
			//  to faster but slightly less accurate digamma estimator change
			assertEquals(expectedMi, mi, 0.0000001);
			assertEquals(expectedMi, miFromEmptyVectors, 0.0000001);
		}
	}
}
