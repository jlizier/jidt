package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.ActiveInfoStorageCalculatorViaMutualInfo;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class ActiveInfoStorageMultiVariateGaussianTester extends TestCase {

	public void testBasicFunctionality() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();

		int[] ks = {1,2,3,4};
		int[] taus = {1,2,3,4};
		
		for (int kIndex = 0; kIndex < ks.length; kIndex++) {
			int k = ks[kIndex];
			for (int tauIndex = 0; tauIndex < taus.length; tauIndex++) {
				int tau = taus[tauIndex];
		
				ActiveInfoStorageCalculatorMultiVariateGaussian ais =
						new ActiveInfoStorageCalculatorMultiVariateGaussian();
				ais.initialise(data[0].length, k, tau);
				ais.setObservations(data);

				assertEquals(data.length -  (k-1)*tau - 1,
						ais.getNumObservations());
				
				@SuppressWarnings("unused")
				double result = ais.computeAverageLocalOfObservations();
			}
		}
	}

	public void testHandlesOneDimension() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		double[] column1 = MatrixUtils.selectColumn(data, 0);
		double[][] column1In2D = MatrixUtils.selectColumns(data, new int[] {0});

		int[] ks = {1,2,3,4};
		int[] taus = {1,2,3,4};
		
		for (int kIndex = 0; kIndex < ks.length; kIndex++) {
			int k = ks[kIndex];
			for (int tauIndex = 0; tauIndex < taus.length; tauIndex++) {
				int tau = taus[tauIndex];
		
				ActiveInfoStorageCalculatorMultiVariateGaussian ais =
						new ActiveInfoStorageCalculatorMultiVariateGaussian();
				ais.initialise(1, k, tau);
				ais.setObservations(column1);
				assertEquals(data.length -  (k-1)*tau - 1,
						ais.getNumObservations());
				double resultFromColumn1 = ais.computeAverageLocalOfObservations();
				
				ais.initialise(1, k, tau);
				ais.setObservations(column1In2D);
				assertEquals(data.length -  (k-1)*tau - 1,
						ais.getNumObservations());
				double resultFrom2D = ais.computeAverageLocalOfObservations();
				
				assertEquals(resultFromColumn1, resultFrom2D, 0.00000001);
			}
		}
	}

	public void testAutoEmbeddingAIS() throws Exception {
		System.out.println("Start AIS multivariate autoembedding test (Gaussian).");
		
		checkAutoEmbeddingLength(false, false);
		checkAutoEmbeddingLength(false, true);
		checkAutoEmbeddingLength(true, false);
		checkAutoEmbeddingLength(true, true);
	}
	
	public void checkAutoEmbeddingLength(boolean source1Has3, boolean source2Has3) throws Exception {
		
		// Generate multivariate data
		int dataLength = 5000;
		RandomGenerator rg = new RandomGenerator();
		double[][] data = rg.generateNormalData(dataLength, 2, 0, 1);
		
		for (int i=3; i < dataLength; i++) {
			data[i][0] = 0.2*data[i-2][0] + 0.2*data[i-1][0] + data[i][0];
			if (source1Has3) {
				data[i][0] += 0.2*data[i-3][0];
			}
			data[i][1] = 0.2*data[i-2][1] + 0.2*data[i-1][1] + data[i][1];
			if (source2Has3) {
				data[i][1] += 0.2*data[i-3][1];
			}
		}
		
		int correctK = (source1Has3 || source2Has3) ? 3 : 2;
		
		// Instantiate calculator and set search bounds
		ActiveInfoStorageCalculatorMultiVariateGaussian aisCalc = 
					new ActiveInfoStorageCalculatorMultiVariateGaussian();
		// Can't search larger than the max we expect because fluctuations can make the selection larger
		aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_AUTO_EMBED_METHOD,
				ActiveInfoStorageCalculatorViaMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS);
		// We only check up to the correct k -- selections beyond this are stochastic, so the algorithm
		//  could choose a larger value if we let it, and that's not wrong.
		aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_K_SEARCH_MAX, Integer.toString(correctK));
		aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_TAU_SEARCH_MAX, "1");
		aisCalc.setDebug(true);
		
		// Run optimisation
		aisCalc.initialise(2);
		aisCalc.setObservations(data);
		int optimisedK = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculatorViaMutualInfo.K_PROP_NAME));
		
		// Test that answer was correct
		assertEquals(correctK, optimisedK);
		
		double estimate = aisCalc.computeAverageLocalOfObservations();
		System.out.printf("Computed average was %.5f from %d samples (embedding length optimisedK = %d\n",
				estimate, aisCalc.getNumObservations(), optimisedK);
	}
}
