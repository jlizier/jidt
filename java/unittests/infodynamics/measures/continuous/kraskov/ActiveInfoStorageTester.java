package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import junit.framework.TestCase;

public class ActiveInfoStorageTester extends TestCase {

	protected String NUM_THREADS_TO_USE_DEFAULT = MutualInfoCalculatorMultiVariateKraskov.USE_ALL_THREADS;
	protected String NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;
	
	public void testEmbedding() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();

		int[] ks = {1,2,3,4};
		int[] taus = {1,2,3,4};
		
		for (int kIndex = 0; kIndex < ks.length; kIndex++) {
			int k = ks[kIndex];
			for (int tauIndex = 0; tauIndex < taus.length; tauIndex++) {
				int tau = taus[tauIndex];
		
				ActiveInfoStorageCalculatorKraskov ais = new ActiveInfoStorageCalculatorKraskov();
				ais.initialise(k, tau);
				ais.setObservations(MatrixUtils.selectColumn(data, 0));

				assertEquals(data.length -  (k-1)*tau - 1,
						ais.getNumObservations());
			}
		}
	}
	
	public void testAutoEmbeddingRagwitz() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/SFI-heartRate_breathVol_bloodOx.txt");
		double[][] data = afr.getDouble2DMatrix();

		ActiveInfoStorageCalculatorKraskov ais = new ActiveInfoStorageCalculatorKraskov();
		// ais.setDebug(true);
		
		// Use one thread to test first:
		ais.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS, "1");
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ);
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_K_SEARCH_MAX, "4");
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_TAU_SEARCH_MAX, "2");
		ais.initialise();
		ais.setObservations(MatrixUtils.selectColumn(data, 0));
		int optimisedK = Integer.parseInt(ais.getProperty(ActiveInfoStorageCalculator.K_PROP_NAME));
		int optimisedTau = Integer.parseInt(ais.getProperty(ActiveInfoStorageCalculator.TAU_PROP_NAME));
		double aisOptimisedSingleThread = ais.computeAverageLocalOfObservations();
		System.out.println("AIS was " + aisOptimisedSingleThread + " for k=" + optimisedK +
				" and tau=" + optimisedTau + " optimised over kNNs=" +
				ais.getProperty(ActiveInfoStorageCalculatorKraskov.PROP_RAGWITZ_NUM_NNS));
		
		// Test that the answer was k=3, tau=1 for this data set
		assertEquals(3, optimisedK);
		assertEquals(1, optimisedTau);
		// Test that kNNs are equal to that used by the MI calculator when we have not set this
		assertEquals(ais.getProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_K),
				ais.getProperty(ActiveInfoStorageCalculatorKraskov.PROP_RAGWITZ_NUM_NNS));
		
		// Test that we get the same answer by a multi-threaded approach
		// Use one thread to test first:
		ais.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS, 
				MutualInfoCalculatorMultiVariateKraskov.USE_ALL_THREADS);
		ais.initialise();
		ais.setObservations(MatrixUtils.selectColumn(data, 0));
		double aisOptimisedMultiThread = ais.computeAverageLocalOfObservations();
		assertEquals(aisOptimisedSingleThread, aisOptimisedMultiThread, 0.00000001);

		// Test that we get the same answer by setting these parameters
		ais = new ActiveInfoStorageCalculatorKraskov();
		ais.initialise(optimisedK, optimisedTau);
		ais.setObservations(MatrixUtils.selectColumn(data, 0));
		double aisManualParamSetting = ais.computeAverageLocalOfObservations();
		assertEquals(aisOptimisedSingleThread, aisManualParamSetting, 0.00000001);
		
		// Finally, test that we can use a different number of kNNs to the MI calculator
		ais = new ActiveInfoStorageCalculatorKraskov();
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_RAGWITZ_NUM_NNS, "8");
		ais.initialise(optimisedK, optimisedTau);
		ais.setObservations(MatrixUtils.selectColumn(data, 0));
		ais.computeAverageLocalOfObservations();
		assertEquals(8, Integer.parseInt(ais.getProperty(ActiveInfoStorageCalculatorKraskov.PROP_RAGWITZ_NUM_NNS)));
	}
	
}
