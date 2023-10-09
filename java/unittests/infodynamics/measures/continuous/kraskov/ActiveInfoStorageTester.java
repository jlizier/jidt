package infodynamics.measures.continuous.kraskov;

import java.util.Arrays;

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
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
		// Select data points 2350:3550
		data = MatrixUtils.selectRows(data, 2349, 3550-2350+1);

		ActiveInfoStorageCalculatorKraskov ais = new ActiveInfoStorageCalculatorKraskov();
		// ais.setDebug(true);
		
		// Use one thread to test first:
		ais.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS, "1");
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ);
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_K_SEARCH_MAX, "4");
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_TAU_SEARCH_MAX, "2");
		ais.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0"); // Need consistency of results for unit test
		ais.initialise();
		ais.setObservations(MatrixUtils.selectColumn(data, 0));
		int optimisedK = Integer.parseInt(ais.getProperty(ActiveInfoStorageCalculator.K_PROP_NAME));
		int optimisedTau = Integer.parseInt(ais.getProperty(ActiveInfoStorageCalculator.TAU_PROP_NAME));
		double aisOptimisedSingleThread = ais.computeAverageLocalOfObservations();
		System.out.println("AIS was " + aisOptimisedSingleThread + " for k=" + optimisedK +
				" and tau=" + optimisedTau + " optimised over kNNs=" +
				ais.getProperty(ActiveInfoStorageCalculatorKraskov.PROP_RAGWITZ_NUM_NNS));
		
		// Test that the answer was k=2, tau=1 for this data set (k=3 for full data set)
		assertEquals(2, optimisedK);
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
		ais.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0"); // Need consistency of results for unit test
		ais.initialise(optimisedK, optimisedTau);
		ais.setObservations(MatrixUtils.selectColumn(data, 0));
		double aisManualParamSetting = ais.computeAverageLocalOfObservations();
		assertEquals(aisOptimisedSingleThread, aisManualParamSetting, 0.00000001);
		
		// Test that it works if we supply a validity vector:
		ais = new ActiveInfoStorageCalculatorKraskov();
		ais.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0"); // Need consistency of results for unit test
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ);
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_K_SEARCH_MAX, "4");
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_TAU_SEARCH_MAX, "2");
		ais.initialise();
		boolean[] validity = new boolean[data.length];
		Arrays.fill(validity, true);
		ais.setObservations(MatrixUtils.selectColumn(data, 0), validity);
		double aisWithValidity = ais.computeAverageLocalOfObservations();
		assertEquals(aisOptimisedSingleThread, aisWithValidity, 0.00000001);
		assertEquals(optimisedK, Integer.parseInt(ais.getProperty(ActiveInfoStorageCalculator.K_PROP_NAME)));
		assertEquals(optimisedTau, Integer.parseInt(ais.getProperty(ActiveInfoStorageCalculator.TAU_PROP_NAME)));

		// Finally, test that we can use a different number of kNNs to the MI calculator
		ais = new ActiveInfoStorageCalculatorKraskov();
		ais.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE, "0"); // Need consistency of results for unit test
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_RAGWITZ_NUM_NNS, "10");
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ);
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_K_SEARCH_MAX, "4");
		ais.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_TAU_SEARCH_MAX, "4");
		ais.initialise();
		ais.setObservations(MatrixUtils.selectColumn(data, 0));
		double differentNNResult = ais.computeAverageLocalOfObservations();
		assertEquals(10, Integer.parseInt(ais.getProperty(ActiveInfoStorageCalculatorKraskov.PROP_RAGWITZ_NUM_NNS)));
		System.out.printf("Confirmed that we can change the number of nearest neighbours for Ragwitz optimisation, using " +
				"10 neighbours we get k=%s, tau=%s, ais=%.3f\n", ais.getProperty(ActiveInfoStorageCalculator.K_PROP_NAME),
				ais.getProperty(ActiveInfoStorageCalculator.TAU_PROP_NAME), differentNNResult);
	}
	
	/**
	 * Test that observationSetIndices and observationStartTimePoints 
	 * for the underlying MI estimator are written properly
	 * 
	 * @throws Exception 
	 */
	@SuppressWarnings("static-access")
	public void testObservationSetIndices() throws Exception {
		
		int timeSteps = 1000;
		
		ActiveInfoStorageCalculatorKraskov aisCalc = 
				new ActiveInfoStorageCalculatorKraskov();
		RandomGenerator rg = new RandomGenerator();
		double[] data = rg.generateNormalData(timeSteps, 0, 1);
		aisCalc.initialise();
		
		// First check that for a simple single observation set everything works:
		aisCalc.setObservations(data);

		int[] observationSetIds = aisCalc.getObservationSetIndices();
		int[] timeSeriesIndices = aisCalc.getObservationTimePoints();
		assert(observationSetIds.length == timeSteps - 1);
		for (int t = 0; t < timeSteps - 1; t++) {
			assertEquals(observationSetIds[t], 0);
			assertEquals(timeSeriesIndices[t], 1 + t); // Should be offset by 1 for k history
		}
		
		// Now add the same one twice:
		aisCalc.initialise();
		aisCalc.startAddObservations();
		aisCalc.addObservations(data);
		aisCalc.addObservations(data);
		aisCalc.finaliseAddObservations();
		observationSetIds = aisCalc.getObservationSetIndices();
		timeSeriesIndices = aisCalc.getObservationTimePoints();
		assert(observationSetIds.length == 2*timeSteps - 2);
		for (int t = 0; t < timeSteps - 1; t++) {
			assertEquals(observationSetIds[t], 0);
			assertEquals(timeSeriesIndices[t], 1 + t);
		}
		for (int t = 0; t < timeSteps - 1; t++) {
			assertEquals(observationSetIds[timeSteps - 1 + t], 1);
			assertEquals(timeSeriesIndices[timeSteps - 1 + t], 1 + t);
		}
		
		// Now add NUM_SEGMENTS segments:
		int NUM_SEGMENTS = 10;
		aisCalc.setProperty(aisCalc.K_PROP_NAME, "3");
		aisCalc.setProperty(aisCalc.TAU_PROP_NAME, "3");
		aisCalc.initialise();
		aisCalc.startAddObservations();
		for (int r = 0; r < NUM_SEGMENTS; r++) {
			aisCalc.addObservations(rg.generateNormalData(timeSteps, 0, 1));
		}
		aisCalc.finaliseAddObservations();
		observationSetIds = aisCalc.getObservationSetIndices();
		timeSeriesIndices = aisCalc.getObservationTimePoints();
		int k = Integer.valueOf(aisCalc.getProperty(aisCalc.K_PROP_NAME)); // In case we change to autoembedding above
		int tau = Integer.valueOf(aisCalc.getProperty(aisCalc.TAU_PROP_NAME));
		System.out.printf("Embedding dimension %d and delay %d\n", k, tau);
		assert(observationSetIds.length == (timeSteps - (k - 1)*tau - 1) * NUM_SEGMENTS);
		int t = 0;
		for (int r = 0; r < NUM_SEGMENTS; r++) {
			for (int i = 0; i < timeSteps - (k-1)*tau - 1; i++) {
				assertEquals(observationSetIds[t], r);
				assertEquals(timeSeriesIndices[t], (k-1)*tau + 1 + i);
				t++;
			}
		}
	}}
