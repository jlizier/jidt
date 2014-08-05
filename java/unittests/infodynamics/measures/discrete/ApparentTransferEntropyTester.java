package infodynamics.measures.discrete;

import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class ApparentTransferEntropyTester extends TestCase {

	protected RandomGenerator rand = new RandomGenerator();

	public void testNoInfoForPredictableDest() {
		int[] source = rand.generateRandomInts(100, 2);
		int[] dest = source;
		
		TransferEntropyCalculator teCalc =
				new TransferEntropyCalculator(2, 1);
		teCalc.initialise();
		teCalc.addObservations(dest, source);
		double result = teCalc.computeAverageLocalOfObservations();
		assertEquals(0, result, 0.0000001);
	}

	public void testFullyPredictableVariousEmbeddings() {
		int[] source = new int[101];
		int[] dest = new int[101];
		
		for (int t = 0; t < source.length; t++) {
			int phase = t % 4;
			if (phase < 2) {
				dest[t] = 1;
			}
			if (phase % 2 == 0) {
				source[t] = 1;
			}
		}
		
		// with k=1, we won't see the self-predictability of the dest
		//  and will assume all info is in TE
		TransferEntropyCalculator teCalc =
				new TransferEntropyCalculator(2, 1);
		teCalc.initialise();
		teCalc.addObservations(dest, source);
		double result = teCalc.computeAverageLocalOfObservations();
		assertEquals(1, result, 0.0000001);
		
		// with k=2, we will see the self-predictability of the dest
		teCalc = new TransferEntropyCalculator(2, 2);
		teCalc.initialise();
		teCalc.addObservations(dest, source);
		result = teCalc.computeAverageLocalOfObservations();
		assertEquals(0, result, 0.0000001);

		// with k=1, we won't see the self-predictability of the dest
		//  but if we supply the same dest as the source and use a large
		//  source embedding length then we should see it.
		//  Create dest again to ensure we have balanced observations:
		int[] dest2 = new int[102];
		for (int t = 0; t < dest2.length; t++) {
			int phase = t % 4;
			if (phase < 2) {
				dest2[t] = 1;
			}
		}
		teCalc = new TransferEntropyCalculator(2, 1, 2);
		teCalc.initialise();
		teCalc.addObservations(dest2, dest2);
		result = teCalc.computeAverageLocalOfObservations();
		assertEquals(1, result, 0.0000001);
	}
}
