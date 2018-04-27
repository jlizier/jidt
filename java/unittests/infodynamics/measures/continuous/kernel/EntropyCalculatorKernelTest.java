package infodynamics.measures.continuous.kernel;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class EntropyCalculatorKernelTest extends TestCase {

	public void testLocalEntropiesAverage() throws Exception {
		RandomGenerator rg = new RandomGenerator();
		double[] data = rg.generateNormalData(100, 0, 1);
		EntropyCalculatorKernel entcalc = new EntropyCalculatorKernel();
		entcalc.initialise(0.5);
		entcalc.setObservations(data);
		double entropy = entcalc.computeAverageLocalOfObservations();
		double[] localEntropies = entcalc.computeLocalOfPreviousObservations();
		double avgLocal = MatrixUtils.mean(localEntropies);
		// There are many sources of numerical noise in the combination of so many
		//  local entropy calculations here, so we need to leave a wide tolerance
		assertEquals(entropy, avgLocal, 0.000001);
	}
}
