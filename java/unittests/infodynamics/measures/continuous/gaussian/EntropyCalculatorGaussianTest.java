package infodynamics.measures.continuous.gaussian;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class EntropyCalculatorGaussianTest extends TestCase {

	public void testVarianceSetting() {
		EntropyCalculatorGaussian entcalc = new EntropyCalculatorGaussian();
		entcalc.initialise();
		double variance = 3.45567;
		entcalc.setVariance(variance);
		assertEquals(variance, entcalc.variance);
	}

	public void testVarianceCalculation() {
		EntropyCalculatorGaussian entcalc = new EntropyCalculatorGaussian();
		entcalc.initialise();
		RandomGenerator randomGenerator = new RandomGenerator();
		double[] gaussianObservations = randomGenerator.generateNormalData(100, 5, 3);
		double variance = MatrixUtils.stdDev(gaussianObservations);
		variance *= variance;
		entcalc.setObservations(gaussianObservations);
		assertEquals(variance, entcalc.variance);
	}
	
	public void testEntropyCalculation() {
		EntropyCalculatorGaussian entcalc = new EntropyCalculatorGaussian();
		entcalc.initialise();
		double variance = 3.45567;
		entcalc.setVariance(variance);
		double expectedEntropy = 0.5 * Math.log(2.0*Math.PI*Math.E*variance);
		assertEquals(expectedEntropy, entcalc.computeAverageLocalOfObservations());
	}
}
