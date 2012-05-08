package infodynamics.measures.continuous.kraskov;

/**
 * <p>Compute the Mutual Information between two vectors using the Kraskov estimation method.
 * Computes this using the multi-info (or integration) in the marginal spaces, using algorithm 2</p>
 * @see "Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P., Physical Review E 69, (2004) 066138
 * @see http://dx.doi.org/10.1103/PhysRevE.69.066138
 * 
 * @author Joseph Lizier
 */
public class MutualInfoCalculatorMultiVariateKraskovByMulti2 extends
		MutualInfoCalculatorMultiVariateKraskovByMulti {

	/**
	 * 
	 */
	public MutualInfoCalculatorMultiVariateKraskovByMulti2() {
		super();
	}

	@Override
	protected void createMultiInfoCalculators() {
		multiInfoJoint = new MultiInfoCalculatorKraskov2();
		multiInfo1 = new MultiInfoCalculatorKraskov2();
		multiInfo2 =  new MultiInfoCalculatorKraskov2();
	
	}
	
}
