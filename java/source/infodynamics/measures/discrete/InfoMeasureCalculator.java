package infodynamics.measures.discrete;

/**
 * Info theoretic measure calculator base class
 * 
 * Usage:
 * 1. Continuous accumulation of observations before computing :
 *   Call: a. initialise()
 *         b. addObservations() several times over
 *         c. computeLocalFromPreviousObservations() or computeAverageLocalOfObservations()
 * 2. Standalone computation from a single set of observations:
 *   Call: computeLocal() or computeAverageLocal()
 * 
 * @author Joseph Lizier
 * joseph.lizier at gmail.com
 * http://lizier.me/joseph/
 *
 */
public abstract class InfoMeasureCalculator {

	protected double average = 0.0;
	protected double max = 0.0;
	protected double min = 0.0;
	protected double std = 0.0;
	protected int observations = 0;
	protected int base = 0; // number of individual states. Need initialised to 0 for changedSizes

	protected double log_base = 0;
	protected double log_2 = Math.log(2.0);
	protected boolean power_of_2_base = false;
	protected int log_2_base = 0;
	
	protected boolean debug = false;
	
	/**
	 * 
	 * @param blocksize
	 * @param base
	 */
	protected InfoMeasureCalculator(int base) {
		
		this.base = base;
		log_base = Math.log(base);

		if (base < 2) {
			throw new RuntimeException("Can't calculate info theoretic measures for base " + base);
		}
		
		// Check if we've got a power of 2
		power_of_2_base = isPowerOf2(base);
		if (power_of_2_base) {
			log_2_base = (int) Math.round(Math.log(base) / Math.log(2));
		}
	}
	
	/**
	 * Initialise calculator, preparing to take observation sets in
	 * Should be called prior to any of the addObservations() methods.
	 * You can reinitialise without needing to create a new object.
	 */
	public void initialise(){
		average = 0.0;
		max = 0.0;
		min = 0.0;
		std = 0.0;
		observations = 0;
	}
	
	public final double getLastAverage() {
		return average;
	}

	/**
	 * Not declaring this final so that separable calculator
	 *  can throw an exception on it since it does not support it
	 * 
	 * @return
	 */
	public double getLastMax() {
		return max;
	}

	/**
	 * Not declaring this final so that separable calculator
	 *  can throw an exception on it since it does not support it
	 * 
	 * @return
	 */
	public double getLastMin() {
		return min;
	}
	
	public final double getLastStd() {
		return std;
	}
	
	public final int getNumObservations() {
		return observations;
	}
	
	public final static boolean isPowerOf2(int num) {
		int bits = 0;
		int shiftedValue = num;
		for (int b = 0; b < Integer.SIZE; b ++) {
			if ((shiftedValue & 0x01) > 0) {
				// LSB is a 1
				bits++;
				if (bits > 1) {
					// Encountered more than 1 bit set to 1.
					// Is not a power of 2
					return false;
				}
			}
			// Shift a new bit down into the LSB
			shiftedValue = shiftedValue >> 1;
		}
		// Post: num has either 1 bit set to 1 or 0 bits set to 1
		if (bits == 0) {
			return false;
		}
		return true;
	}
}
