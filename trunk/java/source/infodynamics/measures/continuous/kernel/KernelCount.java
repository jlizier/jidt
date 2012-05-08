package infodynamics.measures.continuous.kernel;

/**
 * Structure to hold a kernel count, the total count it was taken from,
 *  and an array of booleans indicating which time points got counted
 * 
 * @author Joseph Lizier
 *
 */
public class KernelCount {
	public int count;
	public int totalObservationsForCount;
	public boolean[] isCorrelatedWithArgument;
	
	public KernelCount(int count, int totalObservationsForCount, boolean[] isCorrelatedWithArgument) {
		this.count = count;
		this.totalObservationsForCount = totalObservationsForCount;
		this.isCorrelatedWithArgument = isCorrelatedWithArgument;
	}
	
	public KernelCount(int count, int totalObservationsForCount) {
		this(count, totalObservationsForCount, null);
	}
}