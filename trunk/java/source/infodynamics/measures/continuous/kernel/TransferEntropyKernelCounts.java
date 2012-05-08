package infodynamics.measures.continuous.kernel;

/**
 * Structure to hold the results of the kernel
 *  estimation for one time point
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com
 *
 */
public class TransferEntropyKernelCounts {
	public int countPast;
	public int countNextPast;
	public int countPastSource;
	public int countNextPastSource;
	
	/**
	 * Create a new TransferEntropyKernelCounts object
	 * 
	 * @param countPast
	 * @param countNextPast
	 * @param countPastSource
	 * @param countNextPastSource
	 */
	public TransferEntropyKernelCounts(int countPast,
			int countNextPast, int countPastSource,
			int countNextPastSource) {
		this.countPast = countPast;
		this.countNextPast = countNextPast;
		this.countPastSource = countPastSource;
		this.countNextPastSource = countNextPastSource;
	}
}

