package infodynamics.measures.continuous;

public interface TransferEntropyCalculator extends ChannelCalculator {

	public static final String K_PROP_NAME = "k_HISTORY";
	
	public void initialise(int k) throws Exception;
	
}
