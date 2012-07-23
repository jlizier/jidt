package infodynamics.measures.continuous;

/**
 * <p>
 * Interface to define computers for transfer entropy between a multi-variate source
 * and destination.
 * </p>
 * 
 * @author Joseph Lizier; joseph.lizier at gmail.com
 *
 */
public interface TransferEntropyCalculatorMultiVariate extends ChannelCalculatorMultiVariate {

	public void initialise(int k, int sourceDimensions, int destDimensions) throws Exception;
	
}
