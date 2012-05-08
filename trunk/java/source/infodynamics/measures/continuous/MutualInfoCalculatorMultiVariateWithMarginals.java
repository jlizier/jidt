package infodynamics.measures.continuous;

public interface MutualInfoCalculatorMultiVariateWithMarginals
	extends MutualInfoCalculatorMultiVariate {

	public double computeAverageJointEntropy();

	public double computeAverageEntropyOfObservation1();

	public double computeAverageEntropyOfObservation2();

	public double computeAverageInfoDistanceOfObservations();

	public double[] computeLocalJointEntropyOfPreviousObservations()
		throws Exception;

	public double[] computeLocalJointEntropyUsingPreviousObservations(double states1[][], double states2[][])
		throws Exception;

	public double[] computeLocalEntropy1OfPreviousObservations();

	public double[] computeLocalEntropy1UsingPreviousObservations(double states[][])
		throws Exception;

	public double[] computeLocalEntropy2OfPreviousObservations();

	public double[] computeLocalEntropy2UsingPreviousObservations(double states[][])
		throws Exception;

	public double[] computeLocalInfoDistanceOfPreviousObservations()
		throws Exception;

	public double[] computeLocalInfoDistanceUsingPreviousObservations(double states1[][], double states2[][])
		throws Exception;

}
