package infodynamics.measures.discrete;

/**
 * Combines functionality for single agents with functionality
 * required in the context of the past.
 * 
 * @author Joseph Lizier
 *
 */
public abstract class SingleAgentMeasureInContextOfPastCalculator extends
		ContextOfPastMeasureCalculator implements SingleAgentMeasure {

	public SingleAgentMeasureInContextOfPastCalculator(int base, int history) {
		super(base, history);
	}

	public final double[][] computeLocal(int[][] states) {
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	public final double[][][] computeLocal(int[][][] states) {
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	public final double computeAverageLocal(int[][] states) {
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	public final double computeAverageLocal(int[][][] states) {
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	public final double[] computeLocalAtAgent(int[][] states, int col) {
		initialise();
		addObservations(states, col);
		return computeLocalFromPreviousObservations(states, col);
	}

	public final double[] computeLocalAtAgent(int[][][] states, int index1, int index2) {
		initialise();
		addObservations(states, index1, index2);
		return computeLocalFromPreviousObservations(states, index1, index2);
	}

	public final double computeAverageLocalAtAgent(int[][] states, int col) {
		initialise();
		addObservations(states, col);
		return computeAverageLocalOfObservations();
	}

}
