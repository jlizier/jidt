package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * Implements conditional mutual information
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
public class ConditionalMutualInformationCalculator extends InfoMeasureCalculator {

	/**
	 * Store the bases for each variable
	 */
	protected int base1;
	protected int base2;
	protected int condBase;
	
	protected int[][][] firstSecondCondCount = null;	// count for (x,y,Cond) tuples
	protected int[][] firstCondCount = null;			// count for (x,Cond) tuples
	protected int[][] secondCondCount = null; // Count for (y,Cond) tuples
	protected int[] condCount = null; // Count for Cond

	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base1 base of first MI variable
	 * @param base2 base of second MI variable
	 * @param condBase base of conditional variable
	 * 
	 * @return
	 */
	public static ConditionalMutualInformationCalculator newInstance(int base1, int base2, int condBase) {
		return new ConditionalMutualInformationCalculator(base1, base2, condBase);
	}

	public ConditionalMutualInformationCalculator(int base1, int base2, int condBase) {

		// Create super object, just with first base
		super(base1);
		
		// Store the bases
		this.base1 = base1;
		this.base2 = base2;
		this.condBase = condBase;
		
		// Create storage for extra counts of observations
		firstSecondCondCount = new int[base1][base2][condBase];
		firstCondCount = new int[base1][condBase];		
		secondCondCount = new int[base2][condBase];
		condCount = new int[condBase];
	}

	/**
	 * Initialise calculator, preparing to take observation sets in
	 * Should be called prior to any of the addObservations() methods.
	 * You can reinitialise without needing to create a new object.
	 *
	 */
	public void initialise(){
		super.initialise();
		
		MatrixUtils.fill(firstSecondCondCount, 0);
		MatrixUtils.fill(firstCondCount, 0);
		MatrixUtils.fill(secondCondCount, 0);
		MatrixUtils.fill(condCount,0);
	}
		
	/**
 	 * Add observations for the given var1,var2,cond tuples of the multi-agent system
 	 *  to our estimates of the pdfs.
	 *
	 * @param var1 values for the first variable
	 * @param var2 values for the second variable
	 * @param cond values for the conditional variable
	 */
	public void addObservations(int var1[], int var2[], int cond[]) {
		int rows = var1.length;
		// increment the count of observations:
		observations += rows;
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			// Add to the count for this particular transition:
			firstSecondCondCount[var1[r]][var2[r]][cond[r]]++;
			firstCondCount[var1[r]][cond[r]]++;
			secondCondCount[var2[r]][cond[r]]++;
			condCount[cond[r]]++;
		}
	}

	/**
 	 * Add observations for the given var1,var2,cond tuples of the multi-agent system
 	 *  to our estimates of the pdfs.
	 *
	 * @param var1 values for the first variable
	 * @param var2 values for the second variable
	 * @param cond values for the conditional variable
	 */
	public void addObservations(int var1[][], int var2[][], int cond[][]) {
		int rows = var1.length;
		int cols = var1[0].length;
		// increment the count of observations:
		observations += rows * cols;
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				// Add to the count for this particular transition:
				firstSecondCondCount[var1[r][c]][var2[r][c]][cond[r][c]]++;
				firstCondCount[var1[r][c]][cond[r][c]]++;
				secondCondCount[var2[r][c]][cond[r][c]]++;
				condCount[cond[r][c]]++;
			}
		}
	}

	/**
	 * Returns the average local conditional MI from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return
	 */
	public double computeAverageLocalOfObservations() {
		double condMi = 0.0;
		double condMiCont = 0.0;

		max = 0;
		min = 0;
		double meanSqLocals = 0;
		for (int condVal = 0; condVal < condBase; condVal++) {
			// compute p(cond)
			// double p_cond = (double) condCount[condVal] / (double) observations;
			for (int var2Val = 0; var2Val < base2; var2Val++) {
				// compute p(var2,cond)
				// double p_var2_cond = (double) seondCondCount[var2Val][condVal] / (double) observations;
				for (int var1Val = 0; var1Val < base1; var1Val++) {
					// compute p(var1,var2,cond)
					double p_var1_var2_cond = (double) firstSecondCondCount[var1Val][var2Val][condVal] / (double) observations;
					// compute p(var1,cond)
					// double p_var1_cond = (double) firstCondCount[var1Val][condVal] / (double) observations;
					// Compute TE contribution:
					if (firstSecondCondCount[var1Val][var2Val][condVal] != 0) {
						/* Double check: should never happen
						if ((sourcePastCount[sourceVal][pastVal] == 0) ||
							(destPastCount[destVal][pastVal] == 0) ||
							(pastCount[pastVal] == 0)) {
							throw new RuntimeException("one subcount was zero!!");
						}
						*/
						
						double logTerm = ((double) firstSecondCondCount[var1Val][var2Val][condVal] / (double) firstCondCount[var1Val][condVal]) /
						 	((double) secondCondCount[var2Val][condVal] / (double) condCount[condVal]);
						double localValue = Math.log(logTerm) / log_2;
						condMiCont = p_var1_var2_cond * localValue;
						if (localValue > max) {
							max = localValue;
						} else if (localValue < min) {
							min = localValue;
						}
						// Add this contribution to the mean 
						//  of the squared local values
						meanSqLocals += condMiCont * localValue;
					} else {
						condMiCont = 0.0;
					}
					condMi += condMiCont;
				}
			}
		}
		
		average = condMi;
		std = Math.sqrt(meanSqLocals - average * average);
		return condMi;
	}
	
	/**
	 * Dump a debug print of the PDFs of our observations
	 */
	public void debugPrintObservations() {
		System.out.println("Var1\tVar2\tCond\tc(1,2,c)\tc(1,c)\tc(2,c)\tc(c)");
		for (int condVal = 0; condVal < condBase; condVal++) {
			// compute p(cond)
			// double p_cond = (double) condCount[condVal] / (double) observations;
			for (int var2Val = 0; var2Val < base2; var2Val++) {
				// compute p(var2,cond)
				// double p_var2_cond = (double) seondCondCount[var2Val][condVal] / (double) observations;
				for (int var1Val = 0; var1Val < base1; var1Val++) {
					// compute p(var1,var2,cond)
					// double p_var1_var2_cond = (double) firstSecondCondCount[var1Val][var2Val][condVal] / (double) observations;
					// compute p(var1,cond)
					// double p_var1_cond = (double) firstCondCount[var1Val][condVal] / (double) observations;
					// Compute TE contribution:
					System.out.println(var1Val + "\t" + var2Val + "\t" + condVal + "\t" +
							firstSecondCondCount[var1Val][var2Val][condVal] + "\t\t" +
							firstCondCount[var1Val][condVal] + "\t" + 
							secondCondCount[var2Val][condVal] + "\t" +
							condCount[condVal]);
				}
			}
		}
	}

	/**
	 * Computes local conditional MI for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 *  
	 * @param var1 values for the first variable
	 * @param var2 values for the second variable
	 * @param cond values for the conditional variable
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int var1[], int var2[], int cond[]){
		int rows = var1.length;

		double[] localCondMi = new double[rows];
		average = 0;
		max = 0;
		min = 0;

		int var1Val, var2Val, condVal;
		double logTerm;
		for (int r = 0; r < rows; r++) {
			var1Val = var1[r];
			var2Val = var2[r];
			condVal = cond[r];
			// Now compute the local value
			logTerm = ((double) firstSecondCondCount[var1Val][var2Val][condVal] / (double) firstCondCount[var1Val][condVal]) /
		 		((double) secondCondCount[var2Val][condVal] / (double) condCount[condVal]);
			localCondMi[r] = Math.log(logTerm) / log_2;
			average += localCondMi[r];
			if (localCondMi[r] > max) {
				max = localCondMi[r];
			} else if (localCondMi[r] < min) {
				min = localCondMi[r];
			}
		}

		average = average/(double) rows;
		
		return localCondMi;
	}

	/**
	 * Compute the significance of obtaining the given average from the given observations,
	 * assuming that the temporal relationship between variable1 and variable2-conditional
	 * was destroyed, while variable2-conditional relationship was retained.
	 * 
	 * @param numPermutationsToCheck number of new orderings of the variable1 to compare against
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) {
		RandomGenerator rg = new RandomGenerator();
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(observations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * Compute the significance of obtaining the given average from the given observations,
	 * assuming that the temporal relationship between variable1 and variable2-conditional
	 * was detroyed, while variable2-conditional relationship was retained.
	 * 
	 * TODO Need to alter the method signature to allow callers to specify
	 * which variable is shuffled. (Note to self: when doing this, will
	 * need to update machine learning code to the new method signature)
	 * 
	 * @param newOrderings the reorderings for variable1 to use
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) {
		
		double actualCondMI = computeAverageLocalOfObservations();
		
		int numPermutationsToCheck = newOrderings.length;
		
		// Reconstruct the observed values of the variables in some order
		int[] var1Values = new int[observations];
		int[] var2Values = new int[observations];
		int[] condValues = new int[observations];
		int t_s = 0;
		for (int val1 = 0; val1 < base1; val1++) {
			for (int val2 = 0; val2 < base2; val2++) {
				for (int condVal = 0; condVal < condBase; condVal++) {
					int numberOfSamples = firstSecondCondCount[val1][val2][condVal];
					MatrixUtils.fill(var1Values, val1, t_s, numberOfSamples);
					MatrixUtils.fill(var2Values, val2, t_s, numberOfSamples);
					MatrixUtils.fill(condValues, condVal, t_s, numberOfSamples);
					t_s += numberOfSamples;
				}
			}
		}
		// We now have arrays of the values that were observed for each
		//  variable, in a random order (well, actually, in order of
		//  increasing joint value of the observations, but this doesn't
		//  matter because:). We will now extract randomly ordered
		//  time series of var1Values, to bootstrap the distribution
		//  of conditional MI values under the null hypothesis.
				
		ConditionalMutualInformationCalculator condMi2 =
				new ConditionalMutualInformationCalculator(base1, base2, condBase);
		condMi2.initialise();
		// Set up the joint counts which remain the same under reordering
		//  of variable 1:
		condMi2.observations = observations;
		condMi2.secondCondCount = secondCondCount;
		condMi2.condCount = condCount;
		int countWhereMIIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the 1st variable
			int[] newData1 = MatrixUtils.extractSelectedTimePoints(var1Values, newOrderings[p]);
			// Compute the required joint probability distributions:
			MatrixUtils.fill(condMi2.firstCondCount, 0);
			MatrixUtils.fill(condMi2.firstSecondCondCount, 0);
			for (int t = 0; t < observations; t++) {
				condMi2.firstCondCount[newData1[t]][condValues[t]]++;
				condMi2.firstSecondCondCount[newData1[t]][var2Values[t]][condValues[t]]++;
			}
			// And get a cond MI value for this realisation of var1Values:
			double newCondMI = condMi2.computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newCondMI;
			if (newCondMI >= actualCondMI) {
				countWhereMIIsMoreSignificantThanOriginal++;
			}
		}
		
		// And return the significance
		measDistribution.pValue = (double) countWhereMIIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualCondMI;
		return measDistribution;
	}

	/**
	 * Standalone routine to 
	 * compute local conditional MI across given variables.
	 * Return a temporal array of local values.
	 * 
	 * @param var1 values for the first variable
	 * @param var2 values for the second variable
	 * @param cond values for the conditional variable
	 * @return
	 */
	public double[] computeLocal(int var1[], int var2[], int cond[]) {
		
		initialise();
		addObservations(var1, var2, cond);
		return computeLocalFromPreviousObservations(var1, var2, cond);
	}

	/**
	 * Standalone routine to 
	 * compute local conditional MI across given variables.
	 * Returns the average
	 * 
	 * @param var1 values for the first variable
	 * @param var2 values for the second variable
	 * @param cond values for the conditional variable
	 * @return
	 */
	public double computeAverageLocal(int var1[], int var2[], int cond[]) {
		
		initialise();
		addObservations(var1, var2, cond);
		return computeAverageLocalOfObservations();
	}
}
