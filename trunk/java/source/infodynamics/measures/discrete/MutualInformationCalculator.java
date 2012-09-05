package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;


/**
 * Usage:
 * 1. Continuous accumulation of observations:
 *   Call: a. initialise()
 *         b. addObservations() several times over
 *         c. computeLocalFromPreviousObservations()
 * 2. Standalone:
 *   Call: localMutualInformation()
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com
 *
 */
public class MutualInformationCalculator extends InfoMeasureCalculator 
	implements ChannelCalculator {

	private int timeDiff = 0;
	private int[][]	jointCount = null; // Count for (i[t-timeDiff], j[t]) tuples
	private int[] iCount = null; // Count for i[t-timeDiff]		
	private int[] jCount = null; // Count for j[t]

	public MutualInformationCalculator(int base, int timeDiff) {
		super(base);
		this.timeDiff = timeDiff;
		jointCount = new int[base][base];
		iCount = new int[base];
		jCount = new int[base];
	}

	/**
	 * Initialise calculator, preparing to take observation sets in
	 *
	 */
	public void initialise(){
		super.initialise();
		MatrixUtils.fill(iCount, 0);
		MatrixUtils.fill(jCount, 0);
		MatrixUtils.fill(jointCount, 0);
	}
	
	/**
	 * Add more observations in to our estimates of the pdfs
	 * Pairs are between the arrays var1 and var2, separated in time by timeDiff (i is first)
	 *
	 */
	public void addObservations(int[] var1, int[] var2) {
		int timeSteps = var1.length;
		// int columns = states[0].length;
		// increment the count of observations:
		observations += (timeSteps - timeDiff); 
		
		// 1. Count the tuples observed
		int iVal, jVal;
		for (int r = timeDiff; r < timeSteps; r++) {
			// Add to the count for this particular pair:
			iVal = var1[r-timeDiff];
			jVal = var2[r];
			jointCount[iVal][jVal]++;
			iCount[iVal]++;
			jCount[jVal]++;
		}
	}

	/**
	 * Add more observations in to our estimates of the pdfs
	 * Pairs are between columns iCol and jCol, separated in time by timeDiff (i is first)
	 *
	 */
	public void addObservations(int states[][], int iCol, int jCol) {
		int rows = states.length;
		// int columns = states[0].length;
		// increment the count of observations:
		observations += (rows - timeDiff); 
		
		// 1. Count the tuples observed
		int iVal, jVal;
		for (int r = timeDiff; r < rows; r++) {
			// Add to the count for this particular pair:
			iVal = states[r-timeDiff][iCol];
			jVal = states[r][jCol];
			jointCount[iVal][jVal]++;
			iCount[iVal]++;
			jCount[jVal]++;
		}
	}
	
	/**
	 * Returns the average local mutual information storage from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return
	 */
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		double miCont = 0.0;

		max = 0;
		min = 0;
		double meanSqLocals = 0;
		if (debug) {
			System.out.println("i\tj\tp_i\tp_j\tp_joint\tlocal");
		}
		for (int i = 0; i < base; i++) {
			// compute p_i
			double probi = (double) iCount[i] / (double) observations;
			for (int j = 0; j < base; j++) {
				// compute p_j
				double probj = (double) jCount[j] / (double) observations;
				// compute p(veci=i, vecj=j)
				double jointProb = (double) jointCount[i][j] / (double) observations;
				// Compute MI contribution:
				if (jointProb * probi * probj > 0.0) {
					double localValue = Math.log(jointProb / (probi * probj)) / log_2;
					miCont = jointProb * localValue;
					if (debug) {
						System.out.printf("%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n",
								i, j, probi, probj, jointProb, localValue);
					}
					if (localValue > max) {
						max = localValue;
					} else if (localValue < min) {
						min = localValue;
					}
					// Add this contribution to the mean 
					//  of the squared local values
					meanSqLocals += miCont * localValue;
				} else {
					miCont = 0.0;
				}
				mi += miCont;
			}
		}
		
		average = mi;
		std = Math.sqrt(meanSqLocals - average * average);

		return mi;
	}

	/**
	 * Compute the significance of obtaining the given average from the given observations
	 * 
	 * @param numPermutationsToCheck number of new orderings of the source values to compare against
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) {
		RandomGenerator rg = new RandomGenerator();
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(observations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * Compute the significance of obtaining the given average from the given observations
	 * 
	 * @param newOrderings the reorderings to use
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) {
		double actualMI = computeAverageLocalOfObservations();
		
		int numPermutationsToCheck = newOrderings.length;
		
		// Reconstruct the values of the first and second variables (not necessarily in order)
		int[] iValues = new int[observations];
		int[] jValues = new int[observations];
		int t_i = 0;
		int t_j = 0;
		for (int iVal = 0; iVal < base; iVal++) {
			int numberOfSamplesI = iCount[iVal];
			MatrixUtils.fill(iValues, iVal, t_i, numberOfSamplesI);
			t_i += numberOfSamplesI;
			int numberOfSamplesJ = jCount[iVal];
			MatrixUtils.fill(jValues, iVal, t_j, numberOfSamplesJ);
			t_j += numberOfSamplesJ;
		}
				
		MutualInformationCalculator mi2 = new MutualInformationCalculator(base, timeDiff);
		mi2.initialise();
		mi2.observations = observations;
		mi2.iCount = iCount;
		mi2.jCount = jCount;
		int countWhereMIIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the i variable
			int[] newDataI = MatrixUtils.extractSelectedTimePoints(iValues, newOrderings[p]);
			// compute the joint probability distribution
			MatrixUtils.fill(mi2.jointCount, 0);
			for (int t = 0; t < observations; t++) {
				mi2.jointCount[newDataI[t]][jValues[t]]++;
			}
			// And get an MI value for this realisation:
			double newMI = mi2.computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newMI;
			if (newMI >= actualMI) {
				countWhereMIIsMoreSignificantThanOriginal++;
			}

		}
		
		// And return the significance
		measDistribution.pValue = (double) countWhereMIIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualMI;
		return measDistribution;
	}

	/**
	 * Computes local active information storage for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 *  
	 * @param states
	 * @return
	 */
	double[] computeLocalFromPreviousObservations(int states[][], int iCol, int jCol){
		int rows = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localMI = new double[rows];
		int iVal, jVal;
		double logTerm = 0.0;
		for (int r = timeDiff; r < rows; r++) {
			iVal = states[r-timeDiff][iCol];
			jVal = states[r][jCol];
			logTerm = ( (double) jointCount[iVal][jVal] ) /
			  		  ( (double) jCount[jVal] *
			  			(double) iCount[iVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localMI[r] = Math.log(logTerm) / log_2;
			average += localMI[r];
			if (localMI[r] > max) {
				max = localMI[r];
			} else if (localMI[r] < min) {
				min = localMI[r];
			}
		}
		average = average/(double) observations;
		
		return localMI;
		
	}
	
	/**
	 * Standalone routine to 
	 * compute local active information storage across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * 
	 * @param states - 2D array of states
	 * @return
	 */
	public double[] localMutualInformation(int states[][], int iCol, int jCol) {
		initialise();
		addObservations(states, iCol, jCol);
		return computeLocalFromPreviousObservations(states, iCol, jCol);
		
	}
	
}
