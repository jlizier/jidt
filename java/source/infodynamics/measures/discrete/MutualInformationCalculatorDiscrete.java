/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.measures.discrete;

import infodynamics.utils.AnalyticMeasurementDistribution;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Mutual information (MI) calculator for univariate discrete (int[]) data.</p>
 * 
 * <p>Usage of the class is intended to follow this paradigm:</p>
 * <ol>
 * 		<li>Construct the calculator: {@link #MutualInformationCalculatorDiscrete(int)}
 * 			or {@link #MutualInformationCalculatorDiscrete(int, int)};</li>
 *		<li>Initialise the calculator using {@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			sets of {@link #addObservations(int[], int[])} methods, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average MI: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>local MI values, such as
 * 				{@link #computeLocalFromPreviousObservations(int[], int[])};</li>
 * 				<li>comparison to null distribution, such as
 * 				{@link #computeSignificance()};</li>
 * 				<li>and variants of these.</li>
 * 			</ul>
 * 		</li>
 * 		<li>As an alternative to steps 3 and 4, the user may undertake
 * 			standalone computation from a single set of observations, via
 *  		e.g.: {@link #computeLocal(int[][], int, int)}.</li>
 * 		<li>
 * 		Return to step 2 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991).</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MutualInformationCalculatorDiscrete extends InfoMeasureCalculatorDiscrete 
	implements ChannelCalculatorDiscrete, AnalyticNullDistributionComputer {

	/**
	 * Store the number of symbols for each variable
	 */
	protected int alphabetSize1;
	protected int alphabetSize2;

	protected int timeDiff = 0;
	private EntropyCalculatorDiscrete jointEntropyCalc;
	private EntropyCalculatorDiscrete entropyCalc2;
	private EntropyCalculatorDiscrete entropyCalc1;

	private boolean miComputed = false;
	
	
	/**
	 * Construct a new MI calculator with default bases of 2 and time difference of 0
	 *  between the variables
	 * 
	 * @throws Exception
	 */
	public MutualInformationCalculatorDiscrete() throws Exception {
		this(-1);
	}

	/**
	 * Construct a new MI calculator with default time difference of 0
	 *  between the variables
	 * 
	 * @param base number of symbols for each variable. (same for each)
	 *        E.g. binary variables are in base-2.
	 * @throws Exception
	 */
	public MutualInformationCalculatorDiscrete(int alphabetSize) throws Exception {
		this(alphabetSize, alphabetSize, 0);
	}
	
	// TODO Bring in a MutualInformationCalculatorDiscrete(int alphabetSize1, int alphabetSize2) constructor
	//  but don't do this yet since it will override the previous MutualInformationCalculatorDiscrete(int base, int timeDiff)
	//  constructor present up until v1.4 and that may lead to errors.
	
	/**
	 * Create a new mutual information calculator
	 * 
	 * @param alphabetSize1 number of symbols for first variable.
	 *        E.g. binary variables are in base-2.
	 * @param alphabetSize2 number of symbols for second variable.
	 * @param timeDiff number of time steps across which to compute
	 *   MI for given time series
	 * @throws Exception when timeDiff < 0
	 */
	public MutualInformationCalculatorDiscrete(int alphabetSize1, int alphabetSize2, int timeDiff) throws Exception {
		// Create super object, just with first base
		super(alphabetSize1);

		// For unknown alphabet sizes
		if (alphabetSize1 == -1 || alphabetSize2 == -1) {
			this.entropyCalc1 = new EntropyCalculatorDiscrete();
			this.entropyCalc2 = new EntropyCalculatorDiscrete();
		} 
		// For known alphabet sizes
		else {
			this.entropyCalc1 = new EntropyCalculatorDiscrete(alphabetSize1);
			this.entropyCalc2 = new EntropyCalculatorDiscrete(alphabetSize2);
		}
		this.jointEntropyCalc = new EntropyCalculatorDiscrete();
		
		changeBases(alphabetSize1, alphabetSize2, timeDiff);		
	}

	/**
	 * Common code to be called when bases are changed (does not update arrays though)
	 * 
	 * @param alphabetSize1
	 * @param alphabetSize2
	 * @param timeDiff
	 * @throws Exception 
	 */
	private boolean changeBases(int alphabetSize1, int alphabetSize2, int timeDiff) throws Exception {
		boolean basesChanged = false;
		if ((this.alphabetSize1 != alphabetSize1) || (this.alphabetSize2 != alphabetSize2)) {
			basesChanged = true;
		}
		
		// Store the bases
		this.alphabetSize1 = alphabetSize1;
		this.alphabetSize2 = alphabetSize2;
		
		if (timeDiff < 0) {
			throw new Exception("timeDiff must be >= 0");
		}
		this.timeDiff = timeDiff;
		
		return basesChanged;
	}

	@Override
	public void initialise() {
		try {
			initialise(alphabetSize1, alphabetSize2, timeDiff);
		} catch (Exception e) {
			// The only possible (non runtime) exception here is that the timeDiff was < 0
			// which we've already checked, so we can cast this as a Runtime Exception
			throw new RuntimeException(e);
		}
	}

	/**
	 * Initialise with new bases and time diff
	 * 
	 * @param alphabetSize1
	 * @param alphabetSize2
	 * @param timeDiff
	 */
	public void initialise(int alphabetSize1, int alphabetSize2, int timeDiff) {
		super.initialise(alphabetSize1);
		
		this.entropyCalc1.initialise(alphabetSize1);
		this.entropyCalc2.initialise(alphabetSize2);
		this.jointEntropyCalc.initialise();
		jointEntropyCalc.alphabetSizes = new int[] {alphabetSize1, alphabetSize2};
	}

	/**
	 * {@inheritDoc}
	 * 
	 * Pairs for MI are between the arrays var1 and var2, separated in time by timeDiff
	 * (var1 is first).
	 *
	 */
	@Override
	public void addObservations(int[] var1, int[] var2) {

		// I think this won't work if these aren't the same length...
		// Could make it work with different lengths, unsure if this
		// is ever useful.
		if (var1.length != var2.length) {
			throw new RuntimeException("Array lengths do not match. "
			+ "var1: " + var1.length + ", var2: " + var2.length);
		}

		// adjust var1 observations for the time difference
		int[] temp1 = new int[var1.length - timeDiff];
		for (int i = 0; i < var1.length; i++) {
			if (i < timeDiff) {
				continue;
			}
			temp1[i-timeDiff] = var1[i];
		}

		// adjust var2 osbervations for the time difference
		int[] temp2 = new int[var2.length - timeDiff];
		for (int i = 0; i < var2.length - timeDiff; i++) {
			temp2[i-timeDiff] = var2[i];
		}

		entropyCalc1.addObservations(temp1);
		entropyCalc2.addObservations(temp2);

		// Create the joint observations
		int[][] jointObservations = new int[temp1.length][2];
		for (int i = 0; i < temp1.length; i++) {
			jointObservations[i][0] = temp1[i];
			jointObservations[i][1] = temp2[i];
		}

		// This will throw a runtime error if format is somehow incorrect
		jointEntropyCalc.addObservations(jointObservations);
		observations += (var1.length - timeDiff); 
	}

	public void addObservations(Object[] var1, Object[] var2) {
		// I think this won't work  if these aren't the same length...
		if (var1.length != var2.length) {
			throw new RuntimeException("Array lengths do not match. "
			+ "var1: " + var1.length + ", var2: " + var2.length);
		}

		// adjust var1 observations for the time difference
		Object[] temp1 = new Object[var1.length - timeDiff];
		for (int i = 0; i < var1.length; i++) {
			if (i < timeDiff) {
				continue;
			}
			temp1[i-timeDiff] = var1[i];
		}

		// adjust var2 osbervations for the time difference
		Object[] temp2 = new Object[var2.length - timeDiff];
		for (int i = 0; i < var2.length - timeDiff; i++) {
			temp2[i-timeDiff] = var2[i];
		}

		entropyCalc1.addObservations(temp1);
		entropyCalc2.addObservations(temp2);

		// Create the joint observations
		Object[][] jointObservations = new Object[temp1.length][2];
		for (int i = 0; i < temp1.length; i++) {
			jointObservations[i][0] = temp1[i];
			jointObservations[i][1] = temp2[i];
		}

		// This will throw a runtime error if format is somehow incorrect
		jointEntropyCalc.addObservations(jointObservations);
		observations += (var1.length - timeDiff); 
	}
	
	@Override
	public double computeAverageLocalOfObservations() {
		double entropy1 = entropyCalc1.computeAverageLocalOfObservations();
		double entropy2 = entropyCalc2.computeAverageLocalOfObservations();
		double jointEntropy = jointEntropyCalc.computeAverageLocalOfObservations();

		// computeFromPreviousObservations relies on average, so I've kept that being set.
		average = entropy1 + entropy2 - jointEntropy;
		return average;
	}

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) {
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(observations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	// TODO Implement this method. Old code has been left for reference.
	/*
	 * For Brad's understanding, free to delete this comment.
	 * This is:
	 * 	> calculate MI
	 * 	> shuffle all observations
	 * 	> recalculate MI for each permutation
	 * 	> count number of recalcs > original MI
	 * 	> p = count / numPerms
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) {
		
		// 1/3 temporary fix to get it to compile
		// all original code is commented out below
		return null;
		/*
		double actualMI = computeAverageLocalOfObservations();
		
		int numPermutationsToCheck = newOrderings.length;
		
		// Reconstruct the values of the first and second variables (not necessarily in order)
		int[] iValues = new int[observations];
		int t_i = 0;
		for (int iVal = 0; iVal < alphabetSize1; iVal++) {
			int numberOfSamplesI = iCount[iVal];
			MatrixUtils.fill(iValues, iVal, t_i, numberOfSamplesI);
			t_i += numberOfSamplesI;
		}
		int[] jValues = new int[observations];
		int t_j = 0;
		for (int jVal = 0; jVal < alphabetSize2; jVal++) {
			int numberOfSamplesJ = jCount[jVal];
			MatrixUtils.fill(jValues, jVal, t_j, numberOfSamplesJ);
			t_j += numberOfSamplesJ;
		}
		
		MutualInformationCalculatorDiscrete mi2;
		try {
			mi2 = new MutualInformationCalculatorDiscrete(alphabetSize1, alphabetSize2, timeDiff);
		} catch (Exception e) {
			// The only possible exception is if timeDiff < 0, which 
			// it cannot be. Shut down the JVM
			throw new Error("timeDiff parameter took on value < 0 after being checked at construction");
		}
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
		*/
	}

	@Override
	public AnalyticMeasurementDistribution computeSignificance() {
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}

		int a1 = alphabetSize1;
		int a2 = alphabetSize2;
		if (alphabetSize1 <= -1 || alphabetSize2 == -1) {
			a1 = entropyCalc1.hashedStateCount.keySet().size();
			a2 = entropyCalc2.hashedStateCount.keySet().size();
		}

		return new ChiSquareMeasurementDistribution(average,
				observations,
				(a1 - 1) * (a2 - 1));
	}
	
	/**
	 * Computes local mutual information (or pointwise mutual information)
	 *  for the given (single) specific values, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 *  
	 * @param val1 single specific value of variable 1
	 * @param val2 single specific value of variable 2
	 * @return a local mutual information value for this pair of observations
	 * @throws Exception if this pair were not observed together in the
	 *  previously supplied observations
	 */
	public double computeLocalFromPreviousObservations(int val1, int val2) throws Exception{
		
		// 2/3 temporary fix to get this to compile
		return 0;
		/*
		double logTerm = ( (double) jointCount[val1][val2] ) /
			  		  ( (double) jCount[val2] *
			  			(double) iCount[val1] );
		// Now account for the fact that we've
		//  just used counts rather than probabilities,
		//  and we've got two counts on the bottom
		//  but one count on the top:
		logTerm *= (double) observations;
		double localMI = Math.log(logTerm) / log_2;
		
		return localMI;
		*/
	}
	
	/**
	 * Computes local mutual information (or pointwise mutual information)
	 *  for the given states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 *  
	 * @param var1 new states of variable 1
	 * @param var2 new states of variable 2. Should be same length as var1
	 * @return array of local mutual information values for each
	 *  observation of (var1, var2). Note - if timeDiff > 0, then the
	 *  return length will be var1.length - timeDiff (in line with continuous estimators). 
	 */
	public double[] computeLocalFromPreviousObservations(int[] var1, int[] var2) throws Exception{
		
		// 3/3 temporary fix to get it to compile
		return null;
		/*
		if (var1.length != var2.length) {
			throw new Exception("var1 and var2 must have the same number of observations");
		}
		double[] localMI = new double[var1.length - timeDiff];
		
		double logTerm = 0.0;
		for (int r = timeDiff; r < var1.length; r++) {
			int iVal = var1[r-timeDiff];
			int jVal = var2[r];
			logTerm = ( (double) jointCount[iVal][jVal] ) /
			  		  ( (double) jCount[jVal] *
			  			(double) iCount[iVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localMI[r - timeDiff] = Math.log(logTerm) / log_2;
			average += localMI[r - timeDiff];
			if (localMI[r - timeDiff] > max) {
				max = localMI[r - timeDiff];
			} else if (localMI[r - timeDiff] < min) {
				min = localMI[r - timeDiff];
			}
		}
		average = average/(double) observations;

		return localMI;
		*/
	}

	// TODO -- do not implement, just remove this from the interface
	@Override
	public void addObservations(int[][] states, int sourceIndex, int destIndex) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Unimplemented method 'addObservations'");
	}
}
