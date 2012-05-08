package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.TransferEntropyCalculator;
import infodynamics.measures.continuous.TransferEntropyCommon;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.MeasurementDistribution;
import infodynamics.utils.RandomGenerator;

import java.util.Hashtable;
import java.util.Iterator;


/**
 * <p>Compute the Transfer Entropy using the Kraskov estimation method.<p>
 * <p>Transfer Entropy is defined in:
 * <ul>
 * 		<li>"Measuring information transfer", T. Schreiber, Physical Review E 85, 2,
 * 		(2000), pp. 461-464</li>
 * 		 <li>"Information transfer in continuous processes", A. Kaiser, T. Schreiber,
 * 		Physica D, 166 (2002), pp. 43-62</li>
 * </ul>
 * Interestingly, in the Physica D paper, the authors begin discussing in section 5.2.3.2 how
 *  to use "fixed mass" kernel estimators across all of the underlying entropies, but
 *  state they cannot work out how to correct for the correlations between the counts.
 * </p>
 * 
 * <p>This correction is effectively what was done in Kraskov et al's method for computing the
 *  Mutual information:
 * <ul>
 * 		<li>"Estimating mutual information", Kraskov, A., Stogbauer, H., Grassberger, P.,
 * 		Physical Review E 69, (2004) 066138</li>
 * 		<li>"Synchronization and Interdependence Measures and their Applications to the
 * 		Electroencephalogram of Epilepsy Patients and Clustering of Data",
 *		Alexander Kraskov, Publication Series of the John von Neumann Institute for Computing,
 *		Volume 24, John von Neumann Institute for Computing, J\"{u}lich, Germany, 2004 (PhD Thesis)
 *		</li>
 * </ul>
 * </p>
 * 
 * <p>Now, in Kraskov "Synchronization ..." the author presents a method for estimating
 * Transfer Entropy as a sum of two mutual informations (each estimated by his technique).
 * We implement this estimation of the transfer entropy here.
 * We compute <code>I[ (X',X^k); Y ] - I [X^k;Y ]</code> rather than
 * <code>I[ (Y,X^k); X' ] - I [X^k;X' ]</code> as 
 *  Kraskov found this to be the more accurate option.
 * </p>
 * 
 * <p>
 * Usage:
 * <ol>
 * 		<li>Construct</li>
 * 		<li>SetProperty() for each property</li>
 * 		<li>intialise()</li>
 * 		<li>setObservations(), or [startAddObservations(), addObservations()*, finaliseAddObservations()]
 *   Note: If not using setObservations(), the results from computeLocal or getSignificance
 *    are not likely to be particularly sensible.</li> 
 * 		<li>computeAverageLocalOfObservations() or ComputeLocalOfPreviousObservations()</li>
 * </ol>
 * </p>
 * 
 * @author Joseph Lizier joseph.lizier at gmail.com
 *
 */
public class TransferEntropyCalculatorKraskovByMulti
	extends TransferEntropyCommon implements TransferEntropyCalculator {

	/**
	 * Which Kraskov MI algorithm number to use (1 or 2)
	 */
	protected int kraskovAlgorithmNumber = 2;
	private boolean algChanged = false;
	/**
	 * Storage for the properties ready to pass onto the underlying MI calculators
	 */
	private Hashtable<String,String> props;
	/**
	 * MI calculator for (Past,Next) to Source
	 */
	protected MutualInfoCalculatorMultiVariate mickNextPastToSource;
	/**
	 * MI calculator for Past to Source
	 */
	protected MutualInfoCalculatorMultiVariate mickPastToSource;

	public final static String PROP_KRASKOV_ALG_NUM = "ALG_NUM";

	protected double[][] jointPastVectors;
	protected double[][] jointNextAndPastVectors;
	protected double[][] sourceVectors;

	/**
	 * Construct the calculator
	 *
	 */
	public TransferEntropyCalculatorKraskovByMulti() {
		super();
		props = new Hashtable<String,String>();
		createKraskovMiCalculators();
	}

	public void initialise(int k) throws Exception {
		super.initialise(k); // calls initialise();
	}

	/**
	 * Initialise using default or existing values for k and epsilon
	 */
	public void initialise() throws Exception {
		if (algChanged) {
			// Create the Kraskov MI calculators for the new algorithm number.
			createKraskovMiCalculators();
			algChanged = false;
		}
		// Set the properties for the Kraskov MI calculators
		for (String key : props.keySet()) {
			mickPastToSource.setProperty(key, props.get(key));
			mickNextPastToSource.setProperty(key, props.get(key));
		}
		// Initilalise the Kraskov MI calculators
		initialiseKraskovCalculators();
		addedMoreThanOneObservationSet = false;
	}
	
	/**
	 * Create the underlying Kraskov MI calculators 
	 *
	 */
	protected void createKraskovMiCalculators() {
		if (kraskovAlgorithmNumber == 1) {
			mickPastToSource = new MutualInfoCalculatorMultiVariateKraskovByMulti1();
			mickNextPastToSource = new MutualInfoCalculatorMultiVariateKraskovByMulti1();
		} else {
			// Algorithm 2
			mickPastToSource = new MutualInfoCalculatorMultiVariateKraskovByMulti2();
			mickNextPastToSource = new MutualInfoCalculatorMultiVariateKraskovByMulti2();
		}
	}
	
	protected void initialiseKraskovCalculators() throws Exception {
		mickPastToSource.initialise(k, 1);
		mickNextPastToSource.initialise(k+1, 1);
	}
	
	/**
	 * Sets properties for the calculator.
	 * Valid properties include:
	 * <ul>
	 *  <li>K_PROP_NAME</li>
	 * 	<li>PROP_KRASKOV_ALG_NUM</li>
	 * 	<li>Any valid properties for MutualInfoCalculatorMultiVariateKraskov.setProperty except
	 * 			for MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF</li>
	 * </ul>
	 * One should set MutualInfoCalculatorMultiVariateKraskov.PROP_K here, the number
	 *  of neighbouring points one should count up to in determining the joint kernel size. 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		super.setProperty(propertyName, propertyValue);
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_KRASKOV_ALG_NUM)) {
			int previousAlgNumber = kraskovAlgorithmNumber;
			kraskovAlgorithmNumber = Integer.parseInt(propertyValue);
			if ((kraskovAlgorithmNumber != 1) && (kraskovAlgorithmNumber != 2)) {
				throw new Exception("Kraskov algorithm number (" + kraskovAlgorithmNumber
						+ ") must be either 1 or 2");
			}
			if (kraskovAlgorithmNumber != previousAlgNumber) {
				algChanged = true;
			}
		} else if (propertyName.equalsIgnoreCase(MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF)) {
			// Do nothing, simply prevent this property being set on the 
			//  mutual info calculator.
		} else {
			// Assume it was a property for the MI calculator
			props.put(propertyName, propertyValue);
		}
	}

	/**
	 * Flag that the observations are complete, probability distribution functions can now be built.
	 *
	 */
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[] destination : vectorOfDestinationObservations) {
			totalObservations += destination.length - k;
		}
		jointPastVectors = new double[totalObservations][k];
		jointNextAndPastVectors = new double[totalObservations][k+1];
		sourceVectors = new double[totalObservations][1];

		// Construct the joint vectors from the given observations
		int startObservation = 0;
		Iterator<double[]> iterator = vectorOfDestinationObservations.iterator();
		for (double[] source : vectorOfSourceObservations) {
			double[] destination = iterator.next();
			double[][] currentDestPastVectors = makeJointVectorForPast(destination);
			MatrixUtils.arrayCopy(currentDestPastVectors, 0, 0,
					jointPastVectors, startObservation, 0, currentDestPastVectors.length, k);
			double[][] currentDestNextPastVectors = makeJointVectorForNextPast(destination);
			MatrixUtils.arrayCopy(currentDestNextPastVectors, 0, 0,
					jointNextAndPastVectors, startObservation, 0, currentDestNextPastVectors.length, k + 1);
			try {
				MatrixUtils.copyIntoColumn(sourceVectors, 0, startObservation,
					source, k - 1, source.length - k);
			} catch (Exception e) {
				// The above params should not throw an Exception, so
				//  wrap in a runtime exception
				throw new RuntimeException(e);
			}
			startObservation += destination.length - k;
		}
		
		// Now set the joint vectors in the kernel estimators
		try {
			mickPastToSource.setObservations(jointPastVectors, sourceVectors);
			mickNextPastToSource.setObservations(jointNextAndPastVectors, sourceVectors);
		} catch (Exception e) {
			// The above should not throw an exception since they were constructed here
			//  of the same time length, so wrap in a runtime exception
			throw new RuntimeException(e);
		}

		// Store whether there was more than one observation set:
		addedMoreThanOneObservationSet = vectorOfDestinationObservations.size() > 1;
		
		// And clear the vector of observations
		vectorOfSourceObservations = null;
		vectorOfDestinationObservations = null;
	}

	/**
	 * Returns the transfer entropy calculated as I[ (X',X^k); Y ] - I [X^k;Y ]
	 * using the underlying Kraskov MI calculators 
	 * 
	 * @return the transfer entropy
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		double miPastNextToSource = mickNextPastToSource.computeAverageLocalOfObservations();
		shareDataBetweenUnderlyingCalculators();
		double miPastToSource = mickPastToSource.computeAverageLocalOfObservations();
		lastAverage = miPastNextToSource - miPastToSource;
		if (debug) {
			System.out.println("miPastNextToSource=" + miPastNextToSource + ", miPastToSource=" +
					miPastToSource + " = " + lastAverage);
		}
		return lastAverage;
	}

	/**
	 * If mickPastToSource can utilise anything from mickPastNextToSource after the latter
	 * has run computeAverageLocalOfObservations, arrange that here 
	 *
	 */
	protected void shareDataBetweenUnderlyingCalculators() {
		if (! MutualInfoCalculatorMultiVariateKraskovByMulti.class.isInstance(mickNextPastToSource)) {
			// We don't know of what to share for other calculator types.
			// Subclasses may know and can over-ride this method.
			return;
		}
		MutualInfoCalculatorMultiVariateKraskovByMulti micmvkNextPastToSource =
			(MutualInfoCalculatorMultiVariateKraskovByMulti) mickNextPastToSource;
		MutualInfoCalculatorMultiVariateKraskovByMulti micmvkPastToSource =
			(MutualInfoCalculatorMultiVariateKraskovByMulti) mickPastToSource;
		if (micmvkNextPastToSource.multiInfoJoint.norms != null) {
			// Share the norms already computed with the other mutual info calculator.
			// Just assign the joint norms for it, it will filter them through to the
			//  marginals itself.
			if (micmvkPastToSource.multiInfoJoint.norms == null) {
				micmvkPastToSource.multiInfoJoint.norms = new double[k + 1][][];
				// Point to the norms for the past variables
				for (int v = 0; v < k; v++) {
					micmvkPastToSource.multiInfoJoint.norms[v] =
						micmvkNextPastToSource.multiInfoJoint.norms[1+v];
				}
				// Point to the norms for the source variable
				micmvkPastToSource.multiInfoJoint.norms[k] =
					micmvkNextPastToSource.multiInfoJoint.norms[k + 1];
			}
		}
	}
	
	/**
	 * Returns the local transfer entropy values.
	 * 
	 * Where more than one time series has been added, the array
	 *  contains the local values for each tuple in the order in
	 *  which they were added.
	 * 
	 * If there was only a single time series added, the array
	 *  contains k zero values before the local values.
	 *  (This means the length of the return array is the same
	 *  as the length of the input time series).
	 * 
	 * @return an array of the local transfer entropy values.
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		// Initialise localTE with local MI of PastNextToSource
		double[] localMiPastNextToSource = mickNextPastToSource.computeLocalOfPreviousObservations();
		shareDataBetweenUnderlyingCalculators();
		double[] localMiPastToSource = mickPastToSource.computeLocalOfPreviousObservations();
		MatrixUtils.subtractInPlace(localMiPastNextToSource, localMiPastToSource);
		lastAverage = MatrixUtils.mean(localMiPastNextToSource);
		
		if (debug) {
			if (!addedMoreThanOneObservationSet) {
				for (int i = 0; i < k; k++) {
					System.out.println(0.0);
				}
			}
			for (int t = 0; t < localMiPastNextToSource.length; t++) {
				System.out.println(sourceVectors[t][0] + " -> " + 
					jointNextAndPastVectors[t][0] + ", " +
					jointNextAndPastVectors[t][1] + ": " +
					(localMiPastNextToSource[t]+localMiPastToSource[t]) +
					" - " + localMiPastToSource[t] +
				    " = " + localMiPastNextToSource[t]);
			}
		}
		
		if (!addedMoreThanOneObservationSet) {
			// We need to allow k zeros at the beginning
			int offset = k;
			double[] localTE = new double[localMiPastNextToSource.length + offset];
			System.arraycopy(localMiPastNextToSource, 0, localTE, offset, localMiPastNextToSource.length);
			return localTE;
		} else {
			return localMiPastNextToSource;
		}
	}

	/**
	 * Compute the significance of obtaining the given average TE from the given observations
	 * 
	 * 	This is as per Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128.
	 *
	 * Basically, we shuffle the source observations against the destination tuples.
	 * This keeps the marginal PDFs the same (including the entropy rate of the destination)
	 *  but destroys any correlation between the source and state change of the destination.
	 * 
	 * @param numPermutationsToCheck number of new orderings of the source values to compare against
	 * @return
	 */
	public MeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(totalObservations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * As per {@link computeSignificance(int) computeSignificance()} but supplies
	 *  the re-orderings of the observations of the source variables.
	 * 
	 * Implemented using "virtual" reorderings of the source variable (assuming there
	 *  are a small enough number of time points).
	 * 
	 * @param newOrderings first index is permutation number, i.e. newOrderings[i]
	 * 		is an array of 1 permutation of 0..n-1, where there were n observations.
	 * 		If the length of each permutation in newOrderings
	 * 		is not equal to numObservations, an Exception is thrown. 
	 * @return
	 * @throws Exception
	 */
	public MeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		
		double actualTE = computeAverageLocalOfObservations();
		
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		MeasurementDistribution measDistribution = new MeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the source in the destPastSourceVectors 
			//  and destNextPastSourceVectors vectors
			if (newOrderings[p].length != totalObservations) {
				throw new Exception("permutation " + p +
						" in newOrderings argument was not equal to totalObservations");
			}
			
			// Compute the new component MIs under the "virtual" reordering of the source
			//  variable (it is the second variable for each of these MI calculators, so it
			//  is the one which will be reordered.
			double newMiPastNextToSource = mickNextPastToSource.computeAverageLocalOfObservations(newOrderings[p]);
			double newMiPastToSource = mickPastToSource.computeAverageLocalOfObservations(newOrderings[p]);
			double newTe = newMiPastNextToSource - newMiPastToSource;
			
			measDistribution.distribution[p] = newTe;
			if (newTe >= actualTE) {
				countWhereTeIsMoreSignificantThanOriginal++;
			}
		}
				
		// Return the significance
		measDistribution.pValue = (double) countWhereTeIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualTE;
		return measDistribution;
	}

	/**
	 * As per {@link computeSignificance(int) computeSignificance()} but supplies
	 *  the re-orderings of the observations of the source variables.
	 * Explicitly reorders the source variable instead of allowing the underling MI
	 *  calculators to do it virtually.
	 *  
	 * 
	 * @param newOrderings first index is permutation number, i.e. newOrderings[i]
	 * 		is an array of 1 permutation of 0..n-1, where there were n observations.
	 * 		If the length of each permutation in newOrderings
	 * 		is not equal to numObservations, an Exception is thrown. 
	 * @return
	 * @throws Exception
	 */
	public MeasurementDistribution computeSignificanceExplicitlyReordering(
			int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		
		double actualTE = computeAverageLocalOfObservations();
		
		// Save the relevant source observations here:
		double[][] originalSourceValues = new double[totalObservations][];
		for (int t = 0; t < totalObservations; t++) {
			originalSourceValues[t] = sourceVectors[t];
		}
		
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		MeasurementDistribution measDistribution = new MeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the source in the destPastSourceVectors 
			//  and destNextPastSourceVectors vectors
			if (newOrderings[p].length != totalObservations) {
				throw new Exception("permutation " + p +
						" in newOrderings argument was not equal to totalObservations");
			}
			sourceVectors = MatrixUtils.extractSelectedTimePointsReusingArrays(
					originalSourceValues, newOrderings[p]);
			
			// Make the equivalent operations of intialise
			mickPastToSource.initialise(k, 1);
			mickNextPastToSource.initialise(k+1, 1);
			// Make the equivalent operations of setObservations:
			try {
				mickPastToSource.setObservations(jointPastVectors, sourceVectors);
				mickNextPastToSource.setObservations(jointNextAndPastVectors, sourceVectors);
			} catch (Exception e) {
				// The above should not throw an exception since they were constructed here
				//  of the same time length, so wrap in a runtime exception
				throw new RuntimeException(e);
			}
			// And get a TE value for this realisation:
			double newTe = computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newTe;
			if (newTe >= actualTE) {
				countWhereTeIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the local variables:
		lastAverage = actualTE;
		// Restore the source observations in the joint vectors
		for (int t = 0; t < sourceVectors.length; t++) {
			sourceVectors[t] = originalSourceValues[t];
		}
		// And set the MI estimators back to their previous state
		mickPastToSource.initialise(k, 1);
		mickNextPastToSource.initialise(k+1, 1);
		try {
			mickPastToSource.setObservations(jointPastVectors, sourceVectors);
			mickNextPastToSource.setObservations(jointNextAndPastVectors, sourceVectors);
		} catch (Exception e) {
			// The above should not throw an exception since they were constructed here
			//  of the same time length, so wrap in a runtime exception
			throw new RuntimeException(e);
		}
		
		// And return the significance
		measDistribution.pValue = (double) countWhereTeIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = lastAverage;
		return measDistribution;
	}

	public void setDebug(boolean debug) {
		super.setDebug(debug);
		mickNextPastToSource.setDebug(debug);
		mickPastToSource.setDebug(debug);
	}
}
