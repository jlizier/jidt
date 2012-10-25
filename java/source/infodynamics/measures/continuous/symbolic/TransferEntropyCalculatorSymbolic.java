package infodynamics.measures.continuous.symbolic;

import infodynamics.measures.continuous.TransferEntropyCalculator;
import infodynamics.measures.continuous.TransferEntropyCommon;
import infodynamics.measures.discrete.ApparentTransferEntropyCalculator;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

import java.util.Iterator;
import java.util.Vector;


/**
 * Computing the transfer entropy symbollically
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com
 * @see Mattha\"{u}s Staniek, Klaus Lehnertz, "Symbolic Transfer Entropy", PRL 100, 158101 (2008)
 *
 */
public class TransferEntropyCalculatorSymbolic
		extends TransferEntropyCommon
		implements TransferEntropyCalculator {

	// The calculator used to do the grunt work
	protected ApparentTransferEntropyCalculator teCalc;

	protected int maxEmbeddingLength; // = Max(l,k)
	// l will default to the value of k unless l gets explicitly set
	protected int l = 1;
	protected boolean lWasSet = false;
	
	private int[][] destPermutations;
	// For each permutation index, holds the unique permutation id
	private int[] destPermutationIds;
	// For each possible permutation id, holds the permutation index 
	private int[] idToDestPermutationIndex;
	
	private int[][] sourcePermutations;
	// For each permutation index, holds the unique permutation id
	private int[] sourcePermutationIds;
	// For each possible permutation id, holds the permutation index 
	private int[] idToSourcePermutationIndex;

	// Storage for the computation symbols
	protected Vector<int[]> destSymbolsVector;
	protected Vector<int[]> sourceSymbolsVector;
	
	// Array indices for the 2D array sorted by the first index
	protected static final int VAL_COLUMN = 0;
	protected static final int VAR_NUM_COLUMN = 1;

	public static final String PROP_L = "L";
	
	public TransferEntropyCalculatorSymbolic() {
		// Nothing to do
	}

	public void initialise(int k) throws Exception {
		super.initialise(k); // calls initialise();
	}

	public void initialise() throws Exception {
		if (!lWasSet) {
			l = k;
		}
		
		// The discrete TE calculator will be run with a base of the number of
		//  permutations of the maximum embedding length (out of l,k), and history 1
		//  since the symbols incorporate the past k states.
		maxEmbeddingLength = Math.max(k, l);
		
		// First work out how many permutations of orderings we could have
		RandomGenerator rg = new RandomGenerator();
		destPermutations = rg.generateAllDistinctPerturbations(k);
		// Now generate an int signature for each permutation:
		destPermutationIds = new int[destPermutations.length];
		for (int r = 0; r < destPermutations.length; r++) {
			destPermutationIds[r] = generatePermutationId(destPermutations[r], k);
		}
		// Now we have a list of permutations, each with an ID (which will be up to k^k)
		// Generate a reverse mapping from permutation identifier to permutation id
		idToDestPermutationIndex = new int[MathsUtils.power(k, k)];
		// First initialise all mappings to -1 : this will force an Array lookup error if 
		//  an identifier is not calculated correctly (most of the time)
		for (int i = 0; i < idToDestPermutationIndex.length; i++) {
			idToDestPermutationIndex[i] = -1;
		}
		for (int idIndex = 0; idIndex < destPermutationIds.length; idIndex++) {
			idToDestPermutationIndex[destPermutationIds[idIndex]] = idIndex;
			/* System.out.print("Permutation "+ idIndex + ": ");
			for (int j = 0; j < permutations[idIndex].length; j++) {
				System.out.print(permutations[idIndex][j] + " ");
			}
			System.out.println(" -> id " + permutationIds[idIndex]);
			*/
		}

		// Now generate the permutations for the source
		sourcePermutations = rg.generateAllDistinctPerturbations(l);
		// Now generate an int signature for each permutation:
		sourcePermutationIds = new int[sourcePermutations.length];
		for (int r = 0; r < sourcePermutations.length; r++) {
			sourcePermutationIds[r] = generatePermutationId(sourcePermutations[r], l);
		}
		// Now we have a list of permutations, each with an ID (which will be up to l^l)
		// Generate a reverse mapping from permutation identifier to permutation id
		idToSourcePermutationIndex = new int[MathsUtils.power(l, l)];
		// First initialise all mappings to -1 : this will force an Array lookup error if 
		//  an identifier is not calculated correctly (most of the time)
		for (int i = 0; i < idToSourcePermutationIndex.length; i++) {
			idToSourcePermutationIndex[i] = -1;
		}
		for (int idIndex = 0; idIndex < sourcePermutationIds.length; idIndex++) {
			idToSourcePermutationIndex[sourcePermutationIds[idIndex]] = idIndex;
			/* System.out.print("Permutation "+ idIndex + ": ");
			for (int j = 0; j < permutations[idIndex].length; j++) {
				System.out.print(permutations[idIndex][j] + " ");
			}
			System.out.println(" -> id " + permutationIds[idIndex]);
			*/
		}

		// The discrete calculator only uses a history of 1 here - k is built into
		//  the permutations
		int base = Math.max(destPermutationIds.length, sourcePermutationIds.length);
		teCalc = ApparentTransferEntropyCalculator.newInstance(base, 1);
		teCalc.initialise();
	}

	/**
	 * <p>Set properties for the transfer entropy calculator.
	 * These can include:
	 * <ul>
	 * 		<li>PROP_L</li>
	 * </ul> 
	 * and those for {@link TransferEntropyCommon#setProperty(String,String)}
	 * </p>
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception 
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		super.setProperty(propertyName, propertyValue);
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_L)) {
			l = Integer.parseInt(propertyValue);
			lWasSet = true;
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println("Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[] destination : vectorOfDestinationObservations) {
			totalObservations += destination.length - maxEmbeddingLength;
		}

		destSymbolsVector = new Vector<int[]>();
		sourceSymbolsVector = new Vector<int[]>();
		
		// Construct the symbols from the given observations
		//int startObservation = 0;
		Iterator<double[]> iterator = vectorOfDestinationObservations.iterator();
		for (double[] source : vectorOfSourceObservations) {
			double[] destination = iterator.next();
			// compute the embedding vectors for the destination and source
			double[][] currentDestPastVectors = null;
			double[][] sourceVectors = null;
			try {
				// we don't use the pre-defined makeDestPastVectors method because it skips off 
				//  the last embedding vector
				currentDestPastVectors = MatrixUtils.makeDelayEmbeddingVector(destination, k, k-1, destination.length - k + 1);
				sourceVectors = MatrixUtils.makeDelayEmbeddingVector(source, l, l-1, source.length - l + 1);
			} catch (Exception e) {
				// The parameters for the above call should be fine, so we don't expect to
				//  throw an Exception here - embed in a RuntimeException if it occurs 
				throw new RuntimeException(e);
			}
			
			// Now compute the permutation values for the dest
			double[][] destVariablesAndIndices = new double[k][2];
			int[] destSymbols = new int[currentDestPastVectors.length];
			for (int t = 0; t < currentDestPastVectors.length; t++) {
				// Work out what the order of the embedded variables was here:
				for (int v = 0; v < k; v++) {
					destVariablesAndIndices[v][VAL_COLUMN] = currentDestPastVectors[t][v];
					destVariablesAndIndices[v][VAR_NUM_COLUMN] = v;
				}
				java.util.Arrays.sort(destVariablesAndIndices, FirstIndexComparatorDouble.getInstance());
				// Now the second column contains the order of values here
				double[] permutation = MatrixUtils.selectColumn(destVariablesAndIndices, VAR_NUM_COLUMN);
				int permutationId = generatePermutationId(permutation, k);
				destSymbols[t] = idToDestPermutationIndex[permutationId];
			}

			// Compute the permutation values for the source
			double[][] sourceVariablesAndIndices = new double[l][2];
			int[] sourceSymbols = new int[sourceVectors.length];
			for (int t = 0; t < sourceVectors.length; t++) {
				// Work out what the order of the embedded variables was here:
				for (int v = 0; v < l; v++) {
					sourceVariablesAndIndices[v][VAL_COLUMN] = sourceVectors[t][v];
					sourceVariablesAndIndices[v][VAR_NUM_COLUMN] = v;
				}
				java.util.Arrays.sort(sourceVariablesAndIndices, FirstIndexComparatorDouble.getInstance());
				// Now the second column contains the order of values here
				double[] permutation = MatrixUtils.selectColumn(sourceVariablesAndIndices, VAR_NUM_COLUMN);
				int permutationId = generatePermutationId(permutation, l);
				sourceSymbols[t] = idToSourcePermutationIndex[permutationId];
			}
			
			if (k > l) {
				// l is smaller, so we will have more source embeddings than destination ones.
				// Get rid of the first k - l source embeddings
				sourceSymbols = MatrixUtils.select(sourceSymbols, k - l,
						sourceSymbols.length - (k - l));
			} else if (l > k) {
				// k is smaller, so we will have more dest embeddings than source ones.
				// Get rid of the first l - k source embeddings
				destSymbols = MatrixUtils.select(destSymbols, l - k,
						destSymbols.length - (l - k));
			}
			
			// Add these observations in, and keep them locally
			destSymbolsVector.add(destSymbols);
			sourceSymbolsVector.add(sourceSymbols);
			teCalc.addObservations(destSymbols, sourceSymbols);
			
			if (destination.length - maxEmbeddingLength != destSymbols.length - 1) {
				throw new RuntimeException(
					String.format("Number of observations %d doesn't match what's expected %d",
						destSymbols.length - 1, destination.length - maxEmbeddingLength));
			}
			
			// I don't remember why we were tracking this:
			// (probably before using vectors of observations)
			//startObservation += destSymbols.length - 1;
		}
		
	}

	/**
	 * Generate the unique permutation id for this permutation.
	 * 
	 * @param data
	 * @param base the number of elements being permuted
	 * @return
	 */
	private int generatePermutationId(int[] data, int base) {
		int permutationId = 0;
		for (int c = 0; c < data.length; c++) {
			permutationId *= base;
			permutationId +=  data[c];
		}
		return permutationId;
	}

	/**
	 * Generate the unique permutation id for this ordering of variable ids.
	 * Convert the floating point variable ids into ints first
	 * 
	 * @param data
	 * @param base the number of elements being permuted
	 * @return
	 */
	private int generatePermutationId(double[] ids, int base) {
		int permutationId = 0;
		for (int c = 0; c < ids.length; c++) {
			permutationId *= base;
			permutationId +=  (int) ids[c];
		}
		return permutationId;
	}

	public double computeAverageLocalOfObservations() throws Exception {
		lastAverage = teCalc.computeAverageLocalOfObservations();
		return lastAverage;
	}

	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] locals = new double[totalObservations];
		int currentIndexInLocals = 0;
		Iterator<int[]> iterator = destSymbolsVector.iterator();
		for (int[] sourceSymbols : sourceSymbolsVector) {
			int[] destSymbols = iterator.next();
			double[] theseLocals = teCalc.computeLocalFromPreviousObservations(destSymbols, sourceSymbols);
			System.arraycopy(theseLocals, 1, locals, currentIndexInLocals, theseLocals.length - 1);
			currentIndexInLocals += theseLocals.length - 1;
		}
		lastAverage = MatrixUtils.mean(locals);
		return locals;
	}

	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		return teCalc.computeSignificance(numPermutationsToCheck);
	}

	/**
	 * This method does not actually use the newOrderings supplied,
	 *  which is not strictly what this method is meant to do.
	 * 
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings)
			throws Exception {
		System.out.println("TESymbolic.computeSignificance(): Not using the new orderings supplied");
		return teCalc.computeSignificance(newOrderings.length);
	}
}
