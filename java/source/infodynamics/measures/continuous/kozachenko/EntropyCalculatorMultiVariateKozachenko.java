package infodynamics.measures.continuous.kozachenko;

import infodynamics.measures.continuous.EntropyCalculatorMultiVariate;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MathsUtils;

/**
 * Compute the entropy using the Kozachenko estimation method.
 * See:
 * - "A statistical estimate for the entropy of a random vector"
 * 		Kozachenko, L., Leonenko, N.,
 * 		Problems of Information Transmission, 23 (1987) 9–16
 * - "Estimating mutual information"
 * 		Kraskov, A., Stogbauer, H., Grassberger, P.,
 * 		Physical Review E 69, (2004) 066138
 * - "Measuring Global Behaviour of Multi-Agent Systems from
 *  	Pair-Wise Mutual Information",
 *  	George Mathews, Hugh Durrant-Whyte, and Mikhail Prokopenko
 * 
 * This class computes it exactly as in "Estimating mutual information", i.e. using natural
 *  units and twice the minimum distance.
 *  Implementing this to check if our other implementation was correct.
 * 
 * @author Joseph Lizier
 *
 */
public class EntropyCalculatorMultiVariateKozachenko  
	implements EntropyCalculatorMultiVariate {

	protected boolean debug = false;
	private int totalObservations;
	private int dimensions;
	protected double[][] rawData;
	private double lastEntropy;
	private double[] lastLocalEntropy;
	private boolean isComputed;
	
	public static final double EULER_MASCHERONI_CONSTANT = 0.5772156;
	
		
	public EntropyCalculatorMultiVariateKozachenko() {
		totalObservations = 0;
		dimensions = 0;
		isComputed = false;
		lastLocalEntropy = null;
	}

	public void initialise(int dimensions) {
		this.dimensions = dimensions;
		rawData = null;
		totalObservations = 0;
		isComputed = false;
		lastLocalEntropy = null;
	}

	public void setObservations(double[][] observations) {
		rawData = observations;
		totalObservations = observations.length;
		isComputed = false;
		lastLocalEntropy = null;
	}

	/**
	 * Each row of the data is an observation; each column of
	 *  the row is a new variable in the multivariate observation.
	 * This method signature allows the user to call setObservations for
	 *  joint time series without combining them into a single joint time
	 *  series (we do the combining for them).
	 * 
	 * @param data1
	 * @param data2
	 * @throws Exception When the length of the two arrays of observations do not match.
	 */
	public void setObservations(double[][] data1,
			double[][] data2) throws Exception {
		int timeSteps = data1.length;
		if ((data1 == null) || (data2 == null)) {
			throw new Exception("Cannot have null data arguments");
		}
		if (data1.length != data2.length) {
			throw new Exception("Length of data1 (" + data1.length + ") is not equal to the length of data2 (" +
					data2.length + ")");
		}
		int data1Variables = data1[0].length;
		int data2Variables = data2[0].length;
		double[][] data = new double[timeSteps][data1Variables + data2Variables];
		for (int t = 0; t < timeSteps; t++) {
			System.arraycopy(data1[t], 0, data[t], 0, data1Variables);
			System.arraycopy(data2[t], 0, data[t], data1Variables, data2Variables);
		}
		// Now defer to the normal setObservations method
		setObservations(data);
	}

	/**
	 * Computes average entropy of previously provided observations.
	 * 
	 * @return entropy in natural units
	 */
	public double computeAverageLocalOfObservations() {
		if (isComputed) {
			return lastEntropy;
		}
		double sdTermHere = sdTerm(totalObservations, dimensions);
		double emConstHere = eulerMacheroniTerm(totalObservations);
		double[] minDistance = EuclideanUtils.computeMinEuclideanDistances(rawData);
		double entropy = 0.0;
		if (debug) {
			System.out.println("t,\tminDist,\tlogMinDist,\tsum");
		}
		for (int t = 0; t < rawData.length; t++) {
			entropy += Math.log(2.0 * minDistance[t]);
			if (debug) {
				System.out.println(t + ",\t" + 
						minDistance[t] + ",\t" +
						Math.log(minDistance[t]) + ",\t" +
						entropy);
			}
		}
		// Using natural units
		// entropy /= Math.log(2);
		entropy *= (double) dimensions / (double) totalObservations;
		if (debug) {
			System.out.println("Sum part:   " + entropy);
			System.out.println("Euler part: " + emConstHere);
			System.out.println("Sd term:    " + sdTermHere);
		}
		entropy += emConstHere;
		entropy += sdTermHere;
		lastEntropy = entropy;
		isComputed = true;
		return entropy;
	}

	/**
	 * Computes local entropies of given values, using previously provided observations.
	 * 
	 * @return local entropies in natural units
	 */
	public double[] computeLocalOfPreviousObservations() {
		if (lastLocalEntropy != null) {
			return lastLocalEntropy;
		}

		double sdTermHere = sdTerm(totalObservations, dimensions);
		double emConstHere = eulerMacheroniTerm(totalObservations);
		double constantToAddIn = sdTermHere + emConstHere;
		
		double[] minDistance = EuclideanUtils.computeMinEuclideanDistances(rawData);
		double entropy = 0.0;
		double[] localEntropy = new double[rawData.length];
		if (debug) {
			System.out.println("t,\tminDist,\tlogMinDist,\tlocal,\tsum");
		}
		for (int t = 0; t < rawData.length; t++) {
			localEntropy[t] = Math.log(2.0 * minDistance[t]) * (double) dimensions;
			// using natural units
			// localEntropy[t] /= Math.log(2);
			localEntropy[t] += constantToAddIn;
			entropy += localEntropy[t];
			if (debug) {
				System.out.println(t + ",\t" + 
						minDistance[t] + ",\t" +
						Math.log(minDistance[t]) + ",\t" +
						localEntropy[t] + ",\t" +
						entropy);
			}
		}
		entropy /= (double) totalObservations;
		lastEntropy = entropy;
		lastLocalEntropy = localEntropy;
		return localEntropy;
	}

	public double[] computeLocalUsingPreviousObservations(double[][] states) throws Exception {
		throw new Exception("Local method for other data not implemented");
	}

	public double[] computeLocalUsingPreviousObservations(double[][] states1, double[][] states2) throws Exception {
		throw new Exception("Local method for other data not implemented");
	}

	/**
	 * Returns the value of the Euler-Mascheroni term.
	 * Public for debugging purposes
	 * 
	 * @return
	 */
	public double eulerMacheroniTerm(int N) {
		// Using natural units
		// return EULER_MASCHERONI_CONSTANT / Math.log(2);
		try {
			return -MathsUtils.digamma(1) + MathsUtils.digamma(N);
		} catch (Exception e) {
			// Exception will only be thrown if N < 0
			return 0;
		}
	}
	
	/**
	 * Returns the value of the Sd term
	 * Public for debugging purposes
	 * 
	 * @param numObservations
	 * @param dimensions
	 * @return
	 */
	public double sdTerm(int numObservations, int dimensions) {
		// To compute directly:
		// double unLoggedSdTerm = 
		//	Math.pow(Math.PI/4.0, ((double) dimensions) / 2.0) /	// Brought 2^d term from denominator into Pi term
		//	MathsUtils.gammaOfArgOn2Plus1(dimensions);
		
		// But we need to compute it carefully, to allow the maximum range of dimensions
		// double unLoggedSdTerm = 
		//	1.0 / MathsUtils.gammaOfArgOn2Plus1IncludeDivisor(dimensions, 
		//			Math.pow(Math.PI, ((double) dimensions) / 2.0));
		// Don't include the 2^d in the above divisor, since that makes the divisor < 1, which 
		//  doesn't help at all.
		// unLoggedSdTerm /= Math.pow(2, dimensions);
		// return Math.log(unLoggedSdTerm) / Math.log(2);
		// Using natural units
		// return Math.log(unLoggedSdTerm);
		
		// But even that method falls over by about d = 340.
		// Break down the log into the log of a factorial term and the log of a constant term
		double constantTerm = Math.pow(Math.PI / 4.0, (double) dimensions / 2.0);
		double result = 0.0;
		if (dimensions % 2 == 0) {
			// d even
			// Now take log (1/(d/2)!) = -log (1/(d/2)!) = -sum(d/2 --) {log d/2}
			for (int d = dimensions/2; d > 1; d--) {
				result -= Math.log(d);
			}
		} else {
			// d odd
			constantTerm *= Math.pow(2.0, (double) (dimensions + 1) / 2.0);
			constantTerm /= Math.sqrt(Math.PI);
			// Now take log (1/d!!) = - log (d!!) = - sum(d -= 2) {log d}
			for (int d = dimensions; d > 1; d -= 2) {
				result -= Math.log(d);
			}
		}
		result += Math.log(constantTerm);
		return result;
	}
	
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return lastEntropy;
	}

	public int getNumObservations() {
		return totalObservations;
	}
}
