package infodynamics.utils;

/**
 * 
 * Structure to hold a distribution of info-theoretic measurements,
 *  and a significance value for how an original measurement compared 
 *  with these.
 * 
 * @author Joseph Lizier
 *
 */
public class MeasurementDistribution {

	/**
	 * Distribution of surrogate measurement values
	 */
	public double[] distribution;
	/**
	 * Actual observed value of the measurement
	 */
	public double actualValue;
	/**
	 * Probability that surrogate measurement is greater than
	 *  the observed value
	 */
	public double pValue;
	
	protected boolean computedMean = false;
	protected double meanOfDist;
	protected double stdOfDist;
	
	public MeasurementDistribution(int size) {
		distribution = new double[size];
	}

	public MeasurementDistribution(double[] distribution, double actualValue) {
		this.actualValue = actualValue;
		this.distribution = distribution;
		int countWhereActualIsNotGreater = 0;
		for (int i = 0; i < distribution.length; i++) {
			if (distribution[i] >= actualValue) {
				countWhereActualIsNotGreater++;
			}
		}
		pValue = (double) countWhereActualIsNotGreater / (double) distribution.length;
	}
	
	/*
	public double computeGaussianSignificance() {
		// Need to conpute the significance based on the assumption of 
		//  an underlying Gaussian distribution.
		// Use the t distribution for analysis, since we have a finite
		//  number of samples to comptue the mean and std from.
		return 0;
	}
	*/
	
	public double getTSscore() {
		if (! computedMean) {
			meanOfDist = MatrixUtils.mean(distribution);
			stdOfDist = MatrixUtils.stdDev(distribution, meanOfDist);
			computedMean = true;
		}
		double t = (actualValue - meanOfDist) / stdOfDist;
		return t;
	}
	
	public double getMeanOfDistribution() {
		if (! computedMean) {
			meanOfDist = MatrixUtils.mean(distribution);
			stdOfDist = MatrixUtils.stdDev(distribution, meanOfDist);
			computedMean = true;
		}
		return meanOfDist;
	}

	public double getStdOfDistribution() {
		if (! computedMean) {
			meanOfDist = MatrixUtils.mean(distribution);
			stdOfDist = MatrixUtils.stdDev(distribution, meanOfDist);
			computedMean = true;
		}
		return stdOfDist;
	}
}
