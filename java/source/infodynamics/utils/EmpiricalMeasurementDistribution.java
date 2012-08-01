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
public class EmpiricalMeasurementDistribution extends MeasurementDistribution {

	/**
	 * Distribution of surrogate measurement values
	 */
	public double[] distribution;
	/**
	 * Whether the mean of the surrogate measurement distribution has
	 * been computed
	 */
	protected boolean computedMean = false;
	/**
	 * Computed mean of the surrogate measurement distribution
	 */
	protected double meanOfDist;
	/**
	 * Computed mean of the surrogate measurement distribution
	 */
	protected double stdOfDist;
	
	public EmpiricalMeasurementDistribution(int size) {
		super(); // Creating the super class with mean and pValue 0
		// These value will be filled out by the caller later.
		distribution = new double[size];
	}

	public EmpiricalMeasurementDistribution(double[] distribution, double actualValue) {
		super(actualValue, 0); // Using pValue = 0 temporarily ...
		this.distribution = distribution;
		int countWhereActualIsNotGreater = 0;
		for (int i = 0; i < distribution.length; i++) {
			if (distribution[i] >= actualValue) {
				countWhereActualIsNotGreater++;
			}
		}
		pValue = (double) countWhereActualIsNotGreater / (double) distribution.length;
	}
	
	// TODO Compute the significance under the assumption of a Gaussian distribution
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
