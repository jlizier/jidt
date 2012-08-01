package infodynamics.utils;

/**
 * <p>Structure to hold a distribution of info-theoretic measurements,
 *  and a significance value for how an original measurement compared 
 *  with these.</p>
 * 
 * <p>While in theory this class could be directly used, it is it's
 * children which are intended to be used.</p>
 * 
 * @author Joseph Lizier
 *
 */
public class MeasurementDistribution {
	/**
	 * Actual observed value of the measurement
	 */
	public double actualValue;
	/**
	 * Probability that surrogate measurement is greater than
	 *  the observed value.
	 * (Small pValue means that the observed value is highly significant).
	 */
	public double pValue;

	/**
	 * Allow empty constructor for internal use only when actualValue and
	 * pValue will be supplied later
	 */
	protected MeasurementDistribution() {
	}
	
	/**
	 * Construct with supplied actual value and p-value for it.
	 * 
	 * @param actualValue actual observed value
	 * @param pValue p-value that the surrogate measurement is larger
	 *  than the observed value
	 */
	public MeasurementDistribution(double actualValue, double pValue) {
		this.actualValue = actualValue;
		this.pValue = pValue;
	}
}
