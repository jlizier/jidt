package infodynamics.utils;

/**
 * Generic class to represent analytic distributions of info theoretic measurements under
 * some null hypothesis of a relationship between the variables.
 * This class is not designed to be creatable - it's children (which represent
 * various concrete distributions) are those which should be created.
 * 
 * @author Joseph Lizier
 *
 */
public class AnalyticMeasurementDistribution extends MeasurementDistribution {

	protected AnalyticMeasurementDistribution(double actualValue, double pValue) {
		super(actualValue, pValue);
	}
}
