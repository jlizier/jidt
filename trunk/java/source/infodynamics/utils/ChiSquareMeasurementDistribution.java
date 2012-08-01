package infodynamics.utils;

/**
 * Class to represent analytic distributions of info theoretic measurements under
 * some null hypothesis of a relationship between the variables, where that
 * distribution is a Chi Square distribution.
 *
 * @author Joseph Lizier
 *
 */
public class ChiSquareMeasurementDistribution extends
		AnalyticMeasurementDistribution {

	protected int degreesOfFreedom;
	
	public ChiSquareMeasurementDistribution(double actualValue, int degreesOfFreedom) {
		super(actualValue, 1 - MathsUtils.chiSquareCdf(actualValue, degreesOfFreedom));
		this.degreesOfFreedom = degreesOfFreedom;
	}
}
