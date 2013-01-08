package infodynamics.utils;


/**
 * Calculators implementing this interface must provide a method to compute
 *  the statistical significance of their measurement, returning an analytically
 *  determined distribution.
 * 
 * @author Joseph Lizier
 *
 */
public interface AnalyticNullDistributionComputer {

	public AnalyticMeasurementDistribution computeSignificance() throws Exception;
}
