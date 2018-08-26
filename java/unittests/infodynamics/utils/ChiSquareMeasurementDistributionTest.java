/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2017, Joseph T. Lizier
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

package infodynamics.utils;

import junit.framework.TestCase;

/**
 * Test functionality of the utility functions in ChiSquareMeasurementDistribution
 * 
 * @author Joseph Lizier.
 *
 */
public class ChiSquareMeasurementDistributionTest extends TestCase {

	public void testPValueAndEstimateCorrespondence() {
		double[] actualValues = new double[] {0.000001, 0.000002, 0.000005, 0.00001, 0.00002, 0.00005,
					0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
		int[] numObs = new int[] {100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000};
		int[] degFreedom = new int[] {1, 2, 3, 5, 10, 20};
		for (int a = 0; a < actualValues.length; a++) {
			for (int d = 0; d < degFreedom.length; d++) {
				for (int n = 0; n < numObs.length; n++) {
					ChiSquareMeasurementDistribution csmDist = new ChiSquareMeasurementDistribution(actualValues[a], numObs[n], degFreedom[d]);
					System.out.printf("N=%d, degFree=%d, measuredValue=%.6f: pValue=%.6f\n", numObs[n], degFreedom[d], actualValues[a], csmDist.pValue);
					// Check that we get the same p value if we stick this estimate in afterwards:
					assertEquals(csmDist.pValue, csmDist.computePValueForGivenEstimate(actualValues[a]), 0.000001);
					// Check that if we ask for an inverse on this pValue we get the same estimate back:
					// Can easily get rounding errors here when the pValue approaches 1 or 0, so we won't test these for now:
					if ((csmDist.pValue < 0.99) && (csmDist.pValue > 0.000001)) {
						assertEquals(actualValues[a], csmDist.computeEstimateForGivenPValue(csmDist.pValue), 0.000001);
					}
				}
			}
		}
	}
	
	public void testAnalyticMeanAndStd() {
		int numValues = 2000;
		double[] actualValues = new double[] {0.000001, 0.000002};
		int[] numObs = new int[] {100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000};
		int[] degFreedom = new int[] {1, 2, 3, 5, 10, 20};
		for (int a = 0; a < actualValues.length; a++) {
			for (int d = 0; d < degFreedom.length; d++) {
				for (int n = 0; n < numObs.length; n++) {
					ChiSquareMeasurementDistribution csmDist = new ChiSquareMeasurementDistribution(actualValues[a], numObs[n], degFreedom[d]);
					// System.out.printf("N=%d, degFree=%d, measuredValue=%.6f\n", numObs[n], degFreedom[d], actualValues[a]);
					// Compute the mean of the null measurement distribution analytically and semi-empirically:
					double meanOfDist = 0.0, meanSqrs = 0.0;
					for (int i = 0; i < numValues; i++) {
						double contribution = csmDist.computeEstimateForGivenPValue((1.0*i+0.5)/(double)numValues);
						meanOfDist += contribution;
						meanSqrs += contribution * contribution;
					}
					meanOfDist /= numValues;
					// Check that we get the same mean value if we call the analytic method:
					assertEquals(meanOfDist, csmDist.getMeanOfDistribution(), 0.00001);
					meanSqrs /= numValues;
					double std = Math.sqrt(meanSqrs - meanOfDist*meanOfDist);
					// check that we get the same std deviation if we call the analytic method:
					// This one has a little more numerical error in the empirical one, so we test it
					//  to less decimal places.
					assertEquals(std, csmDist.getStdOfDistribution(), 0.0001);
				}
			}
		}
		
	}
	
	public void testBiasCorrection() {
		double[] actualValues = new double[] {0.000001, 0.000002, 0.000005, 0.00001, 0.00002, 0.00005,
					0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
		int[] numObs = new int[] {100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000};
		int[] degFreedom = new int[] {1, 2, 3, 5, 10, 20};
		for (int a = 0; a < actualValues.length; a++) {
			for (int d = 0; d < degFreedom.length; d++) {
				for (int n = 0; n < numObs.length; n++) {
					// Uncorrected, just to get mean:
					ChiSquareMeasurementDistribution csmDist = new ChiSquareMeasurementDistribution(actualValues[a], numObs[n], degFreedom[d]);
					double meanOfUncorrected = csmDist.getMeanOfDistribution();
					// Now use corrected:
					ChiSquareMeasurementDistribution csmDistCorr = new ChiSquareMeasurementDistribution(actualValues[a] - meanOfUncorrected, numObs[n], degFreedom[d], true);
					assertEquals(csmDist.pValue, csmDistCorr.pValue, 0.0000001);
					// Check that we get the same p value if we stick this estimate in afterwards:
					assertEquals(csmDistCorr.pValue, csmDistCorr.computePValueForGivenEstimate(actualValues[a] - meanOfUncorrected), 0.000001);
					// Check that if we ask for an inverse on this pValue we get the same estimate back:
					// Can easily get rounding errors here when the pValue approaches 1 or 0, so we won't test these for now:
					if ((csmDist.pValue < 0.99) && (csmDist.pValue > 0.000001)) {
						assertEquals(actualValues[a] - meanOfUncorrected, csmDistCorr.computeEstimateForGivenPValue(csmDistCorr.pValue), 0.000001);
					}
				}
			}
		}
	}
}
