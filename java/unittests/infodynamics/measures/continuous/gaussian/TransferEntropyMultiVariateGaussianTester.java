/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
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

package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariate;
import infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariateViaCondMutualInfo;
import infodynamics.utils.RandomGenerator;

public class TransferEntropyMultiVariateGaussianTester
	extends infodynamics.measures.continuous.TransferEntropyMultiVariateAbstractTester {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testLocalsAverageCorrectly() throws Exception {
		
		TransferEntropyCalculatorMultiVariateGaussian teCalc =
				new TransferEntropyCalculatorMultiVariateGaussian();
		
		super.testLocalsAverageCorrectly(teCalc, 2, 100, 1);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		TransferEntropyCalculatorMultiVariateGaussian teCalc =
				new TransferEntropyCalculatorMultiVariateGaussian();
		
		super.testComputeSignificanceDoesntAlterAverage(teCalc, 2, 100, 1);
	}

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testUnivariateSignatureMatchesMultivariate() throws Exception {
		
		TransferEntropyCalculatorMultiVariateGaussian teCalc =
				new TransferEntropyCalculatorMultiVariateGaussian();
		
		super.testUnivariateMatchesMultivariateRoute(teCalc, 100, 1);
	}

	  public void testAutoEmbeddingAIS() throws Exception {
		    System.out.println("Start AIS autoembedding test.");
		    // Generate multivariate data
				RandomGenerator rg = new RandomGenerator();
				double[][] source = rg.generateNormalData(10000, 2, 0, 1);
				double[][] target = rg.generateNormalData(10000, 2, 0, 1);

		    for (int i=3; i < source.length; i++) {
		      source[i][0] = 0.2*source[i-3][0] + 0.2*source[i-2][0] + 0.2*source[i-1][0] + source[i][0];
		      source[i][1] = 0.2*source[i-3][1] + 0.2*source[i-2][1] + 0.2*source[i-1][1] + source[i][1];

		      target[i][0] = 0.2*target[i-3][0] + 0.2*target[i-2][0] + 0.2*target[i-1][0] + target[i][0];
		      target[i][1] = 0.2*target[i-3][1] + 0.2*target[i-2][1] + 0.2*target[i-1][1] + target[i][1];
		    }

		    int correctK = 3;
		    int correctL = 3;

		    // Instantiate calculator and set search bounds
				TransferEntropyCalculatorMultiVariateGaussian teCalc = 
						new TransferEntropyCalculatorMultiVariateGaussian();
				teCalc.setProperty("k", "4");
				teCalc.setProperty(teCalc.PROP_K_SEARCH_MAX, "3");
				teCalc.setProperty(teCalc.PROP_TAU_SEARCH_MAX, "1");
				teCalc.setProperty(teCalc.PROP_AUTO_EMBED_METHOD, teCalc.AUTO_EMBED_METHOD_MAX_CORR_AIS);

		    // Run optimisation
				teCalc.initialise(2, 2);
				teCalc.setObservations(source, target);
				int optimisedK = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorMultiVariate.K_PROP_NAME));
				int optimisedL = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorMultiVariateViaCondMutualInfo.L_PROP_NAME));

		    // Test that answer was correct
				assertEquals(correctK, optimisedK);
				assertEquals(correctL, optimisedL);
		  }

		  public void testAutoEmbeddingTE() throws Exception {
		    System.out.println("Start AIS+TE autoembedding test.");

		    // Generate multivariate data (note that source time series has no memory)
				RandomGenerator rg = new RandomGenerator();
				double[][] source = rg.generateNormalData(10000, 2, 0, 1);
				double[][] target = rg.generateNormalData(10000, 2, 0, 1);

		    for (int i=3; i < source.length; i++) {
		      target[i][0] = 0.2*source[i-2][0] + 0.2*source[i-1][0] + 
		                     0.2*target[i-2][0] + 0.2*target[i-1][0] + target[i][0];

		      target[i][1] = 0.2*source[i-2][1] + 0.2*source[i-1][1] + 
		                     0.2*target[i-2][1] + 0.2*target[i-1][1] + target[i][1];
		    }

		    int correctK = 2;
		    int correctL = 2;

		    // Instantiate calculator and set search bounds
		    TransferEntropyCalculatorMultiVariateGaussian teCalc = 
						new TransferEntropyCalculatorMultiVariateGaussian();
				teCalc.setProperty("k", "4");
				teCalc.setProperty(teCalc.PROP_K_SEARCH_MAX, "2");
				teCalc.setProperty(teCalc.PROP_TAU_SEARCH_MAX, "1");
				teCalc.setProperty(teCalc.PROP_AUTO_EMBED_METHOD, teCalc.AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE);

		    // Run optimisation
				teCalc.initialise(2, 2);
				teCalc.setObservations(source, target);
				int optimisedK = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorMultiVariate.K_PROP_NAME));
				int optimisedL = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorMultiVariateViaCondMutualInfo.L_PROP_NAME));

		    // Test that answer was correct
				assertEquals(correctK, optimisedK);
				assertEquals(correctL, optimisedL);
		  }
}
