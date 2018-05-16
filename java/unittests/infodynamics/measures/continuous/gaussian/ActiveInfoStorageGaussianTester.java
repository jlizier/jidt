/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2018, Pedro Mediano, Joseph T. Lizier
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

import infodynamics.measures.continuous.ActiveInfoStorageCalculatorViaMutualInfo;
import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class ActiveInfoStorageGaussianTester extends TestCase {

	public void testAutoEmbeddingAIS() throws Exception {
		System.out.println("Start AIS autoembedding test (Gaussian).");
		// Generate multivariate data
		RandomGenerator rg = new RandomGenerator();
		double[] source = rg.generateNormalData(5000, 0, 1);
		
		for (int i=3; i < source.length; i++) {
			source[i] = 0.2*source[i-3] + 0.2*source[i-2] + 0.2*source[i-1] + source[i];
		}
		
		int correctK = 3;
		
		// Instantiate calculator and set search bounds
		ActiveInfoStorageCalculatorGaussian aisCalc = 
					new ActiveInfoStorageCalculatorGaussian();
		aisCalc.setProperty("k", "4");
		// Can't search larger than the max we expect because fluctuations can make the selection larger
		aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_K_SEARCH_MAX, "3");
		aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_TAU_SEARCH_MAX, "1");
		aisCalc.setProperty(ActiveInfoStorageCalculatorViaMutualInfo.PROP_AUTO_EMBED_METHOD,
				ActiveInfoStorageCalculatorViaMutualInfo.AUTO_EMBED_METHOD_MAX_CORR_AIS);
		aisCalc.setDebug(true);
		
		// Run optimisation
		aisCalc.initialise();
		aisCalc.setObservations(source);
		int optimisedK = Integer.parseInt(aisCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME));
		
		// Test that answer was correct
		assertEquals(correctK, optimisedK);
		
		double estimate = aisCalc.computeAverageLocalOfObservations();
		System.out.printf("Computed average was %.5f from %d samples\n", estimate, aisCalc.getNumObservations());
	}
}
