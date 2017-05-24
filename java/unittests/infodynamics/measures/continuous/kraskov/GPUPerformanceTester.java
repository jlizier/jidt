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

package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

import java.util.Calendar;
import junit.framework.TestCase;

public class GPUPerformanceTester extends TestCase {

  /**
   * Generic function to benchmark single evaluation on GPU vs CPU in
   * any dataset.
   *
   * @param source
   * @param dest
   * @param test_description
   */
  public void compareGPUPerformance(double[] source, double[] dest,
      String test_description) throws Exception {
    double startTime, cpu_duration, gpu_duration, cpu_val, gpu_val;
    MutualInfoCalculatorMultiVariateKraskov miCalc = 
      new MutualInfoCalculatorMultiVariateKraskov1();
    miCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    

    miCalc.setProperty("USE_GPU", "false");
    startTime = Calendar.getInstance().getTimeInMillis();
    miCalc.initialise(1,1);
    miCalc.setObservations(source, dest);
    cpu_val = miCalc.computeAverageLocalOfObservations();
    cpu_duration = Calendar.getInstance().getTimeInMillis() - startTime;

    miCalc.setProperty("USE_GPU", "true");
    miCalc.initialise(1,1);
    startTime = Calendar.getInstance().getTimeInMillis();
    miCalc.setObservations(source, dest);
    gpu_val = miCalc.computeAverageLocalOfObservations();
    gpu_duration = Calendar.getInstance().getTimeInMillis() - startTime;

    assertEquals(cpu_val, gpu_val, 0.0001);

    System.out.println("GPU Performance test: " + test_description);
    System.out.printf("CPU took %f ms, GPU took %f ms, speed ratio %f",
        cpu_duration, gpu_duration, cpu_duration/gpu_duration);

    return;
  }

  /**
   * Generic function to benchmark surrogates evaluation on GPU vs CPU in
   * any dataset.
   *
   * @param source
   * @param dest
   * @param test_description
   */
  public void compareGPUPerformanceSurrogates(double[] source, double[] dest,
      int nb_surrogates, String test_description) throws Exception {
    double startTime, cpu_duration, gpu_duration, cpu_val, gpu_val;
    MutualInfoCalculatorMultiVariateKraskov miCalc = 
      new MutualInfoCalculatorMultiVariateKraskov1();
    miCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    

    miCalc.setProperty("USE_GPU", "false");
    startTime = Calendar.getInstance().getTimeInMillis();
    miCalc.initialise(1,1);
    miCalc.setObservations(source, dest);
    cpu_val = miCalc.computeSignificance(nb_surrogates).getMeanOfDistribution();
    cpu_duration = Calendar.getInstance().getTimeInMillis() - startTime;

    miCalc.setProperty("USE_GPU", "true");
    miCalc.initialise(1,1);
    startTime = Calendar.getInstance().getTimeInMillis();
    miCalc.setObservations(source, dest);
    gpu_val = miCalc.computeSignificance(nb_surrogates).getMeanOfDistribution();
    gpu_duration = Calendar.getInstance().getTimeInMillis() - startTime;

    assertEquals(cpu_val, gpu_val, 0.0001);

    System.out.println("GPU Performance surrogate test: " + test_description);
    System.out.printf("CPU took %f ms, GPU took %f ms, speed ratio %f", cpu_duration, gpu_duration, cpu_duration/gpu_duration);

    return;
  }

  /**
   * Test in large, low-dimensional data.
   */
  public void testRandomLowDimension() throws Exception {

    MutualInfoCalculatorMultiVariateKraskov miCalc = 
      new MutualInfoCalculatorMultiVariateKraskov1();

    boolean gpuLoaded = true;
    try {
      miCalc.ensureCudaLibraryLoaded();
    } catch (Throwable e) {
      gpuLoaded = false;
    }

    // This will effectively ignore the test if GPU library not loaded properly
    if (!gpuLoaded) {
      return;
    }

    int timeSteps = 10000;
		RandomGenerator rg = new RandomGenerator();
		double[] source = rg.generateNormalData(timeSteps, 0, 1);
		double[] dest = rg.generateNormalData(timeSteps, 0, 1);

    compareGPUPerformance(source, dest, "Random low-dimensional data");

    return;
  }

}

