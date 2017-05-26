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

package infodynamics.demos;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import infodynamics.utils.ParsedProperties;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;

import java.lang.Math;

/**
 * = Example 10 - GPU benchmark script =
 *
 * This class is used to demonstrate how the GPU module is activated and
 * includes a benchmark to test how fast it is, in comparison with its CPU
 * counterpart.
 *
 * If run from Linux, the script <code>plotExample10BenchmarkResults.py</code>
 * will use the Python library matplotlib to plot the results of the benchmark.
 *
 * @author Pedro AM Mediano
 *
 */
public class Example10GPUBenchmark {


  /**
   * Calculate MI for the given multivariate time series using the Kraskov
   * algorithm, either with CPU or GPU implementation.
   *
   * @param src Source variable for MI calculation
   * @param tgt Target variable for MI calculation
   * @param useGPU Whether or not to use the GPU for this calculation
   * @return An array with two doubles. The first component is the time in ms
   * it took to run the calculation and the second component is the actual
   * value of the calculated MI.
   */
  public static double[] runEstimator(double[][] src, double[][] tgt, boolean useGPU) throws Exception {
    MutualInfoCalculatorMultiVariate miCalc = new MutualInfoCalculatorMultiVariateKraskov1();
    if (useGPU) 
      miCalc.setProperty("USE_GPU", "true");
    miCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    miCalc.setProperty("NORMALISE", "false");
    miCalc.setProperty("k", "4");
    miCalc.initialise(src[0].length, tgt[0].length);
    miCalc.setObservations(src, tgt);
    double[] timeAndValue = new double[2];
    long startTime = System.nanoTime();
    timeAndValue[1] = miCalc.computeAverageLocalOfObservations();
    timeAndValue[0] = (System.nanoTime() - startTime)/1000000.0;
    return timeAndValue;
  }

  /**
   * Benchmark on a small predefined dataset, for debugging purposes.
   */
  public static void Benchmark0() throws Exception {
    double[][] src = new double[][] {{0}, {1}, {2}, {3}, {4}, {5}, {6}};
    double[][] tgt = new double[][] {{0}, {0}, {0}, {0}, {0}, {0}, {0}};
    System.out.printf("CPU value: %f\n", runEstimator(src, tgt, false)[1]);
    System.out.printf("GPU value: %f\n", runEstimator(src, tgt, true)[1]);

  }

  /**
   * Benchmark on white Gaussian data of arbitrary dimension.
   */
  public static void Benchmark1(String filename) throws Exception {

    boolean append_to_file = false;
    FileWriter write = new FileWriter(filename, append_to_file);
    PrintWriter print_line = new PrintWriter(write);

    int[] N_vec = new int[] {500, 1000, 2000, 4000};
    int[] D_vec = new int[] {1, 3, 5};

    RandomGenerator rg = new RandomGenerator();
    int nb_repetitions = 5;

    for (int i = 0; i < N_vec.length; i++) {
      for (int j = 0; j < D_vec.length; j++) {

        int N = N_vec[i];
        int D = D_vec[j];

        double[][] src = rg.generateNormalData(N, D, 0, 1);
        double[][] tgt = rg.generateNormalData(N, D, 0, 1);

        double cpuTime = 0, gpuTime = 0;
        for (int r = 0; r < nb_repetitions; r++) {
          // CPU calculation
          double[] cpuVals = runEstimator(src, tgt, false);
          cpuTime += cpuVals[0];

          // GPU calculation
          double[] gpuVals = runEstimator(src, tgt, true);
          gpuTime += gpuVals[0];

          if (Math.abs(cpuVals[1] - gpuVals[1]) > 1e-4) {
            System.out.printf("Values differ. CPU: %f, GPU: %f\n", cpuVals[1], gpuVals[1]);
          }

        }

        print_line.printf("%d\t%d\t%f\t%f\n", N, D,
          cpuTime/nb_repetitions, gpuTime/nb_repetitions);

      }
    }

    print_line.close();

    return;
  }



  /**
   * Benchmark on 2D correlated Gaussian.
   *
   * Reference:
   *
   * http://math.stackexchange.com/questions/446093/generate-correlated-normal-random-variables 
   *
   */
  public static void Benchmark2(String filename) throws Exception {

    boolean append_to_file = false;
    FileWriter write = new FileWriter(filename, append_to_file);
    PrintWriter print_line = new PrintWriter(write);

    int[] N_vec = new int[] {500, 1000, 2000, 4000};
    double[] corr_vec = new double[] {0, 0.2, 0.4, 0.6};

    RandomGenerator rg = new RandomGenerator();
    int nb_repetitions = 5;

    for (int i = 0; i < N_vec.length; i++) {
      for (int j = 0; j < corr_vec.length; j++) {

        int N = N_vec[i];
        double corr = corr_vec[j];

        double[] r1 = rg.generateNormalData(N, 0, 1);
        double[] r2 = rg.generateNormalData(N, 0, 1);
        double[][] src = new double[N][1];
        double[][] tgt = new double[N][1];
        for (int t = 0; t < N; t++) {
          src[t][0] = r1[t];
          tgt[t][0] = corr*src[t][0] + Math.sqrt(1 - corr*corr)*r2[t];
        }

        double cpuTime = 0, gpuTime = 0;
        for (int r = 0; r < nb_repetitions; r++) {
          // CPU calculation
          double[] cpuVals = runEstimator(src, tgt, false);
          cpuTime += cpuVals[0];

          // GPU calculation
          double[] gpuVals = runEstimator(src, tgt, true);
          gpuTime += gpuVals[0];

          if (Math.abs(cpuVals[1] - gpuVals[1]) > 1e-4) {
            System.out.printf("Values differ. CPU: %f, GPU: %f\n", cpuVals[1], gpuVals[1]);
          }
        }

        // Save time in milliseconds
        print_line.printf("%d\t%f\t%f\t%f\n", N, corr,
            cpuTime/nb_repetitions, gpuTime/nb_repetitions);

      }
    }

    print_line.close();

    return;
  }

  /**
   * Run several of the benchmarks above and write the results to a
   * .txt file.
   * 
   * @param args List of file names to save the results of the benchmark.
   */
  public static void main(String[] args) throws Exception {
    // Benchmark0();
    Benchmark1(args[0]);
    Benchmark2(args[1]);
    return;
  }

}
