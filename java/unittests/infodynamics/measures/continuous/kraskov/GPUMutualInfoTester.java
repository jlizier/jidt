package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;

import junit.framework.TestCase;

public class GPUMutualInfoTester extends TestCase {

  /**
   * Test that GPU MI calculations match the CPU version.
   *
   * Only runs if CUDA code has been precompiled.
   *
   * @throws Exception if something goes wrong
   */
  public void testGPUMutualInfo() throws Exception {

    MutualInfoCalculatorMultiVariateKraskov miCalc = 
      new MutualInfoCalculatorMultiVariateKraskov1();
    miCalc.setDebug(true);

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

    // Now starts the actual test. Set up variables and properties
    miCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    ArrayFileReader afr;
    double[][] data;
    double cpu_val, gpu_val;

    // Test for low-dimensional data
		afr = new ArrayFileReader("demos/data/2randomCols-1.txt");
		data = afr.getDouble2DMatrix();

    miCalc.setProperty("USE_GPU", "false");
    miCalc.initialise(1,1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                           MatrixUtils.selectColumn(data, 1));
    cpu_val = miCalc.computeAverageLocalOfObservations();

    miCalc.setProperty("USE_GPU", "true");
    miCalc.initialise(1,1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                           MatrixUtils.selectColumn(data, 1));
    gpu_val = miCalc.computeAverageLocalOfObservations();

    System.out.printf("GPU calculation on 2randomCols data. Expected: %f. Got %f\n", cpu_val, gpu_val);
    assertEquals(cpu_val, gpu_val, 0.00001);

    // Test for slightly higher-dimensional data
		afr = new ArrayFileReader("demos/data/4ColsPairedNoisyDependence-1.txt");
		data = afr.getDouble2DMatrix();

    miCalc.setProperty("USE_GPU", "false");
    miCalc.initialise(2,2);
    miCalc.setObservations(MatrixUtils.selectColumns(data, new int[]{0,1}),
                           MatrixUtils.selectColumns(data, new int[]{2,3}));
    cpu_val = miCalc.computeAverageLocalOfObservations();

    miCalc.setProperty("USE_GPU", "true");
    miCalc.initialise(2,2);
    miCalc.setObservations(MatrixUtils.selectColumns(data, new int[]{0,1}),
                           MatrixUtils.selectColumns(data, new int[]{2,3}));
    gpu_val = miCalc.computeAverageLocalOfObservations();

    System.out.printf("GPU calculation on 4ColsPaired data. Expected: %f. Got %f\n", cpu_val, gpu_val);
    assertEquals(cpu_val, gpu_val, 0.00001);

  }

  /**
   * Test that GPU local MI calculations match the CPU version.
   *
   * Only runs if CUDA code has been precompiled.
   *
   * @throws Exception if something goes wrong
   */
  public void testGPULocalMutualInfo() throws Exception {

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

    // Now starts the actual test. Set up variables and properties
    miCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    ArrayFileReader afr;
    double[][] data;
    double[] cpu_vals, gpu_vals;

    // Test for low-dimensional data
		afr = new ArrayFileReader("demos/data/2randomCols-1.txt");
		data = afr.getDouble2DMatrix();

    miCalc.setProperty("USE_GPU", "false");
    miCalc.initialise(1,1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                           MatrixUtils.selectColumn(data, 1));
    cpu_vals = miCalc.computeLocalOfPreviousObservations();

    miCalc.setProperty("USE_GPU", "true");
    miCalc.initialise(1,1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                           MatrixUtils.selectColumn(data, 1));
    gpu_vals = miCalc.computeLocalOfPreviousObservations();

    for (int i = 0; i < data.length; i++) {
      assertEquals(cpu_vals[i], gpu_vals[i], 0.00001);
    }

    // Test for slightly higher-dimensional data
		afr = new ArrayFileReader("demos/data/4ColsPairedNoisyDependence-1.txt");
		data = afr.getDouble2DMatrix();

    miCalc.setProperty("USE_GPU", "false");
    miCalc.initialise(2,2);
    miCalc.setObservations(MatrixUtils.selectColumns(data, new int[]{0,1}),
                           MatrixUtils.selectColumns(data, new int[]{2,3}));
    cpu_vals = miCalc.computeLocalOfPreviousObservations();

    miCalc.setProperty("USE_GPU", "true");
    miCalc.initialise(2,2);
    miCalc.setObservations(MatrixUtils.selectColumns(data, new int[]{0,1}),
                           MatrixUtils.selectColumns(data, new int[]{2,3}));
    gpu_vals = miCalc.computeLocalOfPreviousObservations();

    for (int i = 0; i < data.length; i++) {
      assertEquals(cpu_vals[i], gpu_vals[i], 0.00001);
    }

  }
}
