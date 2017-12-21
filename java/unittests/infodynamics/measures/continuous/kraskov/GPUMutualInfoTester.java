package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

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
    // miCalc.setDebug(true);

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


  /**
   * Test that GPU MI calculations with given reordering match the CPU version.
   *
   * Tests with the identity reordering (0,1,...N-1) and a random (but fixed)
   * reordering.
   *
   * Only runs if CUDA code has been precompiled.
   *
   * @throws Exception if something goes wrong
   */
  public void testGPUMutualInfoReordering() throws Exception {

    MutualInfoCalculatorMultiVariateKraskov miCalc = 
      new MutualInfoCalculatorMultiVariateKraskov1();
    // miCalc.setDebug(true);

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
    EmpiricalMeasurementDistribution cpu_dist, gpu_dist;

    // Test for low-dimensional data
		afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double [][] data = afr.getDouble2DMatrix();

    RandomGenerator rng = new RandomGenerator();
    int[] perm = rng.generateRandomPerturbations(data.length, 1)[0];
    int[][] newOrderings = new int[][] { MatrixUtils.range(0, data.length-1),
                                         perm};

    miCalc.setProperty("USE_GPU", "false");
    miCalc.initialise(1,1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                           MatrixUtils.selectColumn(data, 1));
    cpu_dist = miCalc.computeSignificance(newOrderings);

    miCalc.setProperty("USE_GPU", "true");
    miCalc.initialise(1,1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                           MatrixUtils.selectColumn(data, 1));
    gpu_dist = miCalc.computeSignificance(newOrderings);

    System.out.printf("GPU calculation on 2randomCols data. CPU got (%f,%f,%f)\n", cpu_dist.actualValue, cpu_dist.distribution[0], cpu_dist.distribution[1]);
    System.out.printf("GPU calculation on 2randomCols data. GPU got (%f,%f,%f)\n", gpu_dist.actualValue, gpu_dist.distribution[0], gpu_dist.distribution[1]);

    assertEquals(cpu_dist.actualValue,     cpu_dist.distribution[0], 0.00001);
    assertEquals(gpu_dist.actualValue,     gpu_dist.distribution[0], 0.00001);
    assertEquals(cpu_dist.actualValue,     gpu_dist.distribution[0], 0.00001);
    assertEquals(gpu_dist.actualValue,     cpu_dist.distribution[0], 0.00001);
  }

  /**
   * Test that GPU MI calculation with random surrogates gives reasonable
   * surrogate distributions and p-values.
   *
   * Only runs if CUDA code has been precompiled.
   *
   * @throws Exception if something goes wrong
   */
  public void testGPURandomPermutation() throws Exception {

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

    // Generate correlated Gaussian data
    RandomGenerator rng = new RandomGenerator();
    double[][] data = rng.generateBivariateNormalData(10000, 0, 1, 0, 1, 0.8);

    // miCalc.setDebug(true);
    miCalc.setProperty("USE_GPU", "true");
    miCalc.initialise(1, 1);
    miCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                           MatrixUtils.selectColumn(data, 1));
    EmpiricalMeasurementDistribution dist =
      miCalc.computeSignificance(10);

    System.out.printf("GPU random surrogates test: actual: %f, mean: %f, std: %f\n", dist.actualValue, dist.getMeanOfDistribution(), dist.getStdOfDistribution());
    assertTrue(dist.actualValue > dist.getMeanOfDistribution() + dist.getStdOfDistribution());

    return;
  }

  /**
   * Test that GPU MI calculation with a non-zero dynamic correlation
   * exclusion window gives correct results.
   *
   * Only runs if CUDA code has been precompiled.
   *
   * @throws Exception if something goes wrong
   */
  public void testGPUExclusionWindow() throws Exception {

    MutualInfoCalculatorMultiVariateKraskov miCalc = 
      new MutualInfoCalculatorMultiVariateKraskov1();
    // miCalc.setDebug(true);

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
    miCalc.setProperty("DYN_CORR_EXCL", "2");
    ArrayFileReader afr;
    double cpu_val, gpu_val;

    // Test for low-dimensional data
		afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double [][] data = afr.getDouble2DMatrix();

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

    System.out.printf("GPU calculation with theiler window. CPU got (%f)\n", cpu_val);
    System.out.printf("GPU calculation with theiler window. GPU got (%f)\n", gpu_val);

    assertEquals(cpu_val, cpu_val, 0.00001);

    return;
  }

}

