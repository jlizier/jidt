package infodynamics.measures.continuous.kraskov;

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

import junit.framework.TestCase;

public class GPUCondMutualInfoTester extends TestCase {

  /**
   * Test that GPU CMI calculations match the CPU version.
   *
   * Only runs if CUDA code has been precompiled.
   *
   * @throws Exception if something goes wrong
   */
  public void testGPUCondMutualInfo() throws Exception {

    ConditionalMutualInfoCalculatorMultiVariateKraskov cmiCalc = 
      new ConditionalMutualInfoCalculatorMultiVariateKraskov1();
    // cmiCalc.setDebug(true);

    boolean gpuLoaded = true;
    try {
      cmiCalc.ensureCudaLibraryLoaded();
    } catch (Throwable e) {
      gpuLoaded = false;
    }

    // This will effectively ignore the test if GPU library not loaded properly
    if (!gpuLoaded) {
      return;
    }

    // Now starts the actual test. Set up variables and properties
    cmiCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    ArrayFileReader afr;
    double[][] data;
    double cpu_val, gpu_val;

    // Test for low-dimensional data
		afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		data = afr.getDouble2DMatrix();

    cmiCalc.setProperty("USE_GPU", "false");
    cmiCalc.initialise(1,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    cpu_val = cmiCalc.computeAverageLocalOfObservations();

    cmiCalc.setProperty("USE_GPU", "true");
    cmiCalc.initialise(1,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    gpu_val = cmiCalc.computeAverageLocalOfObservations();

    System.out.printf("GPU CMI calculation on 4randomCols data. Expected: %f. Got %f\n", cpu_val, gpu_val);
    assertEquals(cpu_val, gpu_val, 0.00001);

    // Test for slightly higher-dimensional data
		afr = new ArrayFileReader("demos/data/4ColsPairedNoisyDependence-1.txt");
		data = afr.getDouble2DMatrix();

    cmiCalc.setProperty("USE_GPU", "false");
    cmiCalc.initialise(2,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumns(data, new int[]{0,1}),
                            MatrixUtils.selectColumns(data, new int[]{2}),
                            MatrixUtils.selectColumns(data, new int[]{3}));
    cpu_val = cmiCalc.computeAverageLocalOfObservations();

    cmiCalc.setProperty("USE_GPU", "true");
    cmiCalc.initialise(2,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumns(data, new int[]{0,1}),
                            MatrixUtils.selectColumns(data, new int[]{2}),
                            MatrixUtils.selectColumns(data, new int[]{3}));
    gpu_val = cmiCalc.computeAverageLocalOfObservations();

    System.out.printf("GPU CMI calculation on 4ColsPaired data. Expected: %f. Got %f\n", cpu_val, gpu_val);
    assertEquals(cpu_val, gpu_val, 0.00001);

  }

  /**
   * Test that GPU local MI calculations match the CPU version.
   *
   * Only runs if CUDA code has been precompiled.
   *
   * @throws Exception if something goes wrong
   */
  public void testGPULocalCondMutualInfo() throws Exception {

    ConditionalMutualInfoCalculatorMultiVariateKraskov cmiCalc = 
      new ConditionalMutualInfoCalculatorMultiVariateKraskov1();

    boolean gpuLoaded = true;
    try {
      cmiCalc.ensureCudaLibraryLoaded();
    } catch (Throwable e) {
      gpuLoaded = false;
    }

    // This will effectively ignore the test if GPU library not loaded properly
    if (!gpuLoaded) {
      return;
    }

    // Now starts the actual test. Set up variables and properties
    cmiCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    ArrayFileReader afr;
    double[][] data;
    double[] cpu_vals, gpu_vals;

    // Test for low-dimensional data
		afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		data = afr.getDouble2DMatrix();

    cmiCalc.setProperty("USE_GPU", "false");
    cmiCalc.initialise(1,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    cpu_vals = cmiCalc.computeLocalOfPreviousObservations();

    cmiCalc.setProperty("USE_GPU", "true");
    cmiCalc.initialise(1,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    gpu_vals = cmiCalc.computeLocalOfPreviousObservations();

    for (int i = 0; i < data.length; i++) {
      assertEquals(cpu_vals[i], gpu_vals[i], 0.00001);
    }

    // Test for slightly higher-dimensional data
		afr = new ArrayFileReader("demos/data/4ColsPairedNoisyDependence-1.txt");
		data = afr.getDouble2DMatrix();

    cmiCalc.setProperty("USE_GPU", "false");
    cmiCalc.initialise(2,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumns(data, new int[]{0,1}),
                            MatrixUtils.selectColumns(data, new int[]{2}),
                            MatrixUtils.selectColumns(data, new int[]{3}));
    cpu_vals = cmiCalc.computeLocalOfPreviousObservations();

    cmiCalc.setProperty("USE_GPU", "true");
    cmiCalc.initialise(2,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumns(data, new int[]{0,1}),
                            MatrixUtils.selectColumns(data, new int[]{2}),
                            MatrixUtils.selectColumns(data, new int[]{3}));
    gpu_vals = cmiCalc.computeLocalOfPreviousObservations();

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

    ConditionalMutualInfoCalculatorMultiVariateKraskov cmiCalc = 
      new ConditionalMutualInfoCalculatorMultiVariateKraskov1();
    // cmiCalc.setDebug(true);

    boolean gpuLoaded = true;
    try {
      cmiCalc.ensureCudaLibraryLoaded();
    } catch (Throwable e) {
      gpuLoaded = false;
    }

    // This will effectively ignore the test if GPU library not loaded properly
    if (!gpuLoaded) {
      return;
    }

    // Now starts the actual test. Set up variables and properties
    cmiCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    ArrayFileReader afr;
    EmpiricalMeasurementDistribution cpu_dist, gpu_dist;

    // Test for low-dimensional data
		afr = new ArrayFileReader("demos/data/4ColsPairedNoisyDependence-1.txt");
		double [][] data = afr.getDouble2DMatrix();

    RandomGenerator rng = new RandomGenerator();
    int[] perm = rng.generateRandomPerturbations(data.length, 1)[0];
    int[][] newOrderings = new int[][] { MatrixUtils.range(0, data.length-1),
                                         perm};

    cmiCalc.setProperty("USE_GPU", "false");
    cmiCalc.initialise(1,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    cpu_dist = cmiCalc.computeSignificance(newOrderings);

    cmiCalc.setProperty("USE_GPU", "true");
    cmiCalc.initialise(1,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    gpu_dist = cmiCalc.computeSignificance(newOrderings);

    System.out.printf("GPU CMI calculation on 4ColsPairedNoisyDependence data. CPU got (%f,%f,%f)\n", cpu_dist.actualValue, cpu_dist.distribution[0], cpu_dist.distribution[1]);
    System.out.printf("GPU CMI calculation on 4ColsPairedNoisyDependence data. GPU got (%f,%f,%f)\n", gpu_dist.actualValue, gpu_dist.distribution[0], gpu_dist.distribution[1]);

    assertEquals(cpu_dist.actualValue,     cpu_dist.distribution[0], 0.00001);
    assertEquals(gpu_dist.actualValue,     gpu_dist.distribution[0], 0.00001);
    assertEquals(cpu_dist.actualValue,     gpu_dist.distribution[0], 0.00001);
    assertEquals(gpu_dist.actualValue,     cpu_dist.distribution[0], 0.00001);

    // Smoke test for reordering with different variables
    cmiCalc.setProperty("USE_GPU", "true");
    cmiCalc.initialise(1,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    gpu_dist = cmiCalc.computeSignificance(1, 10);

    cmiCalc.initialise(1,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    gpu_dist = cmiCalc.computeSignificance(2, 10);

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

    ConditionalMutualInfoCalculatorMultiVariateKraskov cmiCalc = 
      new ConditionalMutualInfoCalculatorMultiVariateKraskov1();

    boolean gpuLoaded = true;
    try {
      cmiCalc.ensureCudaLibraryLoaded();
    } catch (Throwable e) {
      gpuLoaded = false;
    }

    // This will effectively ignore the test if GPU library not loaded properly
    if (!gpuLoaded) {
      return;
    }

    // Now starts the actual test. Set up variables and properties
    cmiCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");

    // Generate correlated Gaussian data
    RandomGenerator rng = new RandomGenerator();
    double[][] data = rng.generateCovariantGaussians(10000, 3, new double[]{0,0,0}, new double[][]{{0, 0.1, 0.2},{0,0,0.3},{0,0,0}});

    // cmiCalc.setDebug(true);
    cmiCalc.setProperty("USE_GPU", "true");
    cmiCalc.initialise(1, 1, 1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    EmpiricalMeasurementDistribution dist =
      cmiCalc.computeSignificance(10);

    System.out.printf("GPU CMI random surrogates test: actual: %f, mean: %f, std: %f\n", dist.actualValue, dist.getMeanOfDistribution(), dist.getStdOfDistribution());
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

    ConditionalMutualInfoCalculatorMultiVariateKraskov cmiCalc = 
      new ConditionalMutualInfoCalculatorMultiVariateKraskov1();
    // cmiCalc.setDebug(true);

    boolean gpuLoaded = true;
    try {
      cmiCalc.ensureCudaLibraryLoaded();
    } catch (Throwable e) {
      gpuLoaded = false;
    }

    // This will effectively ignore the test if GPU library not loaded properly
    if (!gpuLoaded) {
      return;
    }

    // Now starts the actual test. Set up variables and properties
    cmiCalc.setProperty("NOISE_LEVEL_TO_ADD", "0");
    cmiCalc.setProperty("DYN_CORR_EXCL", "2");
    ArrayFileReader afr;
    double cpu_val, gpu_val;

    // Test for low-dimensional data
		afr = new ArrayFileReader("demos/data/4ColsPairedNoisyDependence-1.txt");
		double [][] data = afr.getDouble2DMatrix();

    cmiCalc.setProperty("USE_GPU", "false");
    cmiCalc.initialise(1,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    cpu_val = cmiCalc.computeAverageLocalOfObservations();

    cmiCalc.setProperty("USE_GPU", "true");
    cmiCalc.initialise(1,1,1);
    cmiCalc.setObservations(MatrixUtils.selectColumn(data, 0),
                            MatrixUtils.selectColumn(data, 1),
                            MatrixUtils.selectColumn(data, 2));
    gpu_val = cmiCalc.computeAverageLocalOfObservations();

    System.out.printf("GPU CMI calculation with theiler window. CPU got (%f)\n", cpu_val);
    System.out.printf("GPU CMI calculation with theiler window. GPU got (%f)\n", gpu_val);

    assertEquals(cpu_val, cpu_val, 0.00001);

    return;
  }

}

