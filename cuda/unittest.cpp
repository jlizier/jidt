#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "lest.hpp"

#include "gpuKnnLibrary.h"
#include "gpuMILibrary.h"
#include "gpuCMILibrary.h"
#include "digamma.h"

using lest::approx;

const lest::test specification[] = {

CASE("Basic kNN test in 1D")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 4;
  int dims = 1;
  int k = 1;
  int nchunks = 1;
  int indexes[4];
  float distances[4];
  float pointset[4] = {-1, -1.2, 1, 1.1};

  int err = cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, N, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(indexes[0] == 1);
  EXPECT(indexes[1] == 0);
  EXPECT(indexes[2] == 3);
  EXPECT(indexes[3] == 2);

  EXPECT(distances[0] == approx(0.2));
  EXPECT(distances[1] == approx(0.2));
  EXPECT(distances[2] == approx(0.1));
  EXPECT(distances[3] == approx(0.1));

},

CASE("Basic kNN test in 2D")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 4;
  int dims = 2;
  int k = 1;
  int nchunks = 1;
  int indexes[4];
  float distances[4];
  float pointset[8] = {-1, 0.5, 1.1, 2,
                       -1, 0.5, 1.1, 2};

  int err = cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, N, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(indexes[0] == 1);
  EXPECT(indexes[1] == 2);
  EXPECT(indexes[2] == 1);
  EXPECT(indexes[3] == 2);

  EXPECT(distances[0] == approx(1.5));
  EXPECT(distances[1] == approx(0.6));
  EXPECT(distances[2] == approx(0.6));
  EXPECT(distances[3] == approx(0.9));

},

CASE("Basic kNN test with L2 norm")
{
  int thelier = 0;
  int useMaxNorm = 0;
  int N = 4;
  int dims = 2;
  int k = 1;
  int nchunks = 1;
  int indexes[4];
  float distances[4];

  // Points:  X     Y
  //         -1    -1
  //        0.5   0.5
  //        1.1   1.1
  //          2     2

  float pointset[8] = {-1, 0.5, 1.1, 2,
                       -1, 0.5, 1.1, 2};

  int err = cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, N, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(indexes[0] == 1);
  EXPECT(indexes[1] == 2);
  EXPECT(indexes[2] == 1);
  EXPECT(indexes[3] == 2);

  EXPECT(distances[0] == approx(4.5));
  EXPECT(distances[1] == approx(0.72));
  EXPECT(distances[2] == approx(0.72));
  EXPECT(distances[3] == approx(1.62));

},

CASE("Test kNN test with longer sequences")
{
  int thelier = 0;
  int useMaxNorm = 0;
  int N = 10;
  int dims = 2;
  int k = 1;
  int nchunks = 1;
  int indexes[10];
  float distances[10];
  // This is the same sequence as in the previous test case, padded with a
  // bunch of points very far away.
  float pointset[20] = {-1, 0.5, 1.1, 2, 10, 11, 10.5, -100, -50, 666,
                        -1, 0.5, 1.1, 2, 98, -9, -200, 45.3, -53, 0.1};

  int err = cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, N, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(indexes[0] == 1);
  EXPECT(indexes[1] == 2);
  EXPECT(indexes[2] == 1);
  EXPECT(indexes[3] == 2);

  EXPECT(distances[0] == approx(4.5));
  EXPECT(distances[1] == approx(0.72));
  EXPECT(distances[2] == approx(0.72));
  EXPECT(distances[3] == approx(1.62));

},

CASE("Smoke kNN test with big random dataset")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 1000;
  int dims = 5;
  int k = 4;
  int nchunks = 1;
  int *indexes = (int *) malloc(N * k * sizeof(int));
  float *distances = (float *) malloc(N * k * sizeof(float));
  float *pointset = (float *) malloc(N * dims * sizeof(float));

  int err = cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, N, useMaxNorm);

  free(indexes); free(distances); free(pointset);

  EXPECT(err == 1);

},

CASE("Test kNN with two trials in the same call")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 4;
  int dims = 1;
  int k = 1;
  int nchunks = 2;
  int signalLength = N*nchunks;
  int indexes[8];
  float distances[8];
  float pointset[8] = {5,   6, -5,  -7,
                      50, -50, 60, -70};

  int err = cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, signalLength, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(indexes[0]  == 1);
  EXPECT(indexes[1]  == 0);
  EXPECT(indexes[2]  == 3);
  EXPECT(indexes[3]  == 2);
  EXPECT(indexes[4]  == 2);
  EXPECT(indexes[5]  == 3);
  EXPECT(indexes[6]  == 0);
  EXPECT(indexes[7]  == 1);

  EXPECT(distances[0] == approx(1));
  EXPECT(distances[1] == approx(1));
  EXPECT(distances[2] == approx(2));
  EXPECT(distances[3] == approx(2));
  EXPECT(distances[4] == approx(10));
  EXPECT(distances[5] == approx(20));
  EXPECT(distances[6] == approx(10));
  EXPECT(distances[7] == approx(20));

},

CASE("Test kNN with three trials in the same call")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 4;
  int dims = 1;
  int k = 1;
  int nchunks = 3;
  int signalLength = N*nchunks;
  int indexes[12];
  float distances[12];
  float pointset[12] = {  5,    6,  -5,   -7,
                         50,  -50,  60,  -70,
                        500, -500, 600, -700};

  int err = cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, signalLength, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(indexes[0]   == 1);
  EXPECT(indexes[1]   == 0);
  EXPECT(indexes[2]   == 3);
  EXPECT(indexes[3]   == 2);
  EXPECT(indexes[4]   == 2);
  EXPECT(indexes[5]   == 3);
  EXPECT(indexes[6]   == 0);
  EXPECT(indexes[7]   == 1);
  EXPECT(indexes[8]   == 2);
  EXPECT(indexes[9]   == 3);
  EXPECT(indexes[10]  == 0);
  EXPECT(indexes[11]  == 1);

  EXPECT(distances[0]  == approx(1));
  EXPECT(distances[1]  == approx(1));
  EXPECT(distances[2]  == approx(2));
  EXPECT(distances[3]  == approx(2));
  EXPECT(distances[4]  == approx(10));
  EXPECT(distances[5]  == approx(20));
  EXPECT(distances[6]  == approx(10));
  EXPECT(distances[7]  == approx(20));
  EXPECT(distances[8]  == approx(100));
  EXPECT(distances[9]  == approx(200));
  EXPECT(distances[10] == approx(100));
  EXPECT(distances[11] == approx(200));

},

CASE("Test kNN with two trials of 2D data in the same call")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 4;
  int dims = 2;
  int k = 1;
  int nchunks = 2;
  int signalLength = N*nchunks;
  int indexes[8];
  float distances[8];

  // Points:       X    Y                   y
  //               1    1                   |  o o
  //             1.1    1                   |
  //              -1   -1               ----+----x
  //            -1.2   -1                   |
  //                                  o  o  |

  float pointset[16] = {1, 1.1, -1, -1.2, 1, 1.1, -1, -1.2,
                        1,   1, -1,   -1, 1,   1, -1,   -1};

  int err = cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, signalLength, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(indexes[0]   == 1);
  EXPECT(indexes[1]   == 0);
  EXPECT(indexes[2]   == 3);
  EXPECT(indexes[3]   == 2);
  EXPECT(indexes[4]   == 1);
  EXPECT(indexes[5]   == 0);
  EXPECT(indexes[6]   == 3);
  EXPECT(indexes[7]   == 2);

  EXPECT(distances[0]  == approx(0.1));
  EXPECT(distances[1]  == approx(0.1));
  EXPECT(distances[2]  == approx(0.2));
  EXPECT(distances[3]  == approx(0.2));
  EXPECT(distances[4]  == approx(0.1));
  EXPECT(distances[5]  == approx(0.1));
  EXPECT(distances[6]  == approx(0.2));
  EXPECT(distances[7]  == approx(0.2));
},

CASE("Test kNN with two trials of data with odd dimension")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 4;
  int dims = 3;
  int k = 1;
  int nchunks = 2;
  int signalLength = N*nchunks;
  int indexes[8];
  float distances[8];

  // Points:       X    Y      Z             y
  //               1    1   1.02            |  o o
  //             1.1    1   1.03            |
  //              -1   -1  -1.04        ----+----x
  //            -1.2   -1  -1.05            |
  //                                  o  o  |

  float pointset[24] = {1,  1.1,   -1, -1.2,    1,  1.1,   -1, -1.2,
                        1,    1,   -1,   -1,    1,    1,   -1,   -1,
                     1.02, 1.03, 1.04, 1.05, 1.02, 1.03, 1.04, 1.05};

  int err = cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, signalLength, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(indexes[0]   == 1);
  EXPECT(indexes[1]   == 0);
  EXPECT(indexes[2]   == 3);
  EXPECT(indexes[3]   == 2);
  EXPECT(indexes[4]   == 1);
  EXPECT(indexes[5]   == 0);
  EXPECT(indexes[6]   == 3);
  EXPECT(indexes[7]   == 2);

  EXPECT(distances[0]  == approx(0.1));
  EXPECT(distances[1]  == approx(0.1));
  EXPECT(distances[2]  == approx(0.2));
  EXPECT(distances[3]  == approx(0.2));
  EXPECT(distances[4]  == approx(0.1));
  EXPECT(distances[5]  == approx(0.1));
  EXPECT(distances[6]  == approx(0.2));
  EXPECT(distances[7]  == approx(0.2));
},

CASE("Smoke test for kNN with custom GPU")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 4;
  int dims = 2;
  int k = 1;
  int nchunks = 1;
  int indexes[4];
  float distances[4];
  float pointset[8] = {-1, 0.5, 1.1, 2,
                       -1, 0.5, 1.1, 2};

  int err = cudaFindKnnSetGPU(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, N, useMaxNorm, 0);

  EXPECT(err == 1);

  EXPECT(indexes[0] == 1);
  EXPECT(indexes[1] == 2);
  EXPECT(indexes[2] == 1);
  EXPECT(indexes[3] == 2);

  EXPECT(distances[0] == approx(1.5));
  EXPECT(distances[1] == approx(0.6));
  EXPECT(distances[2] == approx(0.6));
  EXPECT(distances[3] == approx(0.9));

},

CASE("Basic RS test")
{
  int thelier = 0;
  int nchunks = 1;
  int dim = 1;
  int N = 5;
  int useMaxNorm = 1;
  int npoints[5];
  float pointset[] = {0, 2, 2.1, 3.2, 3.5};
  float radii[]    = {1, 1,   1,   1,   1};

  int err = cudaFindRSAll(npoints, pointset, pointset, radii, thelier, nchunks,
              dim, N, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(npoints[0] == 0);
  EXPECT(npoints[1] == 1);
  EXPECT(npoints[2] == 1);
  EXPECT(npoints[3] == 1);
  EXPECT(npoints[4] == 1);

},

CASE("Test RS with different radii")
{
  int thelier = 0;
  int nchunks = 1;
  int dim = 1;
  int N = 5;
  int useMaxNorm = 1;
  int npoints[5];
  float pointset[] = {   0, 2,  2.1, 3.2, 3.5};
  float radii[]    = {2.05, 6, 0.01,   1,   1};

  int err = cudaFindRSAll(npoints, pointset, pointset, radii, thelier, nchunks,
              dim, N, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(npoints[0] == 1);
  EXPECT(npoints[1] == 4);
  EXPECT(npoints[2] == 0);
  EXPECT(npoints[3] == 1);
  EXPECT(npoints[4] == 1);

},

CASE("Test RS with multiple trials in the same call")
{
  int thelier = 0;
  int nchunks = 2;
  int dim = 1;
  int N = 10;
  int useMaxNorm = 1;
  int npoints[10];
  float pointset[] = {0, 2, 2.1, 3.2, 3.5, 0, 2, 2.1, 3.2, 3.5};
  float radii[]    = {1, 1,   1,   1,   1, 1, 1,   1,   1,   1};

  int err = cudaFindRSAll(npoints, pointset, pointset, radii, thelier, nchunks,
              dim, N, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(npoints[0] == 0);
  EXPECT(npoints[1] == 1);
  EXPECT(npoints[2] == 1);
  EXPECT(npoints[3] == 1);
  EXPECT(npoints[4] == 1);
  EXPECT(npoints[5] == 0);
  EXPECT(npoints[6] == 1);
  EXPECT(npoints[7] == 1);
  EXPECT(npoints[8] == 1);
  EXPECT(npoints[9] == 1);

},

CASE("Test RS with multiple trials and different radii")
{
  int thelier = 0;
  int nchunks = 2;
  int dim = 1;
  int N = 10;
  int useMaxNorm = 1;
  int npoints[10];
  float pointset[] = {   0, 2,  2.1, 3.2, 3.5,    0, 2,  2.1, 3.2, 3.5};
  float radii[]    = {2.05, 6, 0.01,   1,   1, 2.05, 6, 0.01,   1,   1};

  int err = cudaFindRSAll(npoints, pointset, pointset, radii, thelier, nchunks,
              dim, N, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(npoints[0] == 1);
  EXPECT(npoints[1] == 4);
  EXPECT(npoints[2] == 0);
  EXPECT(npoints[3] == 1);
  EXPECT(npoints[4] == 1);
  EXPECT(npoints[5] == 1);
  EXPECT(npoints[6] == 4);
  EXPECT(npoints[7] == 0);
  EXPECT(npoints[8] == 1);
  EXPECT(npoints[9] == 1);

},

CASE("Smoke test for RS with custom GPU")
{
  int thelier = 0;
  int nchunks = 1;
  int dim = 1;
  int N = 5;
  int useMaxNorm = 1;
  int npoints[5];
  float pointset[] = {0, 2, 2.1, 3.2, 3.5};
  float radii[]    = {1, 1,   1,   1,   1};

  int err = cudaFindRSAllSetGPU(npoints, pointset, pointset, radii, thelier, nchunks,
              dim, N, useMaxNorm, 0);

  EXPECT(err == 1);

  EXPECT(npoints[0] == 0);
  EXPECT(npoints[1] == 1);
  EXPECT(npoints[2] == 1);
  EXPECT(npoints[3] == 1);
  EXPECT(npoints[4] == 1);

},

CASE("Smoke test of random permutation function")
{
  int N = 30;
  int perm[30];
  for (int i = 0; i < N; i++) {
    perm[i] = i;
  }
  randperm(perm, N);

},

CASE("Statistical test for unbiasedness of random permutation function")
{
  // Initialise permutation
  int N = 3;
  int perm[3];
  for (int i = 0; i < N; i++) {
    perm[i] = i;
  }

  // Sample a bunch of times and store results in a histogram
  int nb_samples = 150000;
  double hist[6] = {0};
  for (int i = 0; i < nb_samples; i++) {
    randperm(perm, N);
    if        ((perm[0] == 0) && (perm[1] == 1) && (perm[2] == 2)){
      hist[0]++;
    } else if ((perm[0] == 0) && (perm[1] == 2) && (perm[2] == 1)){
      hist[1]++;
    } else if ((perm[0] == 1) && (perm[1] == 0) && (perm[2] == 2)){
      hist[2]++;
    } else if ((perm[0] == 1) && (perm[1] == 2) && (perm[2] == 0)){
      hist[3]++;
    } else if ((perm[0] == 2) && (perm[1] == 0) && (perm[2] == 1)){
      hist[4]++;
    } else if ((perm[0] == 2) && (perm[1] == 1) && (perm[2] == 0)){
      hist[5]++;
    } else {
      std::cout << "This permutation is certainly wrong\n";
      EXPECT(0);
    }
  }

  // Check that all perturbations happen more or less equally often
  for (int i = 0; i < 6; i++) {
    EXPECT(hist[i]/((double) nb_samples) == approx(1.0/6.0).epsilon(0.02));
  }

},

CASE("Smoke test of full MI function")
{
  int N = 10;
  int dimx = 1;
  int dimy = 1;
  float source[10] = {0.4, 1, -4, 1, 1, 0.2, 98, 12, 1.2, 1.3};
  float dest[10]   = {-3, 1, 3, -2, 2.1, 8.5, 4.2, 100, 12, 0};
  int k = 2;
  int thelier = 0;
  int nb_surrogates = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;
  float result[3];

  jidt_error_t err = MIKraskov_C(N, source, dimx, dest, dimy, k, thelier,
      nb_surrogates, returnLocals, useMaxNorm, isAlgorithm1, result);

  EXPECT(err == JIDT_SUCCESS);

},

CASE("Test correct pointset arrangement without reorderings")
{
  int N = 10;
  int dimx = 1;
  int dimy = 1;
  int k = 2;
  int thelier = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;

  float source[10] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3};
  float dest[10]   = { -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  12,   0};

  float pointset[20] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3,
                         -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  12,   0};

  float result1[3];
  float result2[3];
  jidt_error_t err;

  err = MIKraskov_C(N, source, dimx, dest, dimy, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1);
  EXPECT(err == JIDT_SUCCESS);

  err = MIKraskovByPointsetChunks(N, source, dimx, dest, dimy, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result2, pointset);
  EXPECT(err == JIDT_SUCCESS);

  EXPECT(result1[0] == approx(result2[0]));
  EXPECT(result1[1] == approx(result2[1]));
  EXPECT(result1[2] == approx(result2[2]));
},

CASE("Test correct pointset arrangement without reorderings in CMI")
{
  int N = 10;
  int dimx = 1;
  int dimy = 1;
  int dimz = 1;
  int k = 2;
  int thelier = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;

  float source[10] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3};
  float dest[10]   = { -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  12,   0};
  float cond[10]   = { -1, 4,  3, -8, 0.3, 2.1, 3.2, 111,  32,   7};

  float pointset[30] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3,
                         -1, 4,  3, -8, 0.3, 2.1, 3.2, 111,  32,   7,
                         -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  12,   0};

  float result1[6];
  float result2[6];
  jidt_error_t err;

  err = CMIKraskov_C(N, source, dimx, dest, dimy, cond, dimz, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1, 1);
  EXPECT(err == JIDT_SUCCESS);

  err = CMIKraskovByPointsetChunks(N, source, dimx, dest, dimy, cond, dimz, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result2, pointset);
  EXPECT(err == JIDT_SUCCESS);

  EXPECT(result1[0] == approx(result2[0]));
  EXPECT(result1[1] == approx(result2[1]));
  EXPECT(result1[2] == approx(result2[2]));
  EXPECT(result1[3] == approx(result2[3]));
  EXPECT(result1[4] == approx(result2[4]));
  EXPECT(result1[5] == approx(result2[5]));
},

CASE("Test correct pointset arrangement in more than one dimension")
{
  int N = 5;
  int dimx = 2;
  int dimy = 2;
  int k = 2;
  int thelier = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;

  // Source points:  X      Y
  //               0.4    0.2
  //                 1     98
  //                -4     12
  //                 1    1.2
  //                 1    1.3
  //
  // Dest points:    X      Y
  //                -3    8.5
  //                 1    4.2
  //                 3    100
  //                -2     12
  //                2.1     0

  float source[10] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3};
  float dest[10]   = { -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  12,   0};

  float pointset[20] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3,
                         -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  12,   0};

  float result1[3];
  float result2[3];
  jidt_error_t err;

  err = MIKraskov_C(N, source, dimx, dest, dimy, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1);
  EXPECT(err == JIDT_SUCCESS);

  err = MIKraskovByPointsetChunks(N, source, dimx, dest, dimy, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result2, pointset);
  EXPECT(err == JIDT_SUCCESS);

  EXPECT(result1[0] == approx(result2[0]));
  EXPECT(result1[1] == approx(result2[1]));
  EXPECT(result1[2] == approx(result2[2]));
},

CASE("Test correct pointset arrangement in more than one dimension for CMI")
{
  int N = 5;
  int dimx = 2;
  int dimy = 2;
  int dimz = 2;
  int k = 2;
  int thelier = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;

  // Source points:  X      Y
  //               0.4    0.2
  //                 1     98
  //                -4     12
  //                 1    1.2
  //                 1    1.3
  //
  // Dest points:    X      Y
  //                -3    8.5
  //                 1    4.2
  //                 3    100
  //                -2     12
  //               2.1      0
  //
  // Cond points:    X      Y
  //                -1    2.1
  //                 4    3.2
  //                 3    111
  //                -8     32
  //               0.3      7

  float source[10] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3};
  float dest[10]   = { -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  12,   0};
  float cond[10]   = { -1, 4,  3, -8, 0.3, 2.1, 3.2, 111,  32,   7};

  float pointset[30] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3,
                         -1, 4,  3, -8, 0.3, 2.1, 3.2, 111,  32,   7,
                         -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  12,   0};

  float result1[6];
  float result2[6];
  jidt_error_t err;

  err = CMIKraskov_C(N, source, dimx, dest, dimy, cond, dimz, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1, 1);
  EXPECT(err == JIDT_SUCCESS);

  err = CMIKraskovByPointsetChunks(N, source, dimx, dest, dimy, cond, dimz, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result2, pointset);
  EXPECT(err == JIDT_SUCCESS);

  EXPECT(result1[0] == approx(result2[0]));
  EXPECT(result1[1] == approx(result2[1]));
  EXPECT(result1[2] == approx(result2[2]));
  EXPECT(result1[3] == approx(result2[3]));
  EXPECT(result1[4] == approx(result2[4]));
  EXPECT(result1[5] == approx(result2[5]));
},

CASE("Test that same sample in repeated chunks gives same result")
{
  int N = 5;
  int dimx = 1;
  int dimy = 1;

  // Sample source and dest data
  float source[5] = {0.4, 1, -4,  1,   1};
  float dest[5]   = { -3, 1,  3, -2, 2.1};

  // Pointset with source and dest repeated twice
  float double_pointset[20] = {0.4, 1, -4,  1,   1, 0.4, 1, -4,  1,   1,
                                -3, 1,  3, -2, 2.1,  -3, 1,  3, -2, 2.1};

  int k = 2;
  int thelier = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;
  float result1[3];
  float result2[2];
  jidt_error_t err;

  err = MIKraskov_C(N, source, dimx, dest, dimy, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1);

  err = MIKraskovByPointsetChunks(N*2, source, dimx, dest, dimy, k, thelier,
      2, returnLocals, useMaxNorm, isAlgorithm1, result2, double_pointset);

  float MI1 = cpuDigamma(k) + cpuDigamma(N) - result1[0]/((double) N);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(result2[0] == approx(MI1));
  EXPECT(result2[0] == approx(result2[1]));

},

CASE("Test that same sample in repeated chunks gives same result in CMI")
{
  int N = 5;
  int dimx = 1;
  int dimy = 1;
  int dimz = 1;

  // Sample source and dest data
  float source[5] = {0.4, 1, -4,  1,   1};
  float dest[5]   = { -3, 1,  3, -2, 2.1};
  float cond[5]   = { -1, 4,  3, -8, 0.3};

  // Pointset with source and dest repeated twice
  float double_pointset[30] = {0.4, 1, -4,  1,   1, 0.4, 1, -4,  1,   1,
                                -1, 4,  3, -8, 0.3,  -1, 4,  3, -8, 0.3,
                                -3, 1,  3, -2, 2.1,  -3, 1,  3, -2, 2.1};

  int k = 2;
  int thelier = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;
  float result1[3];
  float result2[2];
  jidt_error_t err;

  err = CMIKraskov_C(N, source, dimx, dest, dimy, cond, dimz, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1, 1);

  err = CMIKraskovByPointsetChunks(N*2, source, dimx, dest, dimy, cond, dimz, k, thelier,
      2, returnLocals, useMaxNorm, isAlgorithm1, result2, double_pointset);

  float CMI1 = cpuDigamma(k) + result1[0]/((double) N);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(result2[0] == approx(CMI1));
  EXPECT(result2[0] == approx(result2[1]));

},

CASE("Test that same sample of 2D data in repeated chunks gives same result")
{
  int N = 5;
  int dimx = 2;
  int dimy = 2;

  // Source points:  X      Y
  //               0.4    0.2
  //                 1     98
  //                -4     12
  //                 1    1.2
  //                 1    1.3
  //
  // Dest points:    X      Y
  //                -3    8.5
  //                 1    4.2
  //                 3    100
  //                -2     12
  //                2.1     0

  // Sample source and dest data
  float source[10] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3};
  float dest[10]   = { -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  12,   0};

  // Pointset with source and dest repeated twice
  float double_pointset[40] = {0.4,   1,  -4,   1,   1, 0.4,   1,  -4,   1,   1,
                               0.2,  98,  12, 1.2, 1.3, 0.2,  98,  12, 1.2, 1.3,
                                -3,   1,   3,  -2, 2.1,  -3,   1,   3,  -2, 2.1,
                               8.5, 4.2, 100,  12,   0, 8.5, 4.2, 100,  12,   0};

  int k = 2;
  int thelier = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;
  float result1[3];
  float result2[2];
  jidt_error_t err;

  err = MIKraskov_C(N, source, dimx, dest, dimy, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1);

  err = MIKraskovByPointsetChunks(N*2, source, dimx, dest, dimy, k, thelier,
      2, returnLocals, useMaxNorm, isAlgorithm1, result2, double_pointset);

  float MI1 = cpuDigamma(k) + cpuDigamma(N) - result1[0]/((double) N);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(result2[0] == approx(MI1));
  EXPECT(result2[0] == approx(result2[1]));

},

CASE("Test that sample and identity reordering have same MI")
{
  int N = 5;
  int dimx = 1;
  int dimy = 1;
  float source[5] = {0.4, 1, -4,  1,   1};
  float dest[5]   = { -3, 1,  3, -2, 2.1};
  int k = 2;
  int thelier = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;
  float result[2];
  int reorderingsGiven = 1;
  int order[5] = {0, 1, 2, 3, 4};
  int *order_p = order;
  int **reorderings = &order_p;
  jidt_error_t err;

  err = MIKraskovWithReorderings(N, source, dimx, dest, dimy, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result, reorderingsGiven, reorderings);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(result[0] == approx(result[1]));

},

CASE("Test that sample and identity reordering have same CMI")
{
  int N = 5;
  int dimx = 1;
  int dimy = 1;
  int dimz = 1;
  float source[5] = {0.4, 1, -4,  1,   1};
  float dest[5]   = { -3, 1,  3, -2, 2.1};
  float cond[5]   = { -1, 4,  3, -8, 0.3};
  int k = 2;
  int thelier = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;
  float result[2];
  int reorderingsGiven = 1;
  int order[5] = {0, 1, 2, 3, 4};
  int *order_p = order;
  int **reorderings = &order_p;
  jidt_error_t err;

  err = CMIKraskovWithReorderings(N, source, dimx, dest, dimy, cond, dimz, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result, reorderingsGiven, reorderings, 1);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(result[0] == approx(result[1]));

},

CASE("Test identity reordering with more than one dimension")
{
  int N = 5;
  int dimx = 3;
  int dimy = 2;

  // Source points:  X      Y    Z
  //               0.4    0.2    0
  //                 1     98   13
  //                -4     12    7
  //                 1    1.2   -1
  //                 1    1.3    0
  //
  // Dest points:    X      Y
  //                -3    8.5
  //                 1    4.2
  //                 3    100
  //                -2     13
  //                2.1     0

  float source[15] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3, 0, 13, 7, -1, 0};
  float dest[10]   = { -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  13,   0};
  int k = 1;
  int thelier = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;
  float result1[3], result2[2];
  int reorderingsGiven = 1;
  int order[5] = {0, 1, 2, 3, 4};
  int *order_p = order;
  int **reorderings = &order_p;
  jidt_error_t err;

  err = MIKraskovWithReorderings(N, source, dimx, dest, dimy, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1, reorderingsGiven, reorderings);

  EXPECT(err == JIDT_SUCCESS);

  err = MIKraskovWithReorderings(N, source, dimx, dest, dimy, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result2, reorderingsGiven, reorderings);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(result2[0] == approx(result2[1]));

},

CASE("Test identity reordering with larger dataset")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 10;
  int dimx = 3;
  int dimy = 2;
  int k = 4;
  int isAlgorithm1 = 1;
  int returnLocals = 0;

  float *source = (float *) malloc(N * dimx * sizeof(float));
  for (int i = 0; i < N*dimx; i++) { source[i] = rand()/((float) RAND_MAX); };

  float *dest   = (float *) malloc(N * dimy * sizeof(float));
  for (int i = 0; i < N*dimy; i++) { dest[i] = rand()/((float) RAND_MAX); };

  int reorderingsGiven = 1;
  int *order = (int *) malloc(N * sizeof(int));
  for (int i = 0; i < N; i++) { order[i] = i; };
  int **reorderings = &order;

  float result[2];

  jidt_error_t err;
  err = MIKraskovWithReorderings(N, source, dimx, dest, dimy, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result, reorderingsGiven, reorderings);

  free(source); free(dest); free(order);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(result[0] == approx(result[1]));

},

CASE("Test non-identity reordering in 2D")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 5;
  int dimx = 2;
  int dimy = 2;
  int k = 2;
  int isAlgorithm1 = 1;
  int returnLocals = 0;

  // Source points:  X      Y    Z
  //               0.4    0.2    0
  //                 1     98   13
  //                -4     12    7
  //                 1    1.2   -1
  //                 1    1.3    0
  //
  // Dest points:    X      Y
  //                -3    8.5
  //                 1    4.2
  //                 3    100
  //                -2     13
  //                2.1     0

  float source[15] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3, 0, 13, 7, -1, 0};
  float dest[10]   = { -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  13,   0};


  int reorderingsGiven = 1;
  int order[5] = {4, 2, 3, 0, 1};
  int *order_p = order;
  int **reorderings = &order_p;

  float result[2];

  jidt_error_t err;
  err = MIKraskovWithReorderings(N, source, dimx, dest, dimy, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result, reorderingsGiven, reorderings);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(result[0] != approx(result[1]));

},

CASE("Test random surrogates in 2D")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 20;
  int dimx = 2;
  int dimy = 2;
  int k = 2;
  int isAlgorithm1 = 1;
  int returnLocals = 0;


  float *source = (float *) malloc(N * dimx * sizeof(float));
  for (int i = 0; i < N*dimx; i++) { source[i] = rand()/((float) RAND_MAX); };

  float *dest   = (float *) malloc(N * dimy * sizeof(float));
  for (int i = 0; i < N*dimy; i++) { dest[i] = rand()/((float) RAND_MAX); };

  float result[3];

  jidt_error_t err;
  err = MIKraskov_C(N, source, dimx, dest, dimy, k, thelier,
      2, returnLocals, useMaxNorm, isAlgorithm1, result);

  free(source); free(dest);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(result[0] != result[1]);
  EXPECT(result[0] != result[2]);
  EXPECT(result[1] != result[2]);

},

CASE("Test that the first result of calculation with surrogates is the same as without")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 10;
  int dimx = 1;
  int dimy = 1;
  int k = 2;
  int isAlgorithm1 = 1;
  int returnLocals = 0;

  // Source points:  X      Y    Z
  //               0.4    0.2    0
  //                 1     98   13
  //                -4     12    7
  //                 1    1.2   -1
  //                 1    1.3    0
  //
  // Dest points:    X      Y
  //                -3    8.5
  //                 1    4.2
  //                 3    100
  //                -2     13
  //                2.1     0

  float source[15] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3, 0, 13, 7, -1, 0};
  float dest[10]   = { -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  13,   0};

  float result1[3];
  float result2[2];

  jidt_error_t err;
  err = MIKraskov_C(N, source, dimx, dest, dimy, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1);

  EXPECT(err == JIDT_SUCCESS);

  float MI1 = cpuDigamma(k) + cpuDigamma(N) - result1[0]/((double) N);

  err = MIKraskov_C(N, source, dimx, dest, dimy, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result2);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(MI1 == approx(result2[0]));
  EXPECT(result2[0] != result2[1]);
},

CASE("Test that the first result of calculation with surrogates is the same as without in CMI")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 10;
  int dimx = 1;
  int dimy = 1;
  int dimz = 1;
  int k = 2;
  int isAlgorithm1 = 1;
  int returnLocals = 0;

  // Source points:  X      Y    Z
  //               0.4    0.2    0
  //                 1     98   13
  //                -4     12    7
  //                 1    1.2   -1
  //                 1    1.3    0
  //
  // Dest points:    X      Y
  //                -3    8.5
  //                 1    4.2
  //                 3    100
  //                -2     13
  //                2.1     0
  //
  // Cond points:    X      Y
  //                -1    2.1
  //                 4    3.2
  //                 3    111
  //                -8     32
  //               0.3      7

  float source[15] = {0.4, 1, -4,  1,   1, 0.2,  98,  12, 1.2, 1.3, 0, 13, 7, -1, 0};
  float dest[10]   = { -3, 1,  3, -2, 2.1, 8.5, 4.2, 100,  13,   0};
  float cond[10]   = { -1, 4,  3, -8, 0.3, 2.1, 3.2, 111,  32,   7};

  float result1[6];
  float result2[2];

  jidt_error_t err;
  err = CMIKraskov_C(N, source, dimx, dest, dimy, cond, dimz, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1, 1);

  EXPECT(err == JIDT_SUCCESS);

  float CMI1 = cpuDigamma(k) + result1[0]/((double) N);

  err = CMIKraskov_C(N, source, dimx, dest, dimy, cond, dimz, k, thelier,
      1, returnLocals, useMaxNorm, isAlgorithm1, result2, 1);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(CMI1 == approx(result2[0]));
  EXPECT(result2[0] != result2[1]);
},

CASE("Basic digamma sum test")
{
  int N = 2;
  int nx[2] = {10, 12};
  int ny[2] = {2, 5};
  float sumDiGammas;

  int err = cudaSumDigammas(&sumDiGammas, nx, ny, N, 1);

  EXPECT(err == 1);

  float cpu_sumDigammas = 0;
  for (int i = 0; i < N; i++) {
    cpu_sumDigammas += cpuDigamma(nx[i] + 1) + cpuDigamma(ny[i] + 1);
  }

  EXPECT(sumDiGammas == approx(cpu_sumDigammas));
  
},

CASE("Digamma sum test with more data (but still one GPU block)")
{
  int N = 500;
  int *nx = (int *) malloc(N * sizeof(int));
  int *ny = (int *) malloc(N * sizeof(int));
  float sumDiGammas;

  for (int i = 0; i < N; i++) {
    nx[i] = rand() % 300;
    ny[i] = rand() % 300;
  }

  int err = cudaSumDigammas(&sumDiGammas, nx, ny, N, 1);

  EXPECT(err == 1);

  float cpu_sumDigammas = 0;
  for (int i = 0; i < N; i++) {
    cpu_sumDigammas += cpuDigamma(nx[i] + 1) + cpuDigamma(ny[i] + 1);
  }

  free(nx); free(ny);

  EXPECT(sumDiGammas == approx(cpu_sumDigammas));
  
},

CASE("Digamma sum test in single chunk with multiple GPU blocks")
{
  int N = 5000;
  int *nx = (int *) malloc(N * sizeof(int));
  int *ny = (int *) malloc(N * sizeof(int));
  float sumDiGammas;

  for (int i = 0; i < N; i++) {
    nx[i] = rand() % 300;
    ny[i] = rand() % 300;
  }

  int err = cudaSumDigammas(&sumDiGammas, nx, ny, N, 1);

  EXPECT(err == 1);

  float cpu_sumDigammas = 0;
  for (int i = 0; i < N; i++) {
    cpu_sumDigammas += cpuDigamma(nx[i] + 1) + cpuDigamma(ny[i] + 1);
  }

  free(nx); free(ny);

  EXPECT(sumDiGammas == approx(cpu_sumDigammas));
  
},

CASE("Test GPU parallel digamma calculation (no sum) in single block")
{
  int N = 500;
  int *nx = (int *) malloc(N * sizeof(int));
  int *ny = (int *) malloc(N * sizeof(int));
  float *gpu_digammas = (float *) malloc(N * sizeof(float));

  for (int i = 0; i < N; i++) {
    nx[i] = rand() % 300;
    ny[i] = rand() % 300;
  }

  int err = parallelDigammas(gpu_digammas, nx, ny, N);

  EXPECT(err == 1);

  for (int i = 0; i < N; i++) {
    double cpu_digamma = cpuDigamma(nx[i] + 1) + cpuDigamma(ny[i] + 1);
    EXPECT(gpu_digammas[i] == approx(cpu_digamma));
  }

  free(nx); free(ny); free(gpu_digammas);
},

CASE("Test GPU parallel digamma calculation (no sum) in multiple blocks")
{
  int N = 5000;
  int *nx = (int *) malloc(N * sizeof(int));
  int *ny = (int *) malloc(N * sizeof(int));
  float *gpu_digammas = (float *) malloc(N * sizeof(float));

  for (int i = 0; i < N; i++) {
    nx[i] = rand() % 300;
    ny[i] = rand() % 300;
  }

  int err = parallelDigammas(gpu_digammas, nx, ny, N);

  EXPECT(err == 1);

  for (int i = 0; i < N; i++) {
    double cpu_digamma = cpuDigamma(nx[i] + 1) + cpuDigamma(ny[i] + 1);
    EXPECT(gpu_digammas[i] == approx(cpu_digamma));
  }

  free(nx); free(ny); free(gpu_digammas);
},

CASE("Digamma sum with multiple chunks and one block")
{
  int nx[6] = {4, 6, 8, 10, 12, 14};
  int ny[6] = {3, 5, 7,  9, 11, 13};
  int trialLength = 3;
  int nchunks = 2;

  float *gpu_sumDiGammas = (float *) malloc(nchunks * sizeof(float));
  int err = cudaSumDigammas(gpu_sumDiGammas, nx, ny, trialLength, nchunks);
  EXPECT(err == 1);

  for (int i = 0; i < nchunks; i++) {
    float cpu_sumDigammas = 0;
    for (int j = 0; j < trialLength; j++) {
      cpu_sumDigammas += cpuDigamma(nx[trialLength*i + j] + 1)
                          + cpuDigamma(ny[trialLength*i + j] + 1);
    }
    EXPECT(gpu_sumDiGammas[i] == approx(cpu_sumDigammas));
  }

  free(gpu_sumDiGammas);
},

CASE("Digamma sum with multiple chunks and multiple blocks per chunk")
{
  int trialLength = 2000;
  int nchunks = 3;
  int signalLength = trialLength * nchunks;
  int *nx = (int *) malloc(signalLength * sizeof(int));
  int *ny = (int *) malloc(signalLength * sizeof(int));

  for (int i = 0; i < signalLength; i++) {
    nx[i] = rand() % 300;
    ny[i] = rand() % 300;
  }

  float *gpu_sumDiGammas = (float *) malloc(nchunks * sizeof(float));
  int err = cudaSumDigammas(gpu_sumDiGammas, nx, ny, trialLength, nchunks);
  EXPECT(err == 1);

  for (int i = 0; i < nchunks; i++) {
    float cpu_sumDigammas = 0;
    for (int j = 0; j < trialLength; j++) {
      cpu_sumDigammas += cpuDigamma(nx[trialLength*i + j] + 1)
                         + cpuDigamma(ny[trialLength*i + j] + 1);
    }
    EXPECT(gpu_sumDiGammas[i] == approx(cpu_sumDigammas));
  }
  free(nx); free(ny); free(gpu_sumDiGammas);
},

};

int main(int argc, char *argv[]){
  if (int failures = lest::run(specification, argc, argv))
    return failures;
  return std::cout << "ALL PASSED, YEA BOI\n", EXIT_SUCCESS;
}

