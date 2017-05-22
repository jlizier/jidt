#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "lest.hpp"

#include "gpuKnnLibrary.h"
#include "gpuMILibrary.h"
#include "digamma.h"

using lest::approx;

const lest::test specification[] = {

CASE("Basic kNN test")
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

CASE("Test kNN with multiple trials in the same call")
{
  int thelier = 0;
  int useMaxNorm = 1;
  int N = 8;
  int dims = 2;
  int k = 1;
  int nchunks = 2;
  int indexes[8];
  float distances[8];
  float pointset[16] = {-1, 0.5, 1.1, 2, -1.1, -2, -1.6, -1,
                        -1, 0.5, 1.1, 2,   -1, -2, -1.6, -1};

  int err = cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
              nchunks, dims, N, useMaxNorm);

  EXPECT(err == 1);

  EXPECT(indexes[0] == 1);
  EXPECT(indexes[1] == 2);
  EXPECT(indexes[2] == 1);
  EXPECT(indexes[3] == 2);
  EXPECT(indexes[4] == 3);
  EXPECT(indexes[5] == 2);
  EXPECT(indexes[6] == 1);
  EXPECT(indexes[7] == 0);

  EXPECT(distances[0] == approx(1.5));
  EXPECT(distances[1] == approx(0.6));
  EXPECT(distances[2] == approx(0.6));
  EXPECT(distances[3] == approx(0.9));
  EXPECT(distances[4] == approx(0.1));
  EXPECT(distances[5] == approx(0.4));
  EXPECT(distances[6] == approx(0.4));
  EXPECT(distances[7] == approx(0.1));

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

CASE("Test correct pointset arrangement")
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

  printf("===============================\n");
  err = MIKraskov_C(N, source, dimx, dest, dimy, k, thelier,
      0, returnLocals, useMaxNorm, isAlgorithm1, result1);

  printf("===============================\n");
  err = MIKraskovByPointsetChunks(N*2, source, dimx, dest, dimy, k, thelier,
      2, returnLocals, useMaxNorm, isAlgorithm1, result2, double_pointset);

  float MI1 = cpuDigamma(k) + cpuDigamma(N) - result1[0]/((double) N);

  EXPECT(err == JIDT_SUCCESS);
  EXPECT(result2[0] == approx(MI1));
  EXPECT(result2[0] == approx(result2[1]));

},

CASE("Test that two samples with same joints in two chunks give same result")
{
  EXPECT(1 == 1);

},

};

int main(int argc, char *argv[]){
  if (int failures = lest::run(specification, argc, argv))
    return failures;
  return std::cout << "ALL PASSED, YEA BOI\n", EXIT_SUCCESS;
}

