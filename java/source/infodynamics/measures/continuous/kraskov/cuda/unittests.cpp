#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "lest.hpp"

#include "gpuKnnLibrary.h"
#include "gpuMILibrary.h"

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

};

int main(int argc, char *argv[]){
  if (int failures = lest::run(specification, argc, argv))
    return failures;
  return std::cout << "ALL PASSED, YEA BOI\n", EXIT_SUCCESS;
}

