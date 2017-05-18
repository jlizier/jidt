#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "lest.hpp"

using lest::approx;

extern "C" int cudaFindKnn(int* h_bf_indexes, float* h_bf_distances, float* h_pointset,
    float* h_query, int kth, int thelier, int nchunks, int pointdim,
    int signallength, unsigned int useMaxNorm);

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
  float pointset[8] = {-1, 0.5, 1.1, 2, -1, 0.5, 1.1, 2};

  cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
      nchunks, dims, N, useMaxNorm);

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
  float pointset[8] = {-1, 0.5, 1.1, 2, -1, 0.5, 1.1, 2};

  cudaFindKnn(indexes, distances, pointset, pointset, k, thelier,
      nchunks, dims, N, useMaxNorm);

  EXPECT(indexes[0] == 1);
  EXPECT(indexes[1] == 2);
  EXPECT(indexes[2] == 1);
  EXPECT(indexes[3] == 2);

  EXPECT(distances[0] == approx(4.5));
  EXPECT(distances[1] == approx(0.72));
  EXPECT(distances[2] == approx(0.72));
  EXPECT(distances[3] == approx(1.62));

},

  //---------------------------------------------
  // Test kNN with more points far away
  

  //---------------------------------------------
  // Test kNN 


  //---------------------------------------------
  // Test kNN in higher dimensions, padding with zeros


  //---------------------------------------------
  // Test kNN in higher dimensions, padding with zeros

};

int main(int argc, char *argv[]){
  return lest::run( specification, argc, argv );
}

