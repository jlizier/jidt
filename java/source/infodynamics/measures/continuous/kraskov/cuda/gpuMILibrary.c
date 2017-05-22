#include <stdlib.h>
#include <stdio.h>

#include "gpuMILibrary.h"
#include "gpuKnnLibrary.h"
#include "digamma.h"

jidt_error_t MIKraskov_C(int N, float *source, int dimx, float *dest, int dimy,
    int k, int thelier, int nchunks, int returnLocals, int useMaxNorm,
    int isAlgorithm1, float *result) {

  int dims = dimx + dimy;
  int err;
  int i, j;
  const int verbose = 0;
  float *pointset = NULL, *distances = NULL, *radii = NULL;
  int *indexes = NULL, *nx = NULL, *ny = NULL;

  pointset = (float *) malloc(N * dims * sizeof(float));
  for (i = 0; i < N; i++) {
    for (j = 0; j < dimx; j++) {
      register int idx = N*j + i;
      pointset[idx] = source[idx];
    }

    for (j = 0; j < dimy; j++) {
      register int idx = N*j + i;
      pointset[N*dimx + idx] = dest[idx];
    }
  }

  indexes = (int *) malloc(N * k * sizeof(int));
  distances = (float *) malloc(N * k * sizeof(float));
  if (!cudaFindKnn(indexes, distances, pointset, pointset, k,
          thelier, nchunks, dims, N, useMaxNorm)) {
    device_reset();
    return JIDT_ERROR;
  }


  // 3. Find distance R from each point to its k-th neighbour in joint space
  // ======================
  
  radii = (float *) malloc(N * sizeof(float));
  for (i = 0; i < N; i++) {
    radii[i] = distances[N*(k-1) + i];
  }


  // 4. Count points strictly within R in the X-space
  // ======================

  // If we're using algorithm 2 then we need the radius in the marginal space,
  // not the joint (as calculated above)
  if (!isAlgorithm1) {
    findRadiiAlgorithm2(radii, source, indexes, k, dimx, N);
  }

  nx = (int *) malloc(N * sizeof(int));
  if (!cudaFindRSAll(nx, source, source, radii, thelier, nchunks, dimx, N, useMaxNorm)) {
    device_reset();
    return JIDT_ERROR;
  }

  if (!isAlgorithm1) {
    findRadiiAlgorithm2(radii, dest, indexes, k, dimy, N);
  }

  ny = (int *) malloc(N * sizeof(int));
  if (!cudaFindRSAll(ny, dest, dest, radii, thelier, nchunks, dimy, N, useMaxNorm)) {
    device_reset();
    return JIDT_ERROR;
  }


  // 6. Set locals, surrogates or digammas for return
  // ======================

  float sumNx, sumNy, sumDiGammas, digammaK, digammaN;
  int result_size;

  if (nchunks > 1) {
    float digammaK = cpuDigamma(k);
    float digammaN = cpuDigamma(N);
    float trialLength = N/((float) nchunks);
    for (int ii = 0; ii < nchunks; ii++) {
      float sumDiGammas = 0;
      for (i = 0; i < trialLength; i++) {
        sumDiGammas += cpuDigamma(nx[i] + 1) + cpuDigamma(ny[i] + 1);
      }
      result[ii] = digammaK - sumDiGammas/trialLength + digammaN;
    }

  } else if (returnLocals) {
    float digammaK = cpuDigamma(k);
    float digammaN = cpuDigamma(N);
    for (i = 0; i < N; i++) {
      result[i] = digammaK + digammaN - cpuDigamma(nx[i] + 1) - cpuDigamma(ny[i] + 1);
    }

  } else {
    float sumNx = 0;
    float sumNy = 0;
    float sumDiGammas = 0;
    for (i = 0; i < N; i++) {
      sumNx       += nx[i];
      sumNy       += ny[i];
      sumDiGammas += cpuDigamma(nx[i] + 1) + cpuDigamma(ny[i] + 1);
    }

    result[0] = sumDiGammas;
    result[1] = sumNx;
    result[2] = sumNy;

  }

  err = JIDT_SUCCESS;

  FREE(indexes);
  FREE(distances);
  FREE(nx);
  FREE(ny);
  FREE(distances);
  FREE(radii);
  FREE(pointset);

  return err;
}

