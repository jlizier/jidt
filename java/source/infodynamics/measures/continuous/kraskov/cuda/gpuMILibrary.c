#include <stdlib.h>
#include <stdio.h>

#include "gpuMILibrary.h"
#include "gpuKnnLibrary.h"
#include "digamma.h"


/**
 * Make random permutation of perm[].
 *
 * @param perm preallocated and prefilled integer array to be shuffled
 * @param n number of elements in perm
 */
void randperm(int perm[], int n) {
  // Random permutation the order
  for (int i = 0; i < n; i++) {
   int j, t;
   j = rand() % (n-i) + i;
   t = perm[j]; perm[j] = perm[i]; perm[i] = t; // Swap i and j
  }
}


/**
 * Calculate Mutual Information using the KSG algorithm.
 */
jidt_error_t MIKraskov_C(int N, float *source, int dimx, float *dest, int dimy,
    int k, int thelier, int nb_surrogates, int returnLocals, int useMaxNorm,
    int isAlgorithm1, float *result) {

  // Allocate more space if surrogates are requested
  int nchunks = nb_surrogates + 1;
  int dims = dimx + dimy;
  float *pointset = (float *) malloc(N * dims * nchunks * sizeof(float));

  if (nb_surrogates == 0) {
    // If no surrogates requested, copy source and dest as usual
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < dimx; j++) {
        register int idx = N*j + i;
        pointset[idx] = source[idx];
      }

      for (int j = 0; j < dimy; j++) {
        register int idx = N*j + i;
        pointset[N*dimx + idx] = dest[idx];
      }
    }

  } else {
    // If surrogates requested, copy including random permutations
    int perm[N];
    for (int i = 0; i < N; i++) {
      perm[i] = i;
    }

    for (int chunk = 0; chunk < nchunks; chunk++) {
      randperm(perm, N);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < dimx; j++) {
          register int idx = N*j + i;
          pointset[N*dims + idx] = source[idx];
        }

        for (int j = 0; j < dimy; j++) {
          register int idx = N*j + i;
          pointset[N*dims + N*dimx + idx] = dest[idx];
        }
      }
    }
  }

  jidt_error_t err =  MIKraskovByPointsetChunks(N*nchunks, source, dimx,
                                                dest, dimy, k, thelier,
                                                nchunks, returnLocals, useMaxNorm,
                                                isAlgorithm1, result, pointset);
  FREE(pointset);

  return err;
}


jidt_error_t MIKraskovByPointsetChunks(int signalLength, float *source, int dimx,
    float *dest, int dimy, int k, int thelier, int nchunks,
    int returnLocals, int useMaxNorm, int isAlgorithm1, float *result,
    float *pointset) {

  int dims = dimx + dimy;
  int err;
  int trialLength = signalLength/((float) nchunks);
  float *distances = NULL, *radii = NULL;
  int *indexes = NULL, *nx = NULL, *ny = NULL;

  // 1. Re-arrange data and generate shufflings
  // ======================


  // 2. Find nearest neighbours
  // ======================
  indexes = (int *) malloc(signalLength * k * sizeof(int));
  distances = (float *) malloc(signalLength * k * sizeof(float));
  if (!cudaFindKnn(indexes, distances, pointset, pointset, k,
          thelier, nchunks, dims, signalLength, useMaxNorm)) {
    device_reset();
    return JIDT_ERROR;
  }


  // 3. Find distance R from each point to its k-th neighbour in joint space
  // ======================
  
  radii = (float *) malloc(signalLength * sizeof(float));
  for (int i = 0; i < signalLength; i++) {
    radii[i] = distances[signalLength*(k-1) + i];
  }


  // 4. Count points strictly within R in the X-space
  // ======================

  // If we're using algorithm 2 then we need the radius in the marginal space,
  // not the joint (as calculated above)
  if (!isAlgorithm1) {
    findRadiiAlgorithm2(radii, source, indexes, k, dimx, signalLength);
  }

  nx = (int *) malloc(signalLength * sizeof(int));
  if (!cudaFindRSAll(nx, source, source, radii, thelier, nchunks, dimx, signalLength, useMaxNorm)) {
    device_reset();
    return JIDT_ERROR;
  }

  if (!isAlgorithm1) {
    findRadiiAlgorithm2(radii, dest, indexes, k, dimy, signalLength);
  }

  ny = (int *) malloc(signalLength * sizeof(int));
  if (!cudaFindRSAll(ny, dest, dest, radii, thelier, nchunks, dimy, signalLength, useMaxNorm)) {
    device_reset();
    return JIDT_ERROR;
  }


  // 6. Set locals, surrogates or digammas for return
  // ======================

  // Check this for digamma parallelisation
  // https://devtalk.nvidia.com/default/topic/516516/kernel-launch-failure-in-matlab/?offset=2

  if (nchunks > 1) {
    float digammaK = cpuDigamma(k);
    float digammaN = cpuDigamma(trialLength);
    for (int ii = 0; ii < nchunks; ii++) {
      float sumDiGammas = 0;
      for (int i = 0; i < trialLength; i++) {
        sumDiGammas += cpuDigamma(nx[i] + 1) + cpuDigamma(ny[i] + 1);
      }
      result[ii] = digammaK + digammaN - sumDiGammas/trialLength;
    }

  } else if (returnLocals) {
    float digammaK = cpuDigamma(k);
    float digammaN = cpuDigamma(trialLength);
    for (int i = 0; i < trialLength; i++) {
      result[i] = digammaK + digammaN - cpuDigamma(nx[i] + 1) - cpuDigamma(ny[i] + 1);
    }

  } else {
    float sumNx = 0;
    float sumNy = 0;
    float sumDiGammas = 0;
    for (int i = 0; i < trialLength; i++) {
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

  return err;
}

