#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "gpuMILibrary.h"
#include "gpuKnnLibrary.h"
#include "digamma.h"
#include "ctimer.h"


jidt_error_t MIKraskov_C(int N, float *source, int dimx, float *dest, int dimy,
    int k, int thelier, int nb_surrogates, int returnLocals, int useMaxNorm,
    int isAlgorithm1, float *result) {
  return MIKraskovWithReorderings(N, source, dimx, dest, dimy, k, thelier, nb_surrogates,
      returnLocals, useMaxNorm, isAlgorithm1, result, 0, NULL);
}

/**
 * Calculate Mutual Information using the KSG algorithm.
 */
jidt_error_t MIKraskovWithReorderings(int N, float *source, int dimx, float *dest, int dimy,
    int k, int thelier, int nb_surrogates, int returnLocals, int useMaxNorm,
    int isAlgorithm1, float *result, int reorderingsGiven, int **reorderings) {

  CPerfTimer pt = startTimer("Rearranging pointset");

  // Allocate more space if surrogates are requested
  int nchunks = nb_surrogates + 1;
  int dims = dimx + dimy;
  float *pointset = (float *) malloc(N * dims * nchunks * sizeof(float));

  if (nb_surrogates == 0) {
    memcpy(                 pointset, source, N*dimx*sizeof(float));
    memcpy(pointset + nchunks*N*dimx,   dest, N*dimy*sizeof(float));
  }

  if (nb_surrogates > 0) {

     for (int i = 0; i < N; i++) {
       for (int j = 0; j < dimx; j++) {
         pointset[j*N*nchunks+i] = source[N*j+i];
       }

       for (int j = 0; j < dimy; j++) {
         pointset[nchunks*N*dimx + j*N*nchunks + i] = dest[N*j+i];
       }
     }

    // If surrogates requested, copy permutations as well
    int *order;
    int perm[N];
    if (!reorderingsGiven) {
      for (int i = 0; i < N; i++) {
        perm[i] = i;
      }
    }

    for (int s = 0; s < nb_surrogates; s++) {
      if (reorderingsGiven) {
        order = reorderings[s];
      } else {
        randperm(perm, N);
        order = perm;
      }

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < dimx; j++) {
          pointset[(s+1)*N + N*j*nchunks + i] = source[N*j + order[i]];
        }

        for (int j = 0; j < dimy; j++) {
          pointset[nchunks*N*dimx + (s+1)*N + N*j*nchunks + i] = dest[N*j + i];
        }
      }
    }
  }
  stopTimer(pt);

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

  float *d_source, *d_dest, *d_distances, *d_radii, *d_digammas;
  int *d_nx, *d_ny, *d_indexes;

  {
  CPerfTimer pt = startTimer("GPU_warmup");
  gpuWarmUp();
  stopTimer(pt);
  }

  // 1. Allocate space in GPU and transfer memory
  // ======================
  allocateDeviceMemory(signalLength, k, dimx, dimy, &d_source, &d_dest, &d_distances,
      &d_indexes, &d_radii, &d_nx, &d_ny, &d_digammas, pointset);

  // 2. Find nearest neighbours
  // ======================
  {
  CPerfTimer pt = startTimer("kNN_full");
  d_cudaFindKnn(d_indexes, d_distances, d_source, d_source, k,
      thelier, nchunks, dims, signalLength, useMaxNorm);
  stopTimer(pt);
  }

  // 4. Count points strictly within R in the X-space
  // ======================
  {
  CPerfTimer pt = startTimer("RS_full");
  d_cudaFindRSAll(d_nx, d_source, d_source, d_radii, thelier, nchunks, dimx, signalLength, useMaxNorm);
  d_cudaFindRSAll(d_ny, d_dest, d_dest, d_radii, thelier, nchunks, dimy, signalLength, useMaxNorm);
  stopTimer(pt);
  }

  // 6. Set locals, surrogates or digammas for return
  // ======================
  {
  CPerfTimer pt = startTimer("Digammas_full");
  if (returnLocals) {
    float digammaK = cpuDigamma(k);
    float digammaN = cpuDigamma(trialLength);
    float digammas[trialLength];
    d_parallelDigammas(digammas, d_digammas, d_nx, d_ny, signalLength);
    for (int i = 0; i < trialLength; i++) {
      result[i] = digammaK + digammaN - digammas[i];
    }

  } else {

    float digammaK = cpuDigamma(k);
    float digammaN = cpuDigamma(trialLength);
    float sumDigammas[nchunks];
    d_cudaSumDigammas(sumDigammas, d_nx, d_ny, d_digammas, trialLength, nchunks);

    if (nchunks > 1) {
      for (int ii = 0; ii < nchunks; ii++) {
        result[ii] = digammaK + digammaN - sumDigammas[ii]/((float) trialLength);
      }
    } else {
      result[0] = sumDigammas[0];
      result[1] = -1;
      result[2] = -1;
    }

  }
  stopTimer(pt);
  }

  err = JIDT_SUCCESS;

  freeDeviceMemory(d_source);

  return err;
}

