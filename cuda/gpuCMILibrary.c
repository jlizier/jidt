#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "gpuCMILibrary.h"
#include "gpuKnnLibrary.h"
#include "digamma.h"
#include "ctimer.h"


jidt_error_t CMIKraskov_C(int N, float *source, int dimx, float *dest, int dimy,
    float *cond, int dimz, int k, int thelier, int nb_surrogates,
    int returnLocals, int useMaxNorm, int isAlgorithm1, float *result,
    int variableToReorder) {
  return CMIKraskovWithReorderings(N, source, dimx, dest, dimy, cond, dimz,
      k, thelier, nb_surrogates, returnLocals, useMaxNorm, isAlgorithm1, result,
      0, NULL, variableToReorder);
}

/**
 * Calculate Mutual Information using the KSG algorithm.
 */
jidt_error_t CMIKraskovWithReorderings(int N, float *source, int dimx,
    float *dest, int dimy, float *cond, int dimz, int k, int thelier,
    int nb_surrogates, int returnLocals, int useMaxNorm,
    int isAlgorithm1, float *result, int reorderingsGiven, int **reorderings,
    int variableToReorder) {

  CPerfTimer pt = startTimer("Rearranging pointset");

  // Allocate more space if surrogates are requested
  int nchunks = nb_surrogates + 1;
  int dims = dimx + dimy + dimz;
  float *pointset = (float *) malloc(N * dims * nchunks * sizeof(float));

  if (nb_surrogates == 0) {
    memcpy(                pointset, source, N*dimx*sizeof(float));
    memcpy(       pointset + N*dimx,   cond, N*dimz*sizeof(float));
    memcpy(pointset + N*(dimx+dimz),   dest, N*dimy*sizeof(float));
  }

  if (nb_surrogates > 0) {

     for (int i = 0; i < N; i++) {
       for (int j = 0; j < dimx; j++) {
         pointset[j*N*nchunks+i] = source[N*j+i];
       }

       for (int j = 0; j < dimz; j++) {
         pointset[nchunks*N*dimx + j*N*nchunks + i] = cond[N*j+i];
       }

       for (int j = 0; j < dimy; j++) {
         pointset[nchunks*N*(dimx+dimz) + j*N*nchunks + i] = dest[N*j+i];
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
        if (variableToReorder == 1) {
          for (int j = 0; j < dimx; j++) {
            pointset[(s+1)*N + N*j*nchunks + i] = source[N*j + order[i]];
          }

          for (int j = 0; j < dimz; j++) {
            pointset[nchunks*N*dimx + (s+1)*N + N*j*nchunks + i] = cond[N*j + i];
          }

        } else {
          for (int j = 0; j < dimx; j++) {
            pointset[(s+1)*N + N*j*nchunks + i] = source[N*j + i];
          }

          for (int j = 0; j < dimz; j++) {
            pointset[nchunks*N*dimx + (s+1)*N + N*j*nchunks + i] = cond[N*j + order[i]];
          }

        }

        for (int j = 0; j < dimy; j++) {
          pointset[nchunks*N*(dimx+dimz) + (s+1)*N + N*j*nchunks + i] = dest[N*j + i];
        }
      }
    }
  }
  stopTimer(pt);

  jidt_error_t err =  CMIKraskovByPointsetChunks(N*nchunks, source, dimx,
                                                dest, dimy, cond, dimz, k, thelier,
                                                nchunks, returnLocals, useMaxNorm,
                                                isAlgorithm1, result, pointset);
  FREE(pointset);

  return err;
}


jidt_error_t CMIKraskovByPointsetChunks(int signalLength, float *source, int dimx,
    float *dest, int dimy, float *cond, int dimz, int k, int thelier, int nchunks,
    int returnLocals, int useMaxNorm, int isAlgorithm1, float *result,
    float *pointset) {

  int dims = dimx + dimy + dimz;
  int err;
  int trialLength = signalLength/((float) nchunks);

  float *d_source, *d_dest, *d_cond, *d_distances, *d_radii, *d_digammas;
  int *d_nx, *d_ny, *d_nz, *d_indexes;

  {
  CPerfTimer pt = startTimer("GPU_warmup");
  gpuWarmUp();
  stopTimer(pt);
  }

  // 1. Allocate space in GPU and transfer memory
  // ======================
  allocateDeviceMemoryCMI(signalLength, k, dimx, dimy, dimz, &d_source, &d_dest, &d_cond,
      &d_distances, &d_indexes, &d_radii, &d_nx, &d_ny, &d_nz, &d_digammas, pointset);

  // 2. Find nearest neighbours in joint space
  // ======================
  {
  CPerfTimer pt = startTimer("kNN_full");
  d_cudaFindKnn(d_indexes, d_distances, d_source, d_source, k,
      thelier, nchunks, dims, signalLength, useMaxNorm);
  stopTimer(pt);
  }

  // 4. Count points strictly within R in the XZ-, YZ- and Z-spaces
  // ======================
  {
  CPerfTimer pt = startTimer("RS_full");
  d_cudaFindRSAll(d_nx, d_source, d_source, d_radii, thelier, nchunks, dimx + dimz, signalLength, useMaxNorm);
  d_cudaFindRSAll(d_ny, d_cond, d_cond, d_radii, thelier, nchunks, dimy + dimz, signalLength, useMaxNorm);
  d_cudaFindRSAll(d_nz, d_cond, d_cond, d_radii, thelier, nchunks, dimz, signalLength, useMaxNorm);
  stopTimer(pt);
  }

  // 6. Set locals, surrogates or digammas for return
  // ======================
  {
  CPerfTimer pt = startTimer("Digammas_full");
  if (returnLocals) {
    float digammaK = cpuDigamma(k);
    float digammas[trialLength];
    d_parallelDigammasCMI(digammas, d_digammas, d_nx, d_ny, d_nz, signalLength);
    for (int i = 0; i < trialLength; i++) {
      result[i] = digammaK - digammas[i];
    }

  } else {

    float digammaK = cpuDigamma(k);
    float sumDigammas[nchunks];
    d_cudaSumDigammasCMI(sumDigammas, d_nx, d_ny, d_nz, d_digammas, trialLength, nchunks);

    if (nchunks > 1) {
      for (int ii = 0; ii < nchunks; ii++) {
        result[ii] = digammaK - sumDigammas[ii]/((float) trialLength);
      }
    } else {
      // Sign changed to comply with the returnValues processing in the Java
      // KSG CMI calc, which is different from the one in the MI calc.
      result[0] = -1 * sumDigammas[0];
      result[1] = -1;
      result[2] = -1;
      result[3] = -1;
      result[4] = -1;
      result[5] = -1;
    }

  }
  stopTimer(pt);
  }

  err = JIDT_SUCCESS;

  freeDeviceMemory(d_source);

  return err;
}

