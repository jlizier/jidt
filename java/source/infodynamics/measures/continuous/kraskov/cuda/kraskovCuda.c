#include <jni.h>
#include <stdlib.h>
#include <stdio.h>

#include "digamma.c"
#include "gpuKnnLibrary.h"

#define check(ans) { _check((ans), __FILE__, __LINE__); }

#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     infodynamics_measures_continuous_kraskov_MutualInfoCalculatorMultiVariateKraskov
 * Method:    MIKraskov
 * Signature: (I[DI[DIIZZZ)[D
 */
JNIEXPORT jdoubleArray JNICALL
  Java_infodynamics_measures_continuous_kraskov_MutualInfoCalculatorMultiVariateKraskov_MIKraskov(
   JNIEnv *env, jobject thisObj, jint j_N,
   jobjectArray j_sourceArray, jint j_dimx,
   jobjectArray j_destArray,   jint j_dimy,
   jint j_k, jboolean j_returnLocals, jboolean j_useMaxNorm, jboolean j_isAlgorithm1) {

  // Declare variables
  // Variables starting with j_ are Java types
  float *source, *dest, *pointset, *radii, *distances;
  int *nx, *ny, *indexes;
  jdouble *jd_source, *jd_dest, *result;
  unsigned int N, dimx, dimy, dims, k, returnLocals, useMaxNorm, isAlgorithm1, result_size;
  float sumNx, sumNy, sumDiGammas;
  unsigned int i, j;
  long double digammaK, digammaN;
  int thelier = 0, nchunks = 1, err;
  const int verbose = 0;

  if (verbose) printf("Verbose is on.\n");

  // Check that incoming data has correct size
  jsize sourceLength = (*env)->GetArrayLength(env, j_sourceArray);
  jsize destLength = (*env)->GetArrayLength(env, j_destArray);

  if (sourceLength != j_N || destLength != j_N) {
    device_reset();
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Data has wrong length.");
  }

  // 1. Load data into GPU
  // =====================
  //   1a. Copy atomic variables from Java
  N = j_N;
  k = j_k;
  dimx = j_dimx;
  dimy = j_dimy;
  dims = dimx + dimy;
  returnLocals = j_returnLocals ? 1 : 0;
  useMaxNorm   = j_useMaxNorm   ? 1 : 0;
  isAlgorithm1 = j_isAlgorithm1 ? 1 : 0;
  
  if (verbose > 0) {
    printf("Loading data from Java to C.\n");
    printf("N = %i, dimx = %i, dimy = %i\n", N, dimx, dimy);
  }

  //  1b. Get arrays from Java
  source   = (float *) malloc(N * dimx * sizeof(float));
  dest     = (float *) malloc(N * dimy * sizeof(float));
  pointset = (float *) malloc(N * dims * sizeof(float));

  if (NULL == source || NULL == dest || NULL == pointset) {
    printf("Error allocating data.\n");
    device_reset();
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Error allocating data.");
  }

  for (i = 0; i < N; i++) {
    jdoubleArray j_sourceRow = (jdoubleArray) (*env)->GetObjectArrayElement(env, j_sourceArray, i);
    jdoubleArray j_destRow   = (jdoubleArray) (*env)->GetObjectArrayElement(env, j_destArray, i);

    jdouble *sourceRow = (*env)->GetDoubleArrayElements(env, j_sourceRow, NULL);
    jdouble *destRow = (*env)->GetDoubleArrayElements(env, j_destRow, NULL);

    // Data in java are doubles, but GPUs need floats.
    // We have to cast them manually, so we can't memcopy

    // The following two for-loops get two matrices in T-by-D indexing (i.e.
    // first dimension is time, second is variable) and return the data in
    // column-major form
    for (j = 0; j < dimx; j++) {
      register float x = (float) sourceRow[j];
      source[N*j + i]   = x;
      pointset[N*j + i] = x;
    }

    for (j = 0; j < dimy; j++) {
      register float x = (float) destRow[j];
      dest[N*j + i] = x;
      pointset[N*dimx + N*j + i] = x;
    }

    // // The following two for-loops get the same matrices but return arrays
    // // in row-major form. Code is left here for future development purposes.
    // for (j = 0; j < dimx; j++) {
    //   register float x = (float) sourceRow[j];
    //   source[dimx*i + j]   = x;
    //   pointset[dims*i + j] = x;
    // }

    // for (j = 0; j < dimy; j++) {
    //   register float x = (float) destRow[j];
    //   dest[dimy*i + j] = x
    //   pointset[i*dims + dimx + j] = x;
    // }


    (*env)->ReleaseDoubleArrayElements(env, j_sourceRow, sourceRow, 0);
    (*env)->ReleaseDoubleArrayElements(env, j_destRow, destRow, 0);

    (*env)->DeleteLocalRef(env, j_sourceRow);
    (*env)->DeleteLocalRef(env, j_destRow);

    // FIXME: I'm not entirely sure I'm freeing all the memory here. I should
    // check for memory leaks more carefully.

  }

  if (verbose > 1) {
    printf("Flat source array: ");
    for (i = 0; i < N*dimx; i++) printf("%f, ", source[i]);
    printf("\n");
    printf("Flat dest array: ");
    for (i = 0; i < N*dimy; i++) printf("%f, ", dest[i]);
    printf("\n");
    printf("Flat pointset: ");
    for (i = 0; i < N*dims; i++) printf("%f, ", pointset[i]);
    printf("\n");
  }


  // 2. Get indices and distances to the k NN of each point
  // ======================

  indexes = (int *) malloc(N * k * sizeof(int));
  distances = (float *) malloc(N * k * sizeof(float));

  err = cudaFindKnn(indexes, distances, pointset, pointset, k,
      thelier, nchunks, dims, N, useMaxNorm);
  if (err != 1) {
    device_reset();
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Cuda error finding nearest neighbours.");
  }

  if (verbose > 1) {
    printf("Distances: \n");
    for (j = 0; j < N; j++) {
      for (i = 0; i < k; i++) {
        printf("%f\t", distances[j + i*N]);
      }
      printf("\n");
    }
    printf("Indexes: \n");
    for (j = 0; j < N; j++) {
      for (i = 0; i < k; i++) {
        printf("%i\t", indexes[j + i*N]);
      }
      printf("\n");
    }
  }


  // 3. Find distance R from each point to its k-th neighbour in joint space
  // ======================
  
  if (verbose > 0) {
    printf("Getting radii from calculated dists.\n");
  }

  radii = (float *) malloc(N * sizeof(float));
  for (i = 0; i < N; i++) {
    radii[i] = distances[N*(k-1) + i];
  }

  if (verbose > 1) {
    for (i = 0; i < N; i++) {
      printf("%f\t", radii[i]);
    }
    printf("\n");
  }


  // 4. Count points strictly within R in the X-space
  // ======================

  if (verbose > 0) {
    printf("Counting points within radius.\n");
  }

  // If we're using algorithm 2 then we need the radius in the marginal space,
  // not the joint (as calculated above)
  if (!isAlgorithm1) {
    findRadiiAlgorithm2(radii, source, indexes, k, dimx, N);
  }

  nx = (int *) malloc(N * sizeof(int));
  err = cudaFindRSAll(nx, source, source, radii, thelier, nchunks, dimx, N, useMaxNorm);
  if (err != 1) {
    device_reset();
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Cuda error during source range search.");
  }

  if (verbose > 0) {
    printf("X count: ");
    for (i = 0; i < N; i++) printf("%i, ", nx[i]);
    printf("\n\n");
  }

  // 4. Count points strictly within R in the X-space
  // ======================

  if (!isAlgorithm1) {
    findRadiiAlgorithm2(radii, dest, indexes, k, dimy, N);
  }

  ny = (int *) malloc(N * sizeof(int));
  err = cudaFindRSAll(ny, dest, dest, radii, thelier, nchunks, dimy, N, useMaxNorm);
  if (err != 1) {
    device_reset();
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Cuda error during dest range search.");
  }

  if (verbose > 0) {
    printf("Y count: ");
    for (i = 0; i < N; i++) printf("%i, ", ny[i]);
    printf("\n\n");
  }


  // 6. Set locals or digammas for return
  // ======================
  if (verbose > 0) {
    printf("Start calculating digammas\n");
  }

  if (returnLocals) {
    result = (jdouble *) malloc(N * sizeof(jdouble));
    result_size = N;
    digammaK = cpuDigamma(k);
    digammaN = cpuDigamma(N);
    for (i = 0; i < N; i++) {
      result[i] = digammaK + digammaN - cpuDigamma(nx[i] + 1) - cpuDigamma(ny[i] + 1);
    }

  } else {

    result = (jdouble *) malloc(3 * sizeof(jdouble));
    result_size = 3;

    sumNx = 0;
    sumNy = 0;
    sumDiGammas = 0;
    for (i = 0; i < N; i++) {
      sumNx       += nx[i];
      sumNy       += ny[i];
      sumDiGammas += cpuDigamma(nx[i] + 1) + cpuDigamma(ny[i] + 1);

    }

    result[0] = (jdouble) sumDiGammas;
    result[1] = (jdouble) sumNx;
    result[2] = (jdouble) sumNy;

  }


  // 7. Free memory and return
  // =========================
  
  if (verbose > 0) {
    printf("Free memory\n");
  }

  if (source)    free(source);
  if (dest)      free(dest);
  if (pointset)  free(pointset);
  if (radii)     free(radii);
  if (distances) free(distances);
  if (indexes) free(indexes);
  if (nx)   free(nx);
  if (ny)   free(ny);

  // Set Java array for return
  jdoubleArray outJNIArray = (*env)->NewDoubleArray(env, result_size);  // allocate
  if (NULL == outJNIArray) {
    device_reset();
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Error creating return array.");
  }
  (*env)->SetDoubleArrayRegion(env, outJNIArray, 0 , result_size, result);  // copy

  if (result) free(result);
  
  return outJNIArray;


} // End of function MIKraskov1


#ifdef __cplusplus
}
#endif

