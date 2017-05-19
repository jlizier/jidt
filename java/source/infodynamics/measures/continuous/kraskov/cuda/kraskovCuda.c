#include <jni.h>
#include <stdlib.h>
#include <stdio.h>

#include "digamma.c"
#include "gpuKnnLibrary.h"

#define check(ans) { _check((ans), __FILE__, __LINE__); }

typedef enum { JIDT_SUCCESS, JIDT_ERROR } jidt_error_t;

jidt_error_t MIKraskov_C(int N, float *source, int dimx, float *dest, int dimy,
    int k, int thelier, int nchunks, int returnLocals, int useMaxNorm,
    int isAlgorithm1, int nb_surrogates, float *result) {

  int dims = dimx + dimy;
  int err;
  int i, j;
  const int verbose = 0;
  nb_surrogates = 0;

  float *pointset = (float *) malloc(N * dims * sizeof(float));
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

  int *indexes = (int *) malloc(N * k * sizeof(int));
  float *distances = (float *) malloc(N * k * sizeof(float));
  if (!cudaFindKnn(indexes, distances, pointset, pointset, k,
          thelier, nchunks, dims, N, useMaxNorm)) {
    device_reset();
    return JIDT_ERROR;
  }


  // 3. Find distance R from each point to its k-th neighbour in joint space
  // ======================
  
  float *radii = (float *) malloc(N * sizeof(float));
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

  int *nx = (int *) malloc(N * sizeof(int));
  if (!cudaFindRSAll(nx, source, source, radii, thelier, nchunks, dimx, N, useMaxNorm)) {
    device_reset();
    return JIDT_ERROR;
  }

  if (!isAlgorithm1) {
    findRadiiAlgorithm2(radii, dest, indexes, k, dimy, N);
  }

  int *ny = (int *) malloc(N * sizeof(int));
  if (!cudaFindRSAll(ny, dest, dest, radii, thelier, nchunks, dimy, N, useMaxNorm)) {
    device_reset();
    return JIDT_ERROR;
  }


  // 6. Set locals or digammas for return
  // ======================

  float sumNx, sumNy, sumDiGammas, digammaK, digammaN;
  int result_size;

  if (nb_surrogates > 0) {
    printf("Surrogates not supported yet\n");
    return JIDT_ERROR;

  } else if (returnLocals) {
    digammaK = cpuDigamma(k);
    digammaN = cpuDigamma(N);
    for (i = 0; i < N; i++) {
      result[i] = digammaK + digammaN - cpuDigamma(nx[i] + 1) - cpuDigamma(ny[i] + 1);
    }

  } else {
    sumNx = 0;
    sumNy = 0;
    sumDiGammas = 0;
    for (i = 0; i < N; i++) {
      sumNx       += nx[i];
      sumNy       += ny[i];
      sumDiGammas += cpuDigamma(nx[i] + 1) + cpuDigamma(ny[i] + 1);
    }

    result[0] = (float) sumDiGammas;
    result[1] = (float) sumNx;
    result[2] = (float) sumNy;

  }

  return JIDT_SUCCESS;
}


#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     infodynamics_measures_continuous_kraskov_MutualInfoCalculatorMultiVariateKraskov
 * Method:    MIKraskov
 * Signature: (I[DI[DIIZZZI)[D
 */
JNIEXPORT jdoubleArray JNICALL
  Java_infodynamics_measures_continuous_kraskov_MutualInfoCalculatorMultiVariateKraskov_MIKraskov(
   JNIEnv *env, jobject thisObj, jint j_N,
   jobjectArray j_sourceArray, jint j_dimx,
   jobjectArray j_destArray,   jint j_dimy,
   jint j_k, jboolean j_returnLocals, jboolean j_useMaxNorm,
   jboolean j_isAlgorithm1, jint j_nbSurrogates) {

  // Declare variables
  // Variables starting with j_ are Java types
  float *source, *dest, *pointset, *radii, *distances;
  int *nx, *ny, *indexes;
  jdouble *jd_source, *jd_dest;
  unsigned int N, dimx, dimy, dims, k, returnLocals, useMaxNorm, isAlgorithm1, result_size;
  float sumNx, sumNy, sumDiGammas;
  unsigned int i, j;
  long double digammaK, digammaN;
  int thelier = 0, nchunks = 1, err, nb_surrogates;

  // Check that incoming data has correct size
  jsize sourceLength = (*env)->GetArrayLength(env, j_sourceArray);
  jsize destLength = (*env)->GetArrayLength(env, j_destArray);

  if (sourceLength != j_N || destLength != j_N) {
    device_reset();
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Data has wrong length.");
  }

  // Copy variables from Java
  // =====================
  N = j_N;
  k = j_k;
  dimx = j_dimx;
  dimy = j_dimy;
  dims = dimx + dimy;
  returnLocals = j_returnLocals ? 1 : 0;
  useMaxNorm   = j_useMaxNorm   ? 1 : 0;
  isAlgorithm1 = j_isAlgorithm1 ? 1 : 0;
  nb_surrogates = j_nbSurrogates;
  
  source   = (float *) malloc(N * dimx * sizeof(float));
  dest     = (float *) malloc(N * dimy * sizeof(float));

  if (NULL == source || NULL == dest) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Error allocating data.");
  }


  if ((returnLocals || !isAlgorithm1) && (nb_surrogates > 0)) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Surrogates only supported for average MI with KSG1.");
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
      source[N*j + i]   = (float) sourceRow[j];
    }

    for (j = 0; j < dimy; j++) {
      dest[N*j + i] = (float) destRow[j];
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


  // Call C function
  // =========================
  int resultSize;
  if (returnLocals) {
    resultSize = N;
  } else if (nb_surrogates > 0) {
    resultSize = nb_surrogates + 1;
  } else {
    resultSize = 3;
  }
  float *result = (float *) malloc(resultSize * sizeof(float));
  jidt_error_t ret = MIKraskov_C(N, source, dimx, dest, dimy,
                       k, thelier, nchunks, returnLocals, useMaxNorm,
                       isAlgorithm1, nb_surrogates, result);

  if (JIDT_ERROR == ret) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Surrogates only supported for average MI with KSG1.");
  }

  // Free memory and return
  // =========================
  if (source)  free(source);
  if (dest)    free(dest);

  jdouble outCArray[resultSize];
  for (int ii = 0; ii < resultSize; ii++) {
    outCArray[ii] = result[ii];
  }

  // Set Java array for return
  jdoubleArray outJNIArray = (*env)->NewDoubleArray(env, resultSize);  // allocate
  if (NULL == outJNIArray) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Error creating return array.");
  }
  (*env)->SetDoubleArrayRegion(env, outJNIArray, 0 , resultSize, outCArray);  // copy

  if (result) { free(result); }
  
  return outJNIArray;


} // End of function MIKraskov1


#ifdef __cplusplus
}
#endif

