#include <jni.h>
#include <stdlib.h>
#include <stdio.h>

#define check(ans) { _check((ans), __FILE__, __LINE__); }

#include "gpuMILibrary.h"
#include "gpuCMILibrary.h"
#include "ctimer.h"

#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     infodynamics_measures_continuous_kraskov_MutualInfoCalculatorMultiVariateKraskov
 * Method:    MIKraskov
 * Signature: (I[DI[DIIIZZZIZ[I)[D
 */
JNIEXPORT jdoubleArray JNICALL
  Java_infodynamics_measures_continuous_kraskov_MutualInfoCalculatorMultiVariateKraskov_MIKraskov(
   JNIEnv *env, jobject thisObj, jint j_N,
   jobjectArray j_sourceArray, jint j_dimx,
   jobjectArray j_destArray,   jint j_dimy,
   jint j_k, jint j_theiler, jboolean j_returnLocals,
   jboolean j_useMaxNorm, jboolean j_isAlgorithm1, jint j_nbSurrogates,
   jboolean j_reorderingsGiven, jobjectArray j_orderings) {

  // Check that incoming data has correct size
  // =====================
  jsize sourceLength = (*env)->GetArrayLength(env, j_sourceArray);
  jsize destLength = (*env)->GetArrayLength(env, j_destArray);

  // if (sourceLength != j_N || destLength != j_N || (j_N%(j_nbSurrogates+1) != 0)) {
  if (sourceLength != j_N || destLength != j_N) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Data has wrong length.");
  }

  if (!j_isAlgorithm1) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Only algorithm 1 is supported.");
  }

  if ((j_returnLocals || !j_isAlgorithm1) && (j_nbSurrogates > 0)) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Surrogates only supported for average MI with KSG1.");
  }

  // Copy variables from Java
  // =====================
  int N = j_N;
  int k = j_k;
  int dimx = j_dimx;
  int dimy = j_dimy;
  int theiler = j_theiler;
  int returnLocals = j_returnLocals ? 1 : 0;
  int useMaxNorm   = j_useMaxNorm   ? 1 : 0;
  int isAlgorithm1 = j_isAlgorithm1 ? 1 : 0;
  int nb_surrogates = j_nbSurrogates;
  int reorderingsGiven = j_reorderingsGiven ? 1 : 0;
   
  CPerfTimer pt = startTimer("Java array copy");

  float *source = (float *) malloc(N * dimx * sizeof(float));
  float *dest   = (float *) malloc(N * dimy * sizeof(float));

  if (NULL == source || NULL == dest) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Error allocating data.");
  }

  for (int i = 0; i < N; i++) {
    jdoubleArray j_sourceRow = (jdoubleArray) (*env)->GetObjectArrayElement(env, j_sourceArray, i);
    jdoubleArray j_destRow   = (jdoubleArray) (*env)->GetObjectArrayElement(env, j_destArray, i);

    jdouble *sourceRow = (*env)->GetDoubleArrayElements(env, j_sourceRow, NULL);
    jdouble *destRow = (*env)->GetDoubleArrayElements(env, j_destRow, NULL);

    // Data in java are doubles, but GPUs need floats.
    // We have to cast them manually, so we can't memcopy

    // The following two for-loops get two matrices in T-by-D indexing (i.e.
    // first dimension is time, second is variable) and return the data in
    // column-major form
    for (int j = 0; j < dimx; j++) {
      source[N*j + i]   = (float) sourceRow[j];
    }

    for (int j = 0; j < dimy; j++) {
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

  int **reorderings = NULL;
  if (reorderingsGiven) {
    reorderings = (int **) malloc(nb_surrogates * sizeof(int *));
    for (int i = 0; i < nb_surrogates; i++) {
      jintArray j_order = (jdoubleArray) (*env)->GetObjectArrayElement(env, j_orderings, i);
      jint *order = (*env)->GetIntArrayElements(env, j_order, NULL);

      reorderings[i] = (int *) malloc(N * sizeof(int));
      for (int j = 0; j < N; j++) {
        reorderings[i][j] = order[j];
      }

      (*env)->ReleaseIntArrayElements(env, j_order, order, 0);
      (*env)->DeleteLocalRef(env, j_order);
    }
  }

  stopTimer(pt);


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
  jidt_error_t ret;
  if (!reorderingsGiven) {
    ret = MIKraskov_C(N, source, dimx, dest, dimy,
                      k, theiler, nb_surrogates, returnLocals, useMaxNorm,
                      isAlgorithm1, result);

  } else {
    ret = MIKraskovWithReorderings(N, source, dimx, dest, dimy, k, theiler,
                                   nb_surrogates, returnLocals, useMaxNorm,
                                   isAlgorithm1, result, reorderingsGiven,
                                   reorderings);
  }

  if (JIDT_ERROR == ret) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Error in GPU execution.");
  }

  // Free memory and return
  // =========================
  if (source)  free(source);
  if (dest)    free(dest);

  if (reorderingsGiven) {
    for (int i = 0; i < nb_surrogates; i++) {
      if (reorderings[i]) free(reorderings[i]);
    }
    if (reorderings) free (reorderings);
  }

  jdouble outCArray[resultSize];
  for (int i = 0; i < resultSize; i++) {
    outCArray[i] = result[i];
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


} // End of function MIKraskov


/*
 * Class:     infodynamics_measures_continuous_kraskov_ConditionalMutualInfoCalculatorMultiVariateKraskov
 * Method:    CMIKraskov
 * Signature: (I[DI[DI[DIIIZZZIZ[II)[D
 */
JNIEXPORT jdoubleArray JNICALL
  Java_infodynamics_measures_continuous_kraskov_ConditionalMutualInfoCalculatorMultiVariateKraskov_CMIKraskov(
   JNIEnv *env, jobject thisObj, jint j_N,
   jobjectArray j_sourceArray, jint j_dimx,
   jobjectArray j_destArray,   jint j_dimy,
   jobjectArray j_condArray,   jint j_dimz,
   jint j_k, jint j_theiler, jboolean j_returnLocals,
   jboolean j_useMaxNorm, jboolean j_isAlgorithm1, jint j_nbSurrogates,
   jboolean j_reorderingsGiven, jobjectArray j_orderings,
   jint j_variableToReorder) {

  // Check that incoming data has correct size
  // =====================
  jsize sourceLength = (*env)->GetArrayLength(env, j_sourceArray);
  jsize destLength = (*env)->GetArrayLength(env, j_destArray);
  jsize condLength = (*env)->GetArrayLength(env, j_condArray);

  // if (sourceLength != j_N || destLength != j_N || (j_N%(j_nbSurrogates+1) != 0)) {
  if (sourceLength != j_N || destLength != j_N || condLength != j_N) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Data has wrong length.");
  }

  if (!j_isAlgorithm1) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Only algorithm 1 is supported.");
  }

  if ((j_returnLocals || !j_isAlgorithm1) && (j_nbSurrogates > 0)) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Surrogates only supported for average MI with KSG1.");
  }

  // Copy variables from Java
  // =====================
  int N = j_N;
  int k = j_k;
  int dimx = j_dimx;
  int dimy = j_dimy;
  int dimz = j_dimz;
  int theiler = j_theiler;
  int returnLocals = j_returnLocals ? 1 : 0;
  int useMaxNorm   = j_useMaxNorm   ? 1 : 0;
  int isAlgorithm1 = j_isAlgorithm1 ? 1 : 0;
  int nb_surrogates = j_nbSurrogates;
  int reorderingsGiven = j_reorderingsGiven ? 1 : 0;
  int variableToReorder = j_variableToReorder;
   
  CPerfTimer pt = startTimer("Java array copy");

  float *source = (float *) malloc(N * dimx * sizeof(float));
  float *dest   = (float *) malloc(N * dimy * sizeof(float));
  float *cond   = (float *) malloc(N * dimz * sizeof(float));

  if (NULL == source || NULL == dest || NULL == cond) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Error allocating data.");
  }

  for (int i = 0; i < N; i++) {
    jdoubleArray j_sourceRow = (jdoubleArray) (*env)->GetObjectArrayElement(env, j_sourceArray, i);
    jdoubleArray j_destRow   = (jdoubleArray) (*env)->GetObjectArrayElement(env, j_destArray, i);
    jdoubleArray j_condRow   = (jdoubleArray) (*env)->GetObjectArrayElement(env, j_condArray, i);

    jdouble *sourceRow = (*env)->GetDoubleArrayElements(env, j_sourceRow, NULL);
    jdouble *destRow = (*env)->GetDoubleArrayElements(env, j_destRow, NULL);
    jdouble *condRow = (*env)->GetDoubleArrayElements(env, j_condRow, NULL);

    // Data in java are doubles, but GPUs need floats.
    // We have to cast them manually, so we can't memcopy

    // The following for-loops get two matrices in T-by-D indexing (i.e.
    // first dimension is time, second is variable) and return the data in
    // column-major form
    for (int j = 0; j < dimx; j++) {
      source[N*j + i]   = (float) sourceRow[j];
    }

    for (int j = 0; j < dimy; j++) {
      dest[N*j + i] = (float) destRow[j];
    }

    for (int j = 0; j < dimz; j++) {
      cond[N*j + i] = (float) condRow[j];
    }

    (*env)->ReleaseDoubleArrayElements(env, j_sourceRow, sourceRow, 0);
    (*env)->ReleaseDoubleArrayElements(env, j_destRow, destRow, 0);
    (*env)->ReleaseDoubleArrayElements(env, j_condRow, condRow, 0);

    (*env)->DeleteLocalRef(env, j_sourceRow);
    (*env)->DeleteLocalRef(env, j_destRow);
    (*env)->DeleteLocalRef(env, j_condRow);

    // FIXME: I'm not entirely sure I'm freeing all the memory here. I should
    // check for memory leaks more carefully.

  }

  int **reorderings = NULL;
  if (reorderingsGiven) {
    reorderings = (int **) malloc(nb_surrogates * sizeof(int *));
    for (int i = 0; i < nb_surrogates; i++) {
      jintArray j_order = (jdoubleArray) (*env)->GetObjectArrayElement(env, j_orderings, i);
      jint *order = (*env)->GetIntArrayElements(env, j_order, NULL);

      reorderings[i] = (int *) malloc(N * sizeof(int));
      for (int j = 0; j < N; j++) {
        reorderings[i][j] = order[j];
      }

      (*env)->ReleaseIntArrayElements(env, j_order, order, 0);
      (*env)->DeleteLocalRef(env, j_order);
    }
  }

  stopTimer(pt);


  // Call C function
  // =========================
  int resultSize;
  if (returnLocals) {
    resultSize = N;
  } else if (nb_surrogates > 0) {
    resultSize = nb_surrogates + 1;
  } else {
    resultSize = 6;
  }
  float *result = (float *) malloc(resultSize * sizeof(float));
  jidt_error_t ret;
  if (!reorderingsGiven) {
    ret = CMIKraskov_C(N, source, dimx, dest, dimy, cond, dimz,
                      k, theiler, nb_surrogates, returnLocals, useMaxNorm,
                      isAlgorithm1, result, variableToReorder);

  } else {
    ret = CMIKraskovWithReorderings(N, source, dimx, dest, dimy, cond, dimz, k, theiler,
                                   nb_surrogates, returnLocals, useMaxNorm,
                                   isAlgorithm1, result, reorderingsGiven,
                                   reorderings, variableToReorder);
  }

  if (JIDT_ERROR == ret) {
    jclass Exception = (*env)->FindClass(env, "java/lang/Exception");
    (*env)->ThrowNew(env, Exception, "Error in GPU execution.");
  }

  // Free memory and return
  // =========================
  if (source)  free(source);
  if (dest)    free(dest);
  if (cond)    free(cond);

  if (reorderingsGiven) {
    for (int i = 0; i < nb_surrogates; i++) {
      if (reorderings[i]) free(reorderings[i]);
    }
    if (reorderings) free (reorderings);
  }

  jdouble outCArray[resultSize];
  for (int i = 0; i < resultSize; i++) {
    outCArray[i] = result[i];
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


} // End of function CMIKraskov


#ifdef __cplusplus
}
#endif

