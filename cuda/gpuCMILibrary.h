#ifndef GPUCMILIBRARY_H
#define GPUCMILIBRARY_H

#include "gpuKnnLibrary.h"

#define FREE(x) { if (x) free(x); x = NULL; }

#ifdef __cplusplus
extern "C" {
#endif
jidt_error_t CMIKraskovWithReorderings(int N, float *source, int dimx,
    float *dest, int dimy, float *cond, int dimz,
    int k, int thelier, int nb_surrogates, int returnLocals, int useMaxNorm,
    int isAlgorithm1, float *result, int reorderingsGiven, int **reorderings,
    int variableToReorder);

jidt_error_t CMIKraskov_C(int N, float *source, int dimx, float *dest, int dimy,
    float *cond, int dimz, int k, int thelier, int nb_surrogates,
    int returnLocals, int useMaxNorm, int isAlgorithm1, float *result,
    int variableToReorder);

jidt_error_t CMIKraskovByPointsetChunks(int N, float *source, int dimx,
    float *dest, int dimy, float *cond, int dimz, int k, int thelier, int nb_surrogates,
    int returnLocals, int useMaxNorm, int isAlgorithm1, float *result, float *pointset);
#ifdef __cplusplus
}
#endif

#endif

