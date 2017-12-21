#ifndef GPUMILIBRARY_H
#define GPUMILIBRARY_H

#include "gpuKnnLibrary.h"

#define FREE(x) { if (x) free(x); x = NULL; }

#ifdef __cplusplus
extern "C" {
#endif
jidt_error_t MIKraskovWithReorderings(int N, float *source, int dimx, float *dest, int dimy,
    int k, int thelier, int nb_surrogates, int returnLocals, int useMaxNorm,
    int isAlgorithm1, float *result, int reorderingsGiven, int **reorderings);

jidt_error_t MIKraskov_C(int N, float *source, int dimx, float *dest, int dimy,
    int k, int thelier, int nb_surrogates, int returnLocals, int useMaxNorm,
    int isAlgorithm1, float *result);

jidt_error_t MIKraskovByPointsetChunks(int N, float *source, int dimx,
    float *dest, int dimy, int k, int thelier, int nb_surrogates,
    int returnLocals, int useMaxNorm, int isAlgorithm1, float *result,
    float *pointset);
#ifdef __cplusplus
}
#endif

#endif

