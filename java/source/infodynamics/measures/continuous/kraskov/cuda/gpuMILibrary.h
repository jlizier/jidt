#ifndef GPUMILIBRARY_H
#define GPUMILIBRARY_H

#define FREE(x) { if (x) free(x); x = NULL; }

#ifdef __cplusplus
extern "C" {
#endif
typedef enum { JIDT_SUCCESS, JIDT_ERROR } jidt_error_t;

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

void randperm(int perm[], int n);
#ifdef __cplusplus
}
#endif

#endif

