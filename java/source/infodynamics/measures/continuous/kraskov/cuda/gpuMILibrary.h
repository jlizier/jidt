#ifndef GPUMILIBRARY_H
#define GPUMILIBRARY_H

#define FREE(x) { if (x) free(x); x = NULL; }

typedef enum { JIDT_SUCCESS, JIDT_ERROR } jidt_error_t;

jidt_error_t MIKraskov_C(int N, float *source, int dimx, float *dest, int dimy,
    int k, int thelier, int nchunks, int returnLocals, int useMaxNorm,
    int isAlgorithm1, float *result);

#endif

