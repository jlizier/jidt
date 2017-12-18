#ifndef GPUKNNLIBRARY_H
#define GPUKNNLIBRARY_H

#ifdef __cplusplus
extern "C" {
#endif
typedef enum { JIDT_SUCCESS, JIDT_ERROR } jidt_error_t;

int allocateDeviceMemory(int signalLength, int kth, int dimx, int dimy,
    float **source, float **dest, float **distances, int **indexes,
    float **radii, int **nx, int **ny, float **digammas, float *pointset);

int allocateDeviceMemoryCMI(int signalLength, int k, int dimx, int dimy, int dimz,
    float **source, float **dest, float **cond, float **distances, int **indexes,
    float **radii, int **nx, int **ny, int **nz, float **digammas, float *pointset);

int freeDeviceMemory(float *d_pointset);

int cudaFindKnn(int* h_bf_indexes, float* h_bf_distances, float* h_pointset,
    float* h_query, int kth, int thelier, int nchunks, int pointdim,
    int signallength, unsigned int useMaxNorm);

int cudaFindKnnSetGPU(int* h_bf_indexes, float* h_bf_distances,
    float* h_pointset, float* h_query, int kth, int thelier, int nchunks,
    int pointdim, int signallength, unsigned int useMaxNorm, int deviceid);

int cudaFindRSAll(int* h_bf_npointsrange, float* h_pointset,
    float* h_query, float* h_vecradius, int thelier, int nchunks,
    int pointdim, int signallength, unsigned int useMaxNorm);

int cudaFindRSAllSetGPU(int* h_bf_npointsrange, float* h_pointset,
    float* h_query, float* h_vecradius, int thelier, int nchunks,
    int pointdim, int signallength, unsigned int useMaxNorm, int deviceid);

int findRadiiAlgorithm2(float *radii, const float *data, const int *indexes,
    unsigned int k, unsigned int dim, unsigned int N);

int computeSumDigammas(float *sumDiGammas, int *nx, int *ny, unsigned int N);

int parallelDigammas(float *digammas, int *nx, int *ny, int signallength);

int cudaSumDigammas(float *sumDigammas, int *nx, int *ny, int trialLength,
    int nchunks);

int d_cudaFindKnn(int* d_bf_indexes, float* d_bf_distances, float* d_bf_pointset,
    float* d_bf_query, int kth, int thelier, int nchunks, int pointdim,
    int signallength, int useMaxNorm);

int d_cudaFindRSAll(int* d_bf_npointsrange, float* d_bf_pointset, float* d_bf_query,
    float* d_bf_vecradius, int thelier, int nchunks, int pointdim,
    int signallength, int useMaxNorm);

int d_parallelDigammas(float *digammas, float *d_digammas, int *d_nx,
    int *d_ny, int signalLength);

int d_parallelDigammasCMI(float *digammas, float *d_digammas, int *d_nx,
    int *d_ny, int *d_nz, int signalLength);

int cudaBlockReduce(float *sumDigammas, float *d_digammas, int trialLength, int nchunks);

int d_cudaSumDigammas(float *sumDigammas, int *d_nx, int *d_ny,
    float *d_digammas, int trialLength, int nchunks);

int d_cudaSumDigammasCMI(float *sumDigammas, int *d_nx, int *d_ny, int *d_nz,
    float *d_digammas, int trialLength, int nchunks);

void device_reset(void);

void gpuWarmUp(void);

void randperm(int perm[], int n);
#ifdef __cplusplus
}
#endif

#endif

