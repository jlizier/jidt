#ifndef _TEMPLATE_KERNEL_H_
#define _TEMPLATE_KERNEL_H_

#include <stdio.h>

//#include <cpugpuKnn_common.h>
#ifndef INFINITY
#define INFINITY 0x7F800000
#endif

/**
 * Typedef function pointer for arbitrary norm functions.
 */
typedef float (*normFunction_t)(const float* g_uquery, const float* g_vpoint,
                  int pointdim, int signallength);

/**
 * Calculate max norm (L-inf) between two points.
 */
__device__ __host__ float
maxMetricPoints(const float* g_uquery, const float* g_vpoint, int pointdim, int signallength){
  float r_u1;
  float r_v1;
  float r_d1,r_dim=0;

  r_dim=0;
  for(int d=0; d<pointdim; d++){
    r_u1 = *(g_uquery+d*signallength);
    r_v1 = *(g_vpoint+d*signallength);
    r_d1 = fabsf(r_v1 - r_u1);
    r_dim= r_dim < r_d1? r_d1: r_dim;
  }
  return r_dim;
}

/**
 * Calculate squared Euclidean norm (L2) between two points (note we keep the
 * distance squared to avoid taking sqrt all the time).
 */
__device__ __host__ float
squareMetricPoints(const float* g_uquery, const float* g_vpoint, int pointdim, int signallength){
  float r_u1;
  float r_v1;
  float r_d1, r_dim = 0.0f;

  for (int d = 0; d < pointdim; d++) {
    r_u1  = *(g_uquery+d*signallength);
    r_v1  = *(g_vpoint+d*signallength);
    r_d1  = r_v1 - r_u1;
    r_dim += r_d1 * r_d1;
  }
  return r_dim;
}

__device__ normFunction_t pMaxNorm = maxMetricPoints;
__device__ normFunction_t pSquareNorm = squareMetricPoints;

/**
 * Insert point in current list of nearest neighbours.
 */
__device__ float
insertPointKlist(int kth, float distance, int indexv,float* kdistances, int* kindexes){
  int k=0;
  while( (distance>*(kdistances+k)) && (k<kth-1)){k++;}
  //Move value to the next
  for(int k2=kth-1;k2>k;k2--){
    *(kdistances+k2)=*(kdistances+k2-1);
    *(kindexes+k2)=*(kindexes+k2-1);
  }
  //Replace
  *(kdistances+k)=distance;
  *(kindexes+k)=indexv;

  //printf("\n -> Modificacion pila: %.f %.f. New max distance: %.f", *kdistances, *(kdistances+1), *(kdistances+kth-1));
  return *(kdistances+kth-1);
}


/*
 * Main KNN kernel. Find the nearest k neighbours to each point according
 * to the supplied norm function.
 */
__global__ void
kernelKNNshared(const float* g_uquery, const float* g_vpointset,
    int *g_indexes, float* g_distances, const int pointdim,
    const int triallength, const int signallength, const int kth,
    const int exclude, normFunction_t *normFunction) {

  // Shared memory
  extern __shared__ char array[];
  float *kdistances;
  int *kindexes;
  kdistances = (float*)array;
  kindexes = (int*)array+kth*blockDim.x;

  const unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
  const unsigned int itrial = tid / triallength;  //  indextrial

if(tid<signallength){

  for(int k=0;k<kth;k++){
    kdistances[threadIdx.x*kth+k] = INFINITY;
  }

  __syncthreads();

  float r_kdist=INFINITY;
  unsigned int indexi = tid-triallength*itrial;

  for(int t=0; t<triallength; t++){
      int indexu = tid;
      int indexv = (t + itrial*triallength);
      int condition1=indexi-exclude;
      int condition2=indexi+exclude;
      if((t<condition1)||(t>condition2)){
        float temp_dist = normFunction[0](g_uquery+indexu, g_vpointset+indexv,pointdim, signallength);
        if(temp_dist <= r_kdist){
          r_kdist = insertPointKlist(kth,temp_dist,t,kdistances+threadIdx.x*kth,kindexes+threadIdx.x*kth);
        }
      }
      //printf("tid:%d indexes: %d, %d distances: %.f %.f\n",tid, *kindexes, *(kindexes+1), *kdistances, *(kdistances+1));
  }

  __syncthreads();
  // COPY TO GLOBAL MEMORY
  for (int k = 0; k < kth; k++) {
    g_indexes[tid+k*signallength] = kindexes[threadIdx.x*kth+k];
    g_distances[tid+k*signallength] = kdistances[threadIdx.x*kth+k];//*(kdistances+k);
  }
}

}


/*
 * Range search for one data point using bruteforce.
 */
__global__ void
kernelBFRSshared(const float* g_uquery, const float* g_vpointset, int *g_npoints, int pointdim, int triallength, int signallength, int exclude, float radius)
{

  // Shared memory
  extern __shared__ char array[];
  int *s_npointsrange;
  s_npointsrange = (int*)array;

  const unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
  const unsigned int itrial = tid / triallength;  //  indextrial

if(tid<signallength){

  s_npointsrange[threadIdx.x] = 0;
  __syncthreads();


  unsigned int indexi = tid-triallength*itrial;
  for(int t=0; t<triallength; t++){
      int indexu = tid;
      int indexv = (t + itrial*triallength);
      int condition1=indexi-exclude;
      int condition2=indexi+exclude;
      if((t<condition1)||(t>condition2)){
        float temp_dist = maxMetricPoints(g_uquery+indexu, g_vpointset+indexv,pointdim, signallength);
        if(temp_dist <= radius){
          s_npointsrange[threadIdx.x]++;
        }
      }

  }

  __syncthreads();

  // COPY TO GLOBAL MEMORY
  g_npoints[tid] = s_npointsrange[threadIdx.x];

}
}

/*
 * Range search using bruteforce in multiple GPUs
 */
__global__ void
kernelBFRSMultishared(const float* g_uquery, const float* g_vpointset, int *g_npoints, int pointdim, int triallength, int signallength, int exclude, const float* vecradius)
{

    // shared memory
  extern __shared__ char array[];
  int *s_npointsrange;
  s_npointsrange = (int*)array;
    float radius=0;
  const unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
  const unsigned int itrial = tid / triallength;  //  indextrial

if(tid<signallength){

  s_npointsrange[threadIdx.x] = 0;
  __syncthreads();

    radius = *(vecradius+itrial);
  unsigned int indexi = tid-triallength*itrial;
  for(int t=0; t<triallength; t++){
      int indexu = tid;
      int indexv = (t + itrial*triallength);
      int condition1=indexi-exclude;
      int condition2=indexi+exclude;
      if((t<condition1)||(t>condition2)){
        float temp_dist = maxMetricPoints(g_uquery+indexu, g_vpointset+indexv,pointdim, signallength);
        if(temp_dist <= radius){
          s_npointsrange[threadIdx.x]++;
        }
      }

  }

  __syncthreads();
  //printf("\ntid:%d npoints: %d\n",tid, s_npointsrange[threadIdx.x]);
  //COPY TO GLOBAL MEMORY
  g_npoints[tid] = s_npointsrange[threadIdx.x];

}
}


/*
 * Range search for all data points using bruteforce.
 */
__global__ void
kernelBFRSAllshared(const float* g_uquery, const float* g_vpointset,
    int *g_npoints, int pointdim, int triallength, int signallength,
    int exclude, const float* vecradius, normFunction_t *normFunction) {

  // Shared memory
  extern __shared__ char array[];
  int *s_npointsrange;
  s_npointsrange = (int *) array;
  float radius = 0;
  const unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
  const unsigned int itrial = tid / triallength;  //  indextrial

  if(tid<signallength){

    s_npointsrange[threadIdx.x] = 0;
    __syncthreads();

    radius = *(vecradius+tid);
    unsigned int indexi = tid-triallength*itrial;
    for (int t=0; t<triallength; t++){
        // Note: the following two definitions could be swapped depending on
        // the details of surrogate implementation
        int indexu = tid; // old, necessary for shuffled part of the surrogates (i.e. source)
        // int indexu = indexi; // new, admissible for unshuffled part of the surrogates (i.e. dest)
        int indexv = (t + itrial*triallength);
        int condition1=indexi-exclude;
        int condition2=indexi+exclude;
        if((t<condition1)||(t>condition2)){
          float temp_dist = normFunction[0](g_uquery+indexu, g_vpointset+indexv,pointdim, signallength);
          // PEDRO: For KSG algorithm 1 this should be strictly less than R, and in TRENTOOL code it's less or equal. It's a float comparison, so I don't think it matters anyway.
          if(temp_dist < radius){
            s_npointsrange[threadIdx.x]++;
          }
        }

    }

    __syncthreads();

    //COPY TO GLOBAL MEMORY
    g_npoints[tid] = s_npointsrange[threadIdx.x];

  }
}


/**
 * Taken from NVIDIA dev forum:
 * https://devtalk.nvidia.com/default/topic/516516/kernel-launch-failure-in-matlab/
 */
__device__ void digammaXp1(double *y) {

    // double x = *y;
    double x = *y + 1;
    double neginf = -INFINITY;
    const double c = 12,
    digamma1 = -0.57721566490153286,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 */
    s = 1e-6,
    s3 = 1./12,
    s4 = 1./120,
    s5 = 1./252,
    s6 = 1./240,
    s7 = 1./132;
    // s8 = 691./32760,
    // s9 = 1./12,
    // s10 = 3617./8160;
  double result;

  /* Illegal arguments */
  if((x == neginf) || isnan(x)) {
    *y = NAN;
    return;
  }

  /* Singularities */
  if((x <= 0) && (floor(x) == x)) {
    *y = neginf;
    return;
  }

  /* Negative values */

  /* Use the reflection formula (Jeffrey 11.1.6):
   * digamma(-x) = digamma(x+1) + pi*cot(pi*x)
   *
   * This is related to the identity
   * digamma(-x) = digamma(x+1) - digamma(z) + digamma(1-z)
   * where z is the fractional part of x
   * For example:

   * digamma(-3.1) = 1/3.1 + 1/2.1 + 1/1.1 + 1/0.1 + digamma(1-0.1)
   *               = digamma(4.1) - digamma(0.1) + digamma(1-0.1)
   * Then we use
   * digamma(1-z) - digamma(z) = pi*cot(pi*z)
   *

  if(x < 0) {
    *p = digamma(p,1-x) + M_PI/tan(-M_PI*x);
    return;
  }
  */

  /* Use Taylor series if argument <= S */
  if(x <= s) {
      *y = digamma1 - 1/x + trigamma1*x;
      return;
  }

  /* Reduce to digamma(X + N) where (X + N) >= C */
  result = 0;
  while(x < c) {
    result -= 1/x;
    x++;
  }

  /* Use de Moivre's expansion if argument >= C */
  /* This expansion can be computed in Maple via asympt(Psi(x),x) */
  if(x >= c) {
    double r = 1/x, t;
    result += log(x) - 0.5*r;
    r *= r;
#if 0
    result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
#else
    /* this version for lame compilers */
    t = (s5 - r * (s6 - r * s7));
    result -= r * (s3 - r * (s4 - r * t));
#endif
  }

  /* assign the result to the pointer*/
  *y = result;
  return;

}



/**
 * Optimized reduction kernel, obtained from Mark Harris' lecture on
 * GPU optimization:
 * https://docs.nvidia.com/cuda/samples/6_Advanced/reduction/doc/reduction.pdf
 */
template <unsigned int blockSize>
// __device__ void warpReduce(volatile int *sdata, unsigned int tid) {
__device__ void warpReduce(volatile float *sdata, unsigned int tid) {
  if (blockSize >=  64) sdata[tid] += sdata[tid + 32];
  if (blockSize >=  32) sdata[tid] += sdata[tid + 16];
  if (blockSize >=  16) sdata[tid] += sdata[tid +  8];
  if (blockSize >=   8) sdata[tid] += sdata[tid +  4];
  if (blockSize >=   4) sdata[tid] += sdata[tid +  2];
  if (blockSize >=   2) sdata[tid] += sdata[tid +  1];
}

// template <unsigned int blockSize>
// __global__ void reduce6(int *g_idata, int *g_odata, unsigned int n) {
__global__ void reduce6(int *g_nx, int *g_ny, float *g_odata, unsigned int n) {
// extern __shared__ int sdata[];
extern __shared__ float sdata[];
unsigned int tid = threadIdx.x;
const unsigned int blockSize = 512;
unsigned int i = blockIdx.x*(blockSize*2) + tid;
unsigned int gridSize = blockSize*2*gridDim.x;
sdata[tid] = 0;

while (i<n) {
  double dgX1 = (double) g_nx[i];
  double dgY1 = (double) g_ny[i];
  double dgX2 = (double) g_nx[i+blockSize];
  double dgY2 = (double) g_ny[i+blockSize];

  digammaXp1(&dgX1);
  digammaXp1(&dgY1);
  digammaXp1(&dgX2);
  digammaXp1(&dgY2);

  sdata[tid] = dgX1 + dgY1 + dgX2 + dgY2;
  i += gridSize;
}
__syncthreads();

if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
if (tid < 32) warpReduce<blockSize>(sdata, tid);
if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}


__global__ void reduce1(int *g_nx, int *g_ny, float *g_odata, unsigned int N) {
  extern __shared__ float sdata2[];
  // each thread loads one element from global to shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

  if (i < N) {

    double dgX = (double) g_nx[i];
    double dgY = (double) g_ny[i];

    digammaXp1(&dgX);
    digammaXp1(&dgY);

    sdata2[tid] = (float) (dgX + dgY);

  } else {
    sdata2[tid] = 0;
  }
  __syncthreads();

  // do reduction in shared mem
  for (unsigned int s=1; s < blockDim.x; s *= 2) {
    if (tid % (2*s) == 0) {
      sdata2[tid] += sdata2[tid + s];
    }
    __syncthreads();
  }
  // write result for this block to global mem
  if (tid == 0) {
    g_odata[blockIdx.x] = sdata2[0];
  }
}


__global__ void gpuDigammas(float *g_digammas, int *g_nx, int *g_ny, int signallength) {
  const unsigned int i = threadIdx.x + blockDim.x*blockIdx.x;

  if(i < signallength){
    // Fetch n and put it in thread memory
    double dgX = (double) g_nx[i];
    double dgY = (double) g_ny[i];

    // In-place digamma calculation
    digammaXp1(&dgX);
    digammaXp1(&dgY);

    // Copy back to global memory
    g_digammas[i] = (float) (dgX + dgY);
  }
  return;
}


__global__ void gpuDigammasCMI(float *g_digammas, int *g_nx, int *g_ny, int *g_nz, int signallength) {
  const unsigned int i = threadIdx.x + blockDim.x*blockIdx.x;

  if(i < signallength){
    // Fetch n and put it in thread memory
    double dgX = (double) g_nx[i];
    double dgY = (double) g_ny[i];
    double dgZ = (double) g_nz[i];

    // In-place digamma calculation
    digammaXp1(&dgX);
    digammaXp1(&dgY);
    digammaXp1(&dgZ);

    // Copy back to global memory
    g_digammas[i] = (float) (dgX + dgY - dgZ);
  }
  return;
}








#endif // #ifndef _TEMPLATE_KERNEL_H_


