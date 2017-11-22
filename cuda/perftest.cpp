#include <iostream>
#include <random>
#include <chrono>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include "gpuKnnLibrary.h"
#include "gpuMILibrary.h"
#include "digamma.h"
#include "ctimer.h"

int main(int argc, char *argv[]){
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0, 1.0);

  int N = 1000;
  int dimx = 1;
  int dimy = 1;
  int k = 4;
  int thelier = 0;
  int nb_surrogates = 0;
  int returnLocals = 0;
  int useMaxNorm = 1;
  int isAlgorithm1 = 1;
  int resultSize;

  if (argc > 1) {
    nb_surrogates = atoi(argv[1]);
  }

  CPerfTimer pt1 = startTimer("Time generating data");

  float *source = (float *) malloc(N*dimx*sizeof(float));
  float *dest = (float *) malloc(N*dimy*sizeof(float));

  for (int i = 0; i < N*dimx; ++i) {
    source[i] = distribution(generator);
  }
  for (int i = 0; i < N*dimy; ++i) {
    dest[i] = distribution(generator) + source[i];
  }

  if (nb_surrogates == 0) {
    resultSize = 3;
  } else {
    resultSize = nb_surrogates + 1;
  }
  float result[resultSize];
  stopTimer(pt1);

  pt1 = startTimer("MI_full");
  (void) MIKraskov_C(N, source, dimx, dest, dimy, k, thelier,
      nb_surrogates, returnLocals, useMaxNorm, isAlgorithm1, result);
  stopTimer(pt1);

  free(source);
  free(dest);

  return 0;
}

