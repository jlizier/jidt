####################################################
#
# Makefile for GPU Kraskov Mutual Information code
# inside the JIDT library.
#
# (C) Pedro Mediano, 2017
#
####################################################

## Common binaries and flags
#---------------------------
NVCC       ?= nvcc
GCC        ?= g++
CCFLAGS    ?= -O3 -Wall -Wextra -Wno-unused-parameter
NVCCFLAGS  ?= -O3 -D_FORCE_INLINES -Xcompiler -Wall


## Quick check to find 'which' function
#--------------------------------------
ifeq ($(OS),Windows_NT)
  which = $(shell where $1)
else
  which = $(shell which $1)
endif


## Check that NVCC is actually available
#---------------------------------------
ifeq ($(shell which ${NVCC}),)
  $(error No NVCC found. Check that CUDA is installed and that CUDA_PATH is set)
endif


## Find CUDA compute capability for target GPU
#---------------------------------------------
ifneq ($(MAKECMDGOALS),clean)
ifeq (${COMPUTE_CAPABILITY},)
  ifeq ($(wildcard findComputeCapability),)
    $(info Compiling diagnostics script)
    $(shell ${NVCC} -Wno-deprecated-gpu-targets -o findComputeCapability findComputeCapability.cu)
  endif
  COMPUTE_CAPABILITY ?= $(shell ./findComputeCapability)
  $(info CUDA compute capability found: ${COMPUTE_CAPABILITY})
endif
endif
NVCCFLAGS += -arch=sm_${COMPUTE_CAPABILITY}


## Find path to JNI header files
#-------------------------------
# Set these variables if you know where your JNI header files are
# JNI_INCLUDE ?=
# JNI_PLATFORM_INCLUDE ?=
ifeq ($(and $(JNI_INCLUDE),$(JNI_PLATFORM_INCLUDE)),)
  include FindJNI.mk
endif


## Set more flags depending on make arguments
#--------------------------------------------
ifdef DEBUG
  CCFLAGS    += -g
  NVCCFLAGS  += -g -G
endif

ifdef CUSTOMFLAGS
  CCFLAGS    += ${CUSTOMFLAGS}
  NVCCFLAGS  += ${CUSTOMFLAGS}
endif

# Common includes and paths for CUDA. This assumes the CUDA toolkit is in PATH
INCLUDES      := -I. -I./cub -I${JNI_INCLUDE} -I${JNI_PLATFORM_INCLUDE}
LDFLAGS       += -L. -lcuda -lcudart

# Main target rule
all:  libKraskov.so


## Compile device code
#---------------------
gpuKnnLibrary.o: gpuKnnLibrary.c helperfunctions.cu gpuKnnBF_kernel.cu
	${NVCC} ${NVCCFLAGS} ${INCLUDES} -x cu -Xcompiler -fPIC -c gpuKnnLibrary.c

libgpuKnnLibrary.a: gpuKnnLibrary.o
	${AR} -r libgpuKnnLibrary.a gpuKnnLibrary.o


## Compile host code
#-------------------
c_objects = digamma.o gpuMILibrary.o kraskovCuda.o

%.o:	%.c
	${GCC} ${INCLUDES} ${CCFLAGS} -x c -std=c99 -fPIC -c $< -o $@


## Final shared library linking
#------------------------------
libKraskov.so: libgpuKnnLibrary.a $(c_objects)
	${NVCC} ${NVCCFLAGS} ${INCLUDES} -Xcompiler -fPIC -shared -o libKraskov.so $(c_objects) ${LDFLAGS} -lgpuKnnLibrary


## Test binary targets
#------------------------------
test:	unittest perftest

unittest:	libKraskov.so unittest.cpp
	${GCC} ${CCFLAGS} -std=c++11 unittest.cpp -L. -lKraskov -o $@ 

perftest: libKraskov.so perftest.cpp 
	${GCC} ${CCFLAGS} -DTIMER -std=c++11 perftest.cpp -L. -lKraskov -o $@ 

clean:
	rm -f *.o *.so *.a unittest perftest

