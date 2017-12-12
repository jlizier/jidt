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
BIN         =../bin/cuda
UNITBIN     =../unittests/bin/cuda


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
# Set this variable if you know your target compute capability(ies). Example:
# COMPUTE_CAPABILITY = -arch=sm_30
ifneq ($(MAKECMDGOALS),clean)
ifeq (${COMPUTE_CAPABILITY},)
  ifeq ($(wildcard $(BIN)/findComputeCapability),)
    $(info Compiling diagnostics script)
    $(shell ${NVCC} -Wno-deprecated-gpu-targets -o $(BIN)/findComputeCapability findComputeCapability.cu)
  endif
  COMPUTE_CAPABILITY ?= $(shell $(BIN)/findComputeCapability)
  $(info CUDA compute capability found: ${COMPUTE_CAPABILITY})
endif
endif
NVCCFLAGS += ${COMPUTE_CAPABILITY}


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
LDFLAGS       += -L. -L$(BIN) -lcuda -lcudart

.PHONY: all clean test

# Main target rule
all:  $(BIN)/libKraskov.so


## Compile device code
#---------------------
$(BIN)/gpuKnnLibrary.o: gpuKnnLibrary.c helperfunctions.cu gpuKnnBF_kernel.cu
	${NVCC} ${NVCCFLAGS} ${INCLUDES} -x cu -Xcompiler -fPIC -c gpuKnnLibrary.c -o $@

$(BIN)/libgpuKnnLibrary.a: $(BIN)/gpuKnnLibrary.o
	${AR} -r $@ $<


## Compile host code
#-------------------
c_objects = $(addprefix $(BIN)/,digamma.o gpuMILibrary.o gpuCMILibrary.o kraskovCuda.o)

$(BIN)/%.o:	%.c
	${GCC} ${INCLUDES} ${CCFLAGS} -x c -std=c99 -fPIC -c $< -o $@


## Final shared library linking
#------------------------------
$(BIN)/libKraskov.so: $(BIN)/libgpuKnnLibrary.a $(c_objects)
	${NVCC} ${NVCCFLAGS} ${INCLUDES} -Xcompiler -fPIC -shared -o $@ $(c_objects) ${LDFLAGS} -lgpuKnnLibrary


## Test binary targets
#------------------------------
test:	$(UNITBIN)/unittest $(UNITBIN)/perftest

$(UNITBIN)/unittest:	unittest.cpp $(BIN)/libKraskov.so
	${GCC} ${CCFLAGS} -std=c++11 $< -Wl,-rpath,"$(abspath $(BIN))" -L$(BIN) -lKraskov -o $@ 

$(UNITBIN)/perftest: perftest.cpp $(BIN)/libKraskov.so
	${GCC} ${CCFLAGS} -DTIMER -std=c++11 $< -Wl,-rpath,"$(abspath $(BIN))" -L$(BIN) -lKraskov -o $@ 

clean:
	rm -f $(BIN)/* $(UNITBIN)/*

