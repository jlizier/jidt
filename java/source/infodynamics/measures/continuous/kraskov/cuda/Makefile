
## Set these if the CUDA Toolkit is not in your PATH
# Location of the CUDA Toolkit binaries and libraries
# CUDA_PATH       ?= /usr/local/cuda
# CUDA_INC_PATH   ?= $(CUDA_PATH)/include
# CUDA_LIB_PATH   ?= $(CUDA_PATH)/lib
# CUDA_BIN_PATH   ?= $(CUDA_PATH)/bin

# Common binaries and flags. TODO: the architecture flag may not be necessary
NVCC       ?= nvcc
GCC        ?= g++
CCFLAGS    ?= -O3 -Wall -Wextra -Wno-unused-parameter
NVCCFLAGS  ?= -O3 -arch=sm_30 -D_FORCE_INLINES -Xcompiler -Wall
# NVCCFLAGS  ?= -O3

ifdef DEBUG
	CCFLAGS    += -g
	NVCCFLAGS  += -g -G
endif

ifdef CUSTOMFLAGS
	CCFLAGS    += ${CUSTOMFLAGS}
	NVCCFLAGS  += ${CUSTOMFLAGS}
endif


# Common includes and paths for CUDA. This assumes the CUDA toolkit is in PATH
INCLUDES      := -I. -I./cub -I/usr/lib/jvm/java-8-openjdk-amd64/include -I/usr/lib/jvm/java-8-openjdk-amd64/include/linux
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

