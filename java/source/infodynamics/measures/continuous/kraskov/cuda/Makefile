
## Set these if the CUDA Toolkit is not in your PATH
# Location of the CUDA Toolkit binaries and libraries
# CUDA_PATH       ?= /usr/local/cuda
# CUDA_INC_PATH   ?= $(CUDA_PATH)/include
# CUDA_LIB_PATH   ?= $(CUDA_PATH)/lib
# CUDA_BIN_PATH   ?= $(CUDA_PATH)/bin

# Common binaries and flags. TODO: the architecture flag may not be necessary
NVCC       ?= nvcc
GCC        ?= g++
CCFLAGS    ?= -O3
# NVCCFLAGS  ?= -O3 -arch=sm_30
NVCCFLAGS  ?= -O3

ifdef DEBUG
	CCFLAGS    += -g
	NVCCFLAGS  += -g -G
endif

# Common includes and paths for CUDA. This assumes the CUDA toolkit is in PATH
INCLUDES      := -I. -I/usr/lib/jvm/java-7-openjdk-amd64/include
LDFLAGS       += -L. -lcuda -lcudart

# Target rules
all:  libKraskov.so

gpuKnnLibrary.o: gpuKnnLibrary.cu helperfunctions.cu gpuKnnBF_kernel.cu
	${NVCC} ${NVCCFLAGS} ${INCLUDES} -Xcompiler -fPIC -c gpuKnnLibrary.cu

libgpuKnnLibrary.a: gpuKnnLibrary.o
	${AR} -r libgpuKnnLibrary.a gpuKnnLibrary.o

libKraskov.so: libgpuKnnLibrary.a kraskovCuda.c
	${NVCC} ${NVCCFLAGS} ${INCLUDES} -Xcompiler -fPIC -shared -o libKraskov.so kraskovCuda.c ${LDFLAGS} -lgpuKnnLibrary

clean:
	rm -f *.o *.so *.a

