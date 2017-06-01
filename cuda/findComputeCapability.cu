#include <cstdio>
int main() {
  int count = 0;
  if (cudaSuccess != cudaGetDeviceCount(&count)) return -1;
  if (count == 0) return -1;
  if (count == 1) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::printf("-arch=sm_%d%d", prop.major, prop.minor);
  } else {
    for (int device = 0; device < count; ++device) {
      cudaDeviceProp prop;
      if (cudaSuccess == cudaGetDeviceProperties(&prop, device)) {
        std::printf("-gencode arch=compute_%d%d,code=sm_%d%d ", prop.major, prop.minor, prop.major, prop.minor);
      }
    }
  }
  return 0;
}
