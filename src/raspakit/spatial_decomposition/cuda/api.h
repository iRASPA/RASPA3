// The subset of the CUDA driver API (cuda.h) and of NVRTC (nvrtc.h) used by the CUDA backend of the
// spatial-decomposition device step, declared here so that the backend builds without the CUDA toolkit: the two
// libraries are loaded at run time (api.cpp), the executables carry no link-time dependency on them. The types,
// enumerators and prototypes mirror the ABI of CUDA 11.2 and later (the versioned entry points are named by
// their exported symbols, e.g. cuMemAlloc_v2).

#pragma once

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

  typedef int CUresult;
  typedef int CUdevice;
  typedef struct CUctx_st* CUcontext;
  typedef struct CUmod_st* CUmodule;
  typedef struct CUfunc_st* CUfunction;
  typedef struct CUstream_st* CUstream;
  typedef struct CUevent_st* CUevent;
  typedef unsigned long long CUdeviceptr;
  typedef int CUdevice_attribute;
  typedef int CUfunction_attribute;
  typedef int CUjit_option;

  enum
  {
    RASPA_CUDA_SUCCESS = 0,
    RASPA_CUDA_ERROR_NO_DEVICE = 100,
    RASPA_CUDA_ERROR_NOT_FOUND = 500
  };

  enum
  {
    RASPA_CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 1,
    RASPA_CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK = 8,
    RASPA_CU_DEVICE_ATTRIBUTE_CLOCK_RATE = 13,
    RASPA_CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT = 16,
    RASPA_CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR = 75,
    RASPA_CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR = 76,
    RASPA_CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK_OPTIN = 97,
    RASPA_CU_DEVICE_ATTRIBUTE_MEMORY_POOLS_SUPPORTED = 115
  };

  enum
  {
    RASPA_CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 0,
    RASPA_CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES = 1,
    RASPA_CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES = 8
  };

  enum
  {
    RASPA_CU_JIT_INFO_LOG_BUFFER = 3,
    RASPA_CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES = 4,
    RASPA_CU_JIT_ERROR_LOG_BUFFER = 5,
    RASPA_CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES = 6
  };

  enum
  {
    RASPA_CU_STREAM_NON_BLOCKING = 0x1,
    RASPA_CU_EVENT_DISABLE_TIMING = 0x2,
    RASPA_CU_MEMHOSTALLOC_PORTABLE = 0x1
  };

  typedef CUresult (*raspa_cuInit)(unsigned int flags);
  typedef CUresult (*raspa_cuDriverGetVersion)(int* version);
  typedef CUresult (*raspa_cuGetErrorString)(CUresult error, const char** string);
  typedef CUresult (*raspa_cuGetErrorName)(CUresult error, const char** string);
  typedef CUresult (*raspa_cuDeviceGetCount)(int* count);
  typedef CUresult (*raspa_cuDeviceGet)(CUdevice* device, int ordinal);
  typedef CUresult (*raspa_cuDeviceGetName)(char* name, int length, CUdevice device);
  typedef CUresult (*raspa_cuDeviceGetAttribute)(int* value, CUdevice_attribute attribute, CUdevice device);
  typedef CUresult (*raspa_cuDevicePrimaryCtxRetain)(CUcontext* context, CUdevice device);
  typedef CUresult (*raspa_cuDevicePrimaryCtxRelease)(CUdevice device);
  typedef CUresult (*raspa_cuCtxSetCurrent)(CUcontext context);
  typedef CUresult (*raspa_cuCtxGetCurrent)(CUcontext* context);
  typedef CUresult (*raspa_cuStreamCreate)(CUstream* stream, unsigned int flags);
  typedef CUresult (*raspa_cuStreamDestroy)(CUstream stream);
  typedef CUresult (*raspa_cuStreamSynchronize)(CUstream stream);
  typedef CUresult (*raspa_cuEventCreate)(CUevent* event, unsigned int flags);
  typedef CUresult (*raspa_cuEventDestroy)(CUevent event);
  typedef CUresult (*raspa_cuEventRecord)(CUevent event, CUstream stream);
  typedef CUresult (*raspa_cuEventSynchronize)(CUevent event);
  typedef CUresult (*raspa_cuMemAlloc)(CUdeviceptr* pointer, size_t bytes);
  typedef CUresult (*raspa_cuMemFree)(CUdeviceptr pointer);
  typedef CUresult (*raspa_cuMemAllocAsync)(CUdeviceptr* pointer, size_t bytes, CUstream stream);
  typedef CUresult (*raspa_cuMemFreeAsync)(CUdeviceptr pointer, CUstream stream);
  typedef CUresult (*raspa_cuMemHostAlloc)(void** pointer, size_t bytes, unsigned int flags);
  typedef CUresult (*raspa_cuMemFreeHost)(void* pointer);
  typedef CUresult (*raspa_cuMemcpyHtoD)(CUdeviceptr destination, const void* source, size_t bytes);
  typedef CUresult (*raspa_cuMemcpyDtoH)(void* destination, CUdeviceptr source, size_t bytes);
  typedef CUresult (*raspa_cuMemcpyHtoDAsync)(CUdeviceptr destination, const void* source, size_t bytes,
                                              CUstream stream);
  typedef CUresult (*raspa_cuMemcpyDtoHAsync)(void* destination, CUdeviceptr source, size_t bytes, CUstream stream);
  typedef CUresult (*raspa_cuMemcpyDtoDAsync)(CUdeviceptr destination, CUdeviceptr source, size_t bytes,
                                              CUstream stream);
  typedef CUresult (*raspa_cuModuleLoadDataEx)(CUmodule* module, const void* image, unsigned int numberOfOptions,
                                               CUjit_option* options, void** optionValues);
  typedef CUresult (*raspa_cuModuleUnload)(CUmodule module);
  typedef CUresult (*raspa_cuModuleGetFunction)(CUfunction* function, CUmodule module, const char* name);
  typedef CUresult (*raspa_cuFuncGetAttribute)(int* value, CUfunction_attribute attribute, CUfunction function);
  typedef CUresult (*raspa_cuFuncSetAttribute)(CUfunction function, CUfunction_attribute attribute, int value);
  typedef CUresult (*raspa_cuLaunchKernel)(CUfunction function, unsigned int gridX, unsigned int gridY,
                                           unsigned int gridZ, unsigned int blockX, unsigned int blockY,
                                           unsigned int blockZ, unsigned int sharedMemoryBytes, CUstream stream,
                                           void** kernelParameters, void** extra);

  typedef int nvrtcResult;
  typedef struct _nvrtcProgram* nvrtcProgram;

  enum
  {
    RASPA_NVRTC_SUCCESS = 0
  };

  typedef const char* (*raspa_nvrtcGetErrorString)(nvrtcResult result);
  typedef nvrtcResult (*raspa_nvrtcVersion)(int* major, int* minor);
  typedef nvrtcResult (*raspa_nvrtcGetNumSupportedArchs)(int* count);
  typedef nvrtcResult (*raspa_nvrtcGetSupportedArchs)(int* archs);
  typedef nvrtcResult (*raspa_nvrtcCreateProgram)(nvrtcProgram* program, const char* source, const char* name,
                                                  int numberOfHeaders, const char* const* headers,
                                                  const char* const* includeNames);
  typedef nvrtcResult (*raspa_nvrtcDestroyProgram)(nvrtcProgram* program);
  typedef nvrtcResult (*raspa_nvrtcCompileProgram)(nvrtcProgram program, int numberOfOptions,
                                                   const char* const* options);
  typedef nvrtcResult (*raspa_nvrtcGetProgramLogSize)(nvrtcProgram program, size_t* size);
  typedef nvrtcResult (*raspa_nvrtcGetProgramLog)(nvrtcProgram program, char* log);
  typedef nvrtcResult (*raspa_nvrtcGetPTXSize)(nvrtcProgram program, size_t* size);
  typedef nvrtcResult (*raspa_nvrtcGetPTX)(nvrtcProgram program, char* ptx);
  typedef nvrtcResult (*raspa_nvrtcGetCUBINSize)(nvrtcProgram program, size_t* size);
  typedef nvrtcResult (*raspa_nvrtcGetCUBIN)(nvrtcProgram program, char* cubin);

#ifdef __cplusplus
}
#endif
