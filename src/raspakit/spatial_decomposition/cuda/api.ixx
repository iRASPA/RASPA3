module;

#include "api.h"

export module spatial_decomposition_cuda_api;

import std;

/**
 * \brief The CUDA driver API and NVRTC, loaded at run time.
 *
 * The CUDA backend has no link-time dependency on CUDA: libcuda (the driver, present wherever an NVIDIA driver is
 * installed) and libnvrtc (the run-time compiler, from the CUDA toolkit or the `nvidia-cuda-nvrtc` Python wheel)
 * are opened with dlopen on first use and their entry points resolved into this table. Where either is missing
 * the backend reports no device, like the OpenCL ICD loader does without a platform.
 */
export struct CUDAApi
{
  raspa_cuInit cuInit{nullptr};
  raspa_cuDriverGetVersion cuDriverGetVersion{nullptr};
  raspa_cuGetErrorString cuGetErrorString{nullptr};
  raspa_cuGetErrorName cuGetErrorName{nullptr};
  raspa_cuDeviceGetCount cuDeviceGetCount{nullptr};
  raspa_cuDeviceGet cuDeviceGet{nullptr};
  raspa_cuDeviceGetName cuDeviceGetName{nullptr};
  raspa_cuDeviceGetAttribute cuDeviceGetAttribute{nullptr};
  raspa_cuDevicePrimaryCtxRetain cuDevicePrimaryCtxRetain{nullptr};
  raspa_cuDevicePrimaryCtxRelease cuDevicePrimaryCtxRelease{nullptr};
  raspa_cuCtxSetCurrent cuCtxSetCurrent{nullptr};
  raspa_cuCtxGetCurrent cuCtxGetCurrent{nullptr};
  raspa_cuStreamCreate cuStreamCreate{nullptr};
  raspa_cuStreamDestroy cuStreamDestroy{nullptr};
  raspa_cuStreamSynchronize cuStreamSynchronize{nullptr};
  raspa_cuEventCreate cuEventCreate{nullptr};
  raspa_cuEventDestroy cuEventDestroy{nullptr};
  raspa_cuEventRecord cuEventRecord{nullptr};
  raspa_cuEventSynchronize cuEventSynchronize{nullptr};
  raspa_cuMemAlloc cuMemAlloc{nullptr};
  raspa_cuMemFree cuMemFree{nullptr};
  raspa_cuMemAllocAsync cuMemAllocAsync{nullptr};  ///< may be null (driver below 11.2)
  raspa_cuMemFreeAsync cuMemFreeAsync{nullptr};    ///< may be null (driver below 11.2)
  raspa_cuMemHostAlloc cuMemHostAlloc{nullptr};
  raspa_cuMemFreeHost cuMemFreeHost{nullptr};
  raspa_cuMemcpyHtoD cuMemcpyHtoD{nullptr};
  raspa_cuMemcpyDtoH cuMemcpyDtoH{nullptr};
  raspa_cuMemcpyHtoDAsync cuMemcpyHtoDAsync{nullptr};
  raspa_cuMemcpyDtoHAsync cuMemcpyDtoHAsync{nullptr};
  raspa_cuMemcpyDtoDAsync cuMemcpyDtoDAsync{nullptr};
  raspa_cuModuleLoadDataEx cuModuleLoadDataEx{nullptr};
  raspa_cuModuleUnload cuModuleUnload{nullptr};
  raspa_cuModuleGetFunction cuModuleGetFunction{nullptr};
  raspa_cuFuncGetAttribute cuFuncGetAttribute{nullptr};
  raspa_cuFuncSetAttribute cuFuncSetAttribute{nullptr};
  raspa_cuLaunchKernel cuLaunchKernel{nullptr};

  raspa_nvrtcGetErrorString nvrtcGetErrorString{nullptr};
  raspa_nvrtcVersion nvrtcVersion{nullptr};
  raspa_nvrtcGetNumSupportedArchs nvrtcGetNumSupportedArchs{nullptr};  ///< may be null (NVRTC below 11.2)
  raspa_nvrtcGetSupportedArchs nvrtcGetSupportedArchs{nullptr};        ///< may be null (NVRTC below 11.2)
  raspa_nvrtcCreateProgram nvrtcCreateProgram{nullptr};
  raspa_nvrtcDestroyProgram nvrtcDestroyProgram{nullptr};
  raspa_nvrtcCompileProgram nvrtcCompileProgram{nullptr};
  raspa_nvrtcGetProgramLogSize nvrtcGetProgramLogSize{nullptr};
  raspa_nvrtcGetProgramLog nvrtcGetProgramLog{nullptr};
  raspa_nvrtcGetPTXSize nvrtcGetPTXSize{nullptr};
  raspa_nvrtcGetPTX nvrtcGetPTX{nullptr};
  raspa_nvrtcGetCUBINSize nvrtcGetCUBINSize{nullptr};
  raspa_nvrtcGetCUBIN nvrtcGetCUBIN{nullptr};

  /// The error text of a driver call.
  std::string errorString(CUresult result) const;
};

/// The loaded API, or nullptr when the driver or NVRTC could not be loaded (see cudaApiError); loads on first
/// call, thread-safe.
export const CUDAApi* cudaApi();
/// Why the API is unavailable (empty when it is).
export std::string cudaApiError();
