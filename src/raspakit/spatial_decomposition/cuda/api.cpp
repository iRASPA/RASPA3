module;

#include "api.h"

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <dlfcn.h>
#endif

module spatial_decomposition_cuda_api;

import std;

namespace
{
#if defined(_WIN32)
using LibraryHandle = HMODULE;
LibraryHandle openLibrary(const char* name) { return LoadLibraryA(name); }
void* lookup(LibraryHandle library, const char* name)
{
  return reinterpret_cast<void*>(GetProcAddress(library, name));
}
constexpr const char* driverNames[] = {"nvcuda.dll"};
constexpr const char* nvrtcNames[] = {"nvrtc64_130_0.dll", "nvrtc64_120_0.dll"};
#else
using LibraryHandle = void*;
LibraryHandle openLibrary(const char* name) { return dlopen(name, RTLD_NOW | RTLD_LOCAL); }
void* lookup(LibraryHandle library, const char* name) { return dlsym(library, name); }
#if defined(__APPLE__)
constexpr const char* driverNames[] = {"libcuda.dylib"};
constexpr const char* nvrtcNames[] = {"libnvrtc.dylib"};
#else
constexpr const char* driverNames[] = {"libcuda.so.1", "libcuda.so"};
constexpr const char* nvrtcNames[] = {"libnvrtc.so.13", "libnvrtc.so.12", "libnvrtc.so"};
#endif
#endif

/// Directories where the `nvidia-cuda-nvrtc` Python wheel installs libnvrtc: nvidia/cuda_nvrtc/lib below the
/// site-packages that hold this library (found from the address of this function).
std::vector<std::string> wheelDirectories()
{
  std::vector<std::string> directories;
#if !defined(_WIN32)
  Dl_info info{};
  if (dladdr(reinterpret_cast<void*>(&wheelDirectories), &info) != 0 && info.dli_fname != nullptr)
  {
    std::filesystem::path path = std::filesystem::path(info.dli_fname).parent_path();
    for (int level = 0; level < 4 && !path.empty(); ++level)
    {
      const std::filesystem::path candidate = path / "nvidia" / "cuda_nvrtc" / "lib";
      std::error_code error;
      if (std::filesystem::is_directory(candidate, error)) directories.push_back(candidate.string());
      path = path.parent_path();
    }
  }
#endif
  return directories;
}

LibraryHandle openFirst(std::span<const char* const> names, const std::vector<std::string>& directories)
{
  for (const char* name : names)
  {
    if (LibraryHandle handle = openLibrary(name)) return handle;
  }
  for (const std::string& directory : directories)
  {
    for (const char* name : names)
    {
      const std::string path = (std::filesystem::path(directory) / name).string();
      if (LibraryHandle handle = openLibrary(path.c_str())) return handle;
    }
  }
  return nullptr;
}

template <typename F>
bool resolve(LibraryHandle library, const char* name, F& function, std::string& missing)
{
  function = reinterpret_cast<F>(lookup(library, name));
  if (function == nullptr)
  {
    if (!missing.empty()) missing += ", ";
    missing += name;
    return false;
  }
  return true;
}

template <typename F>
void resolveOptional(LibraryHandle library, const char* name, F& function)
{
  function = reinterpret_cast<F>(lookup(library, name));
}

struct Loader
{
  CUDAApi api{};
  bool loaded{false};
  std::string error{};

  Loader()
  {
    const std::vector<std::string> directories = wheelDirectories();

    LibraryHandle driver = openFirst(driverNames, directories);
    if (driver == nullptr)
    {
      error = "the CUDA driver library (libcuda) is not installed";
      return;
    }

    LibraryHandle nvrtc = nullptr;
    if (const char* path = std::getenv("RASPA_NVRTC_LIBRARY"); path != nullptr && *path != '\0')
    {
      nvrtc = openLibrary(path);
    }
    if (nvrtc == nullptr) nvrtc = openFirst(nvrtcNames, directories);
    if (nvrtc == nullptr)
    {
      error =
          "the CUDA run-time compiler (libnvrtc) is not installed: install the CUDA toolkit or the "
          "nvidia-cuda-nvrtc Python package, or point RASPA_NVRTC_LIBRARY at the library";
      return;
    }

    std::string missing;
    resolve(driver, "cuInit", api.cuInit, missing);
    resolve(driver, "cuDriverGetVersion", api.cuDriverGetVersion, missing);
    resolve(driver, "cuGetErrorString", api.cuGetErrorString, missing);
    resolve(driver, "cuGetErrorName", api.cuGetErrorName, missing);
    resolve(driver, "cuDeviceGetCount", api.cuDeviceGetCount, missing);
    resolve(driver, "cuDeviceGet", api.cuDeviceGet, missing);
    resolve(driver, "cuDeviceGetName", api.cuDeviceGetName, missing);
    resolve(driver, "cuDeviceGetAttribute", api.cuDeviceGetAttribute, missing);
    resolve(driver, "cuDevicePrimaryCtxRetain", api.cuDevicePrimaryCtxRetain, missing);
    resolve(driver, "cuDevicePrimaryCtxRelease_v2", api.cuDevicePrimaryCtxRelease, missing);
    resolve(driver, "cuCtxSetCurrent", api.cuCtxSetCurrent, missing);
    resolve(driver, "cuCtxGetCurrent", api.cuCtxGetCurrent, missing);
    resolve(driver, "cuStreamCreate", api.cuStreamCreate, missing);
    resolve(driver, "cuStreamDestroy_v2", api.cuStreamDestroy, missing);
    resolve(driver, "cuStreamSynchronize", api.cuStreamSynchronize, missing);
    resolve(driver, "cuEventCreate", api.cuEventCreate, missing);
    resolve(driver, "cuEventDestroy_v2", api.cuEventDestroy, missing);
    resolve(driver, "cuEventRecord", api.cuEventRecord, missing);
    resolve(driver, "cuEventSynchronize", api.cuEventSynchronize, missing);
    resolve(driver, "cuMemAlloc_v2", api.cuMemAlloc, missing);
    resolve(driver, "cuMemFree_v2", api.cuMemFree, missing);
    resolveOptional(driver, "cuMemAllocAsync", api.cuMemAllocAsync);
    resolveOptional(driver, "cuMemFreeAsync", api.cuMemFreeAsync);
    resolve(driver, "cuMemHostAlloc", api.cuMemHostAlloc, missing);
    resolve(driver, "cuMemFreeHost", api.cuMemFreeHost, missing);
    resolve(driver, "cuMemcpyHtoD_v2", api.cuMemcpyHtoD, missing);
    resolve(driver, "cuMemcpyDtoH_v2", api.cuMemcpyDtoH, missing);
    resolve(driver, "cuMemcpyHtoDAsync_v2", api.cuMemcpyHtoDAsync, missing);
    resolve(driver, "cuMemcpyDtoHAsync_v2", api.cuMemcpyDtoHAsync, missing);
    resolve(driver, "cuMemcpyDtoDAsync_v2", api.cuMemcpyDtoDAsync, missing);
    resolve(driver, "cuModuleLoadDataEx", api.cuModuleLoadDataEx, missing);
    resolve(driver, "cuModuleUnload", api.cuModuleUnload, missing);
    resolve(driver, "cuModuleGetFunction", api.cuModuleGetFunction, missing);
    resolve(driver, "cuFuncGetAttribute", api.cuFuncGetAttribute, missing);
    resolve(driver, "cuFuncSetAttribute", api.cuFuncSetAttribute, missing);
    resolve(driver, "cuLaunchKernel", api.cuLaunchKernel, missing);

    resolve(nvrtc, "nvrtcGetErrorString", api.nvrtcGetErrorString, missing);
    resolve(nvrtc, "nvrtcVersion", api.nvrtcVersion, missing);
    resolveOptional(nvrtc, "nvrtcGetNumSupportedArchs", api.nvrtcGetNumSupportedArchs);
    resolveOptional(nvrtc, "nvrtcGetSupportedArchs", api.nvrtcGetSupportedArchs);
    resolve(nvrtc, "nvrtcCreateProgram", api.nvrtcCreateProgram, missing);
    resolve(nvrtc, "nvrtcDestroyProgram", api.nvrtcDestroyProgram, missing);
    resolve(nvrtc, "nvrtcCompileProgram", api.nvrtcCompileProgram, missing);
    resolve(nvrtc, "nvrtcGetProgramLogSize", api.nvrtcGetProgramLogSize, missing);
    resolve(nvrtc, "nvrtcGetProgramLog", api.nvrtcGetProgramLog, missing);
    resolve(nvrtc, "nvrtcGetPTXSize", api.nvrtcGetPTXSize, missing);
    resolve(nvrtc, "nvrtcGetPTX", api.nvrtcGetPTX, missing);
    resolve(nvrtc, "nvrtcGetCUBINSize", api.nvrtcGetCUBINSize, missing);
    resolve(nvrtc, "nvrtcGetCUBIN", api.nvrtcGetCUBIN, missing);

    if (!missing.empty())
    {
      error = std::format("the CUDA libraries lack the entry points {}", missing);
      return;
    }

    const CUresult result = api.cuInit(0);
    if (result != RASPA_CUDA_SUCCESS)
    {
      error = std::format("cuInit failed ({})", api.errorString(result));
      return;
    }
    loaded = true;
  }
};

const Loader& loader()
{
  static const Loader instance{};
  return instance;
}
}  // namespace

std::string CUDAApi::errorString(CUresult result) const
{
  const char* name = nullptr;
  const char* text = nullptr;
  if (cuGetErrorName != nullptr) cuGetErrorName(result, &name);
  if (cuGetErrorString != nullptr) cuGetErrorString(result, &text);
  return std::format("{}: {}", name != nullptr ? name : "CUDA error", text != nullptr ? text : "unknown");
}

const CUDAApi* cudaApi()
{
  const Loader& instance = loader();
  return instance.loaded ? &instance.api : nullptr;
}

std::string cudaApiError() { return loader().error; }
