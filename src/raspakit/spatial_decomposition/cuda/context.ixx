module;

export module spatial_decomposition_cuda_context;

import std;

import spatial_decomposition_device_context;

/// Whether a CUDA device is usable on this machine (loads the driver and NVRTC on first use).
export bool cudaAvailable();
/// Name of the selected CUDA device (empty when none).
export std::string cudaDeviceName();
/// Why no CUDA device is usable (empty when one is): the missing library or the driver error.
export std::string cudaUnavailableReason();
/// The DeviceContext on the selected CUDA device (own stream on the primary context); throws when none.
export std::unique_ptr<DeviceContext> createCUDAContext();

/// The CUDA C++ definitions of the kernel dialect (dialect_source.cpp), prepended to the shared kernel sources.
extern const char* const cudaKernelDialect;
