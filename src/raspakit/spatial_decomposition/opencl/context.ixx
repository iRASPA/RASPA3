module;

export module spatial_decomposition_opencl_context;

import std;

import spatial_decomposition_device_context;

/// Whether an OpenCL device is available (initializes the OpenCL runtime on first use).
export bool openclAvailable();
/// Name of the OpenCL device (empty when none).
export std::string openclDeviceName();
/// The DeviceContext on the OpenCL device of the `opencl` module (own in-order command queue); throws when none.
export std::unique_ptr<DeviceContext> createOpenCLContext();
