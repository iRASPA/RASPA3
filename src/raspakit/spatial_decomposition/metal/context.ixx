module;

export module spatial_decomposition_metal_context;

import std;

import spatial_decomposition_device_context;

/// Whether a Metal device is available on this machine.
export bool metalAvailable();
/// Name of the default Metal device (empty when none).
export std::string metalDeviceName();
/// The DeviceContext on the default Metal device (own command queue); throws when none.
export std::unique_ptr<DeviceContext> createMetalContext();

/// The MSL definitions of the kernel dialect (dialect_source.cpp), prepended to the shared kernel sources.
extern const char* const metalKernelDialect;
