module;

export module cbmc_device_kernels;

/**
 * \file kernel_sources.ixx
 * \brief The device kernel of the multiple-try CBMC grower, as source text shared by the backends.
 *
 * Written in the kernel dialect of spatial_decomposition_device_kernels (kernel_sources.ixx there documents the
 * vocabulary), so that the one source compiles with the OpenCL, Metal and CUDA dialect headers. The host-side
 * mirror of its parameter struct (CBMC::DeviceGrowParameters) lives in environment.ixx with a static_assert on
 * its size; the two must be edited together.
 */

/// The external energies of a batch of rigid trial molecules against a static background (grow_kernel_source.cpp).
export extern const char* const cbmcDeviceGrowKernelSource;
