module;

export module spatial_decomposition_device_kernels;

/**
 * \file kernel_sources.ixx
 * \brief The device kernels of the spatial-decomposition MD step, as source text shared by the backends.
 *
 * The three sources (pair_kernel_source.cpp, mesh_kernel_source.cpp, bonded_kernel_source.cpp) hold the physics
 * once. They are written in a small *kernel dialect*: plain C with the address spaces, work-item queries,
 * barriers, vector constructors and bit casts behind macros, which the backend that compiles them defines in a
 * header it prepends to the source (OpenCL C 1.2: opencl/dialect_source.cpp). A backend for another API
 * provides its own header; the kernels themselves do not change.
 *
 * The vocabulary:
 *
 *   KERNEL                         entry-point qualifier
 *   KERNEL_GROUP_SIZE(n)           entry point with a required work-group size of n work-items
 *   DEVICE_FUNCTION                qualifier of a non-entry device function
 *   GLOBAL / CONSTANT / LOCAL      address-space qualifiers of pointers and arrays
 *   RESTRICT                       the restrict qualifier
 *   GLOBAL_ID() / LOCAL_ID()       1D work-item indices (uint)
 *   GROUP_ID() / LOCAL_SIZE() / GLOBAL_SIZE()
 *   LOCAL_BARRIER()                work-group barrier with a local-memory fence
 *   FLOAT2(...) / FLOAT3(...) / FLOAT4(...) / INT2(...)   vector constructors (all component forms)
 *   AS_UINT(x) / AS_FLOAT(x)       bit casts between float and uint
 *   ATOMIC_ADD_INT(pointer, value) atomic add on a GLOBAL int
 *
 * Everything else (float2/3/4 and int2 types with .xyz swizzles, uint, size_t, the C math library, dot, length,
 * cross, clamp, rint, erf, exp, min/max, popcount-free bit arithmetic, `#pragma unroll`, `#define`) is the common
 * subset of OpenCL C and the Metal Shading Language; a CUDA header additionally has to supply the swizzles.
 *
 * The host-side mirrors of the parameter structs the kernels read (Parameters, BuildParameters, MeshParameters,
 * BondedParameters, Term, MoleculeInfo) live in backend.ixx and bonded_topology.ixx with static_asserts on their
 * sizes; the two sides must be edited together.
 */

/// The list build, the pruning and the pair kernel (pair_kernel_source.cpp).
export extern const char* const deviceKernelPairSource;
/// The particle-mesh Ewald sum (mesh_kernel_source.cpp); the B-spline order is the compile-time MESH_ORDER.
export extern const char* const deviceKernelMeshSource;
/// The per-molecule terms (bonded_kernel_source.cpp).
export extern const char* const deviceKernelBondedSource;
