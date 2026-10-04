module;

export module spatial_decomposition_device_kernels;

/**
 * \file kernel_sources.ixx
 * \brief The device kernels of the spatial-decomposition MD step, as source text shared by the backends.
 *
 * The four sources (pair_kernel_source.cpp, mesh_kernel_source.cpp, bonded_kernel_source.cpp,
 * resident_kernel_source.cpp) hold the physics once. They are written in a small *kernel dialect*: plain C with
 * the address spaces, work-item queries, barriers, vector constructors and bit casts behind macros, which the
 * backend that compiles them defines in a header it prepends to the source (OpenCL C 1.2:
 * opencl/dialect_source.cpp, Metal: metal/dialect_source.cpp, CUDA/NVRTC: cuda/dialect_source.cpp). A backend
 * for another API provides its own header; the kernels themselves do not change.
 *
 * The vocabulary:
 *
 *   KERNEL                         entry-point qualifier
 *   KERNEL_GROUP_SIZE(n)           entry point with a required work-group size of n work-items
 *   KERNEL_INDEX_ARGS              last token of every kernel parameter list (no comma before it): the languages
 *                                  without implicit work-item queries declare the index parameters here
 *   VALUE_ARG(type, name)          a scalar kernel argument passed by value (set with DeviceArg::value)
 *   DEVICE_FUNCTION                qualifier of a non-entry device function
 *   GLOBAL / CONSTANT / LOCAL / PRIVATE   address-space qualifiers of pointers (PRIVATE: pointers to work-item
 *                                  memory, i.e. function parameters pointing at local variables)
 *   LOCAL_DECL(declaration)        a work-group (local-memory) variable or array declared in a kernel body
 *   LOCAL_ARG(type, name)          a kernel parameter pointing at local memory of the size given by the launch
 *                                  (DeviceArg::local); declared after the buffer and value arguments
 *   LOCAL_ARG_BIND(type, name)     first statements of a kernel with LOCAL_ARG parameters, one per parameter in
 *                                  order: the languages that size local memory at the launch (CUDA) bind the
 *                                  pointer here; the others expand to nothing
 *   RESTRICT                       the restrict qualifier
 *   GLOBAL_ID() / LOCAL_ID()       1D work-item indices (uint); valid in kernel bodies only, device functions
 *   GROUP_ID() / LOCAL_SIZE() / GLOBAL_SIZE()   receive them as parameters
 *   LOCAL_BARRIER()                work-group barrier with a local-memory fence
 *   FLOAT2(...) / FLOAT3(...) / FLOAT4(...) / INT2(...)   vector constructors (all component forms)
 *   AS_UINT(x) / AS_FLOAT(x)       bit casts between float and uint
 *   ATOMIC_ADD_INT(pointer, value) atomic add on a GLOBAL int
 *   erf / erfc                     single-precision; supplied by the header where the language lacks them
 *
 * Everything else (float2/3/4 and int2 types with .xyz swizzles, uint, size_t, the C math library, dot, length,
 * cross, clamp, rint, exp, min/max, FLT_MAX, popcount-free bit arithmetic, `#pragma unroll`, `#define`) is the
 * common subset of OpenCL C and the Metal Shading Language; the CUDA header (cuda/dialect_source.cpp)
 * additionally supplies the vector types with their swizzles and operators, and the geometric functions.
 * Kernel buffer arguments are bound by their position among the buffer and value arguments, local-memory
 * arguments by their position among the local-memory arguments (DeviceContext::launch).
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
/// The resident velocity-Verlet integrator in double-float arithmetic (resident_kernel_source.cpp); compile with
/// DeviceMath::Strict.
export extern const char* const deviceKernelResidentSource;
