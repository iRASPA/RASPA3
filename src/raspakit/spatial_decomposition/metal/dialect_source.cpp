module;

module spatial_decomposition_metal_context;

// The Metal Shading Language side of the kernel dialect of spatial_decomposition_device_kernels
// (device/kernel_sources.ixx): prepended to the shared kernel sources by MetalContext::compileKernel.
//
// MSL has no implicit work-item queries: every kernel declares them through KERNEL_INDEX_ARGS. Buffer and
// by-value arguments take the buffer indices 0.. in declaration order (MetalContext::launch binds them by the
// same count), threadgroup pointer arguments the threadgroup indices 0... MSL has no erf/erfc either: the ones
// below are single precision to ~1e-7 (Taylor series below 1, the Chebyshev fit of Numerical Recipes above).
const char* const metalKernelDialect = R"MSL(
#include <metal_stdlib>
using namespace metal;

#define KERNEL kernel
#define KERNEL_GROUP_SIZE(n) [[max_total_threads_per_threadgroup(n)]] kernel
#define KERNEL_INDEX_ARGS , uint _raspa_gid [[thread_position_in_grid]], uint _raspa_lid [[thread_index_in_threadgroup]], uint _raspa_grp [[threadgroup_position_in_grid]], uint _raspa_lsz [[threads_per_threadgroup]], uint _raspa_gsz [[threads_per_grid]]
#define VALUE_ARG(type, name) constant type& name
#define DEVICE_FUNCTION static inline
#define GLOBAL device
#define CONSTANT constant
#define LOCAL threadgroup
#define LOCAL_DECL(declaration) threadgroup declaration
#define LOCAL_ARG(type, name) threadgroup type* name
#define LOCAL_ARG_BIND(type, name)
#define PRIVATE thread
#define RESTRICT
#define GLOBAL_ID() (_raspa_gid)
#define LOCAL_ID() (_raspa_lid)
#define GROUP_ID() (_raspa_grp)
#define LOCAL_SIZE() (_raspa_lsz)
#define GLOBAL_SIZE() (_raspa_gsz)
#define LOCAL_BARRIER() threadgroup_barrier(mem_flags::mem_threadgroup)
#define FLOAT2(...) float2(__VA_ARGS__)
#define FLOAT3(...) float3(__VA_ARGS__)
#define FLOAT4(...) float4(__VA_ARGS__)
#define INT2(...) int2(__VA_ARGS__)
#define AS_UINT(x) as_type<uint>(x)
#define AS_FLOAT(x) as_type<float>(x)
#define ATOMIC_ADD_INT(pointer, value) atomic_fetch_add_explicit(reinterpret_cast<device atomic_int*>(pointer), (value), memory_order_relaxed)

// erfc(x) for x >= 0: Chebyshev fit (Numerical Recipes erfcc), fractional error below 1.2e-7 everywhere
static inline float raspa_erfc_positive(float x)
{
  const float t = 1.0f / (1.0f + 0.5f * x);
  const float poly =
      -1.26551223f +
      t * (1.00002368f +
           t * (0.37409196f +
                t * (0.09678418f +
                     t * (-0.18628806f +
                          t * (0.27886807f + t * (-1.13520398f + t * (1.48851587f + t * (-0.82215223f + t * 0.17087277f))))))));
  return t * exp(-x * x + poly);
}

// erf(x) for |x| <= 1: Taylor series, 11 terms (the omitted terms are below 1e-8 at x = 1)
static inline float raspa_erf_small(float x)
{
  const float x2 = x * x;
  float term = x;
  float sum = x;
  for (int n = 1; n <= 10; ++n)
  {
    term *= -x2 / (float)n;
    sum += term / (float)(2 * n + 1);
  }
  return 1.1283791670955126f * sum;
}

static inline float erf(float x)
{
  const float a = fabs(x);
  const float value = a <= 1.0f ? raspa_erf_small(a) : 1.0f - raspa_erfc_positive(a);
  return copysign(value, x);
}

static inline float erfc(float x)
{
  const float a = fabs(x);
  const float value = a <= 1.0f ? 1.0f - raspa_erf_small(a) : raspa_erfc_positive(a);
  return x >= 0.0f ? value : 2.0f - value;
}
)MSL";
