module;

module spatial_decomposition_cuda_context;

// The CUDA C++ side of the kernel dialect of spatial_decomposition_device_kernels (device/kernel_sources.ixx):
// prepended to the shared kernel sources by CUDAContext::compileKernel and compiled by NVRTC.
//
// CUDA has no address-space qualifiers on pointers (GLOBAL / CONSTANT / LOCAL / PRIVATE expand to nothing), sizes
// local memory at the launch (one `extern __shared__` block, carved by LOCAL_ARG_BIND from the byte offsets the
// context appends to the kernel arguments), and its built-in vector types have neither operators nor swizzles:
// the dialect substitutes its own float2 / float3 / float4 with the layout of the OpenCL and Metal types where it
// matters (float2: 8 bytes, float4: 16 bytes aligned; float3 never crosses the host boundary) and supplies the
// operators, constructors, geometric functions and vector-wide math the kernels use. Kernels are `extern "C"` so
// that cuModuleGetFunction finds them by their plain names.
const char* const cudaKernelDialect = R"CUDA(
typedef unsigned int uint;
#ifndef FLT_MAX
#define FLT_MAX 3.402823466e+38f
#endif

extern __shared__ __align__(16) unsigned char _raspa_shared[];

#define KERNEL extern "C" __global__
#define KERNEL_GROUP_SIZE(n) extern "C" __global__ __launch_bounds__(n)
#define KERNEL_INDEX_ARGS
#define VALUE_ARG(type, name) const type name
#define DEVICE_FUNCTION static __device__ __forceinline__
#define GLOBAL
#define CONSTANT
#define LOCAL
#define LOCAL_DECL(declaration) __shared__ declaration
#define LOCAL_ARG(type, name) const unsigned int _raspa_local_##name
#define LOCAL_ARG_BIND(type, name) type* name = reinterpret_cast<type*>(_raspa_shared + _raspa_local_##name)
#define PRIVATE
#define RESTRICT __restrict__
#define GLOBAL_ID() ((uint)(blockIdx.x * blockDim.x + threadIdx.x))
#define LOCAL_ID() ((uint)threadIdx.x)
#define GROUP_ID() ((uint)blockIdx.x)
#define LOCAL_SIZE() ((uint)blockDim.x)
#define GLOBAL_SIZE() ((uint)(gridDim.x * blockDim.x))
#define LOCAL_BARRIER() __syncthreads()
#define FLOAT2(...) raspa_make_float2(__VA_ARGS__)
#define FLOAT3(...) raspa_make_float3(__VA_ARGS__)
#define FLOAT4(...) raspa_make_float4(__VA_ARGS__)
#define INT2(...) make_int2(__VA_ARGS__)
#define AS_UINT(x) __float_as_uint(x)
#define AS_FLOAT(x) __uint_as_float(x)
#define ATOMIC_ADD_INT(pointer, value) atomicAdd((int*)(pointer), (int)(value))

// the vector types of the dialect (plain aggregates: usable in __shared__ arrays and in the parameter structs).
// NVRTC has no anonymous structs in unions, so the .xyz swizzle (read-only in the kernels) is a member function
// reached through the macro on its name.
struct __align__(8) raspa_float2
{
  float x, y;
};
struct raspa_float3
{
  float x, y, z;
};
struct __align__(16) raspa_float4
{
  float x, y, z, w;
  __device__ __forceinline__ raspa_float3 xyz() const
  {
    raspa_float3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
  }
};
#define xyz xyz()
#define float2 raspa_float2
#define float3 raspa_float3
#define float4 raspa_float4

__device__ __forceinline__ float2 raspa_make_float2(float x, float y)
{
  float2 v;
  v.x = x;
  v.y = y;
  return v;
}
__device__ __forceinline__ float2 raspa_make_float2(float s) { return raspa_make_float2(s, s); }
__device__ __forceinline__ float3 raspa_make_float3(float x, float y, float z)
{
  float3 v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}
__device__ __forceinline__ float3 raspa_make_float3(float s) { return raspa_make_float3(s, s, s); }
__device__ __forceinline__ float3 raspa_make_float3(float2 a, float z) { return raspa_make_float3(a.x, a.y, z); }
__device__ __forceinline__ float3 raspa_make_float3(float x, float2 a) { return raspa_make_float3(x, a.x, a.y); }
__device__ __forceinline__ float4 raspa_make_float4(float x, float y, float z, float w)
{
  float4 v;
  v.x = x;
  v.y = y;
  v.z = z;
  v.w = w;
  return v;
}
__device__ __forceinline__ float4 raspa_make_float4(float s) { return raspa_make_float4(s, s, s, s); }
__device__ __forceinline__ float4 raspa_make_float4(float3 a, float w) { return raspa_make_float4(a.x, a.y, a.z, w); }
__device__ __forceinline__ float4 raspa_make_float4(float x, float3 a) { return raspa_make_float4(x, a.x, a.y, a.z); }
__device__ __forceinline__ float4 raspa_make_float4(float2 a, float2 b) { return raspa_make_float4(a.x, a.y, b.x, b.y); }
__device__ __forceinline__ float4 raspa_make_float4(float2 a, float z, float w) { return raspa_make_float4(a.x, a.y, z, w); }
__device__ __forceinline__ float4 raspa_make_float4(float x, float y, float2 b) { return raspa_make_float4(x, y, b.x, b.y); }

// arithmetic: component-wise between vectors, broadcast with a scalar
#define RASPA_VEC2_OPS(OP)                                                                                   \
  __device__ __forceinline__ float2 operator OP(float2 a, float2 b) { return raspa_make_float2(a.x OP b.x, a.y OP b.y); } \
  __device__ __forceinline__ float2 operator OP(float2 a, float s) { return raspa_make_float2(a.x OP s, a.y OP s); }      \
  __device__ __forceinline__ float2 operator OP(float s, float2 b) { return raspa_make_float2(s OP b.x, s OP b.y); }      \
  __device__ __forceinline__ float2& operator OP##=(float2& a, float2 b) { a = a OP b; return a; }                  \
  __device__ __forceinline__ float2& operator OP##=(float2& a, float s) { a = a OP s; return a; }
#define RASPA_VEC3_OPS(OP)                                                                                   \
  __device__ __forceinline__ float3 operator OP(float3 a, float3 b) { return raspa_make_float3(a.x OP b.x, a.y OP b.y, a.z OP b.z); } \
  __device__ __forceinline__ float3 operator OP(float3 a, float s) { return raspa_make_float3(a.x OP s, a.y OP s, a.z OP s); }      \
  __device__ __forceinline__ float3 operator OP(float s, float3 b) { return raspa_make_float3(s OP b.x, s OP b.y, s OP b.z); }      \
  __device__ __forceinline__ float3& operator OP##=(float3& a, float3 b) { a = a OP b; return a; }                  \
  __device__ __forceinline__ float3& operator OP##=(float3& a, float s) { a = a OP s; return a; }
#define RASPA_VEC4_OPS(OP)                                                                                   \
  __device__ __forceinline__ float4 operator OP(float4 a, float4 b) { return raspa_make_float4(a.x OP b.x, a.y OP b.y, a.z OP b.z, a.w OP b.w); } \
  __device__ __forceinline__ float4 operator OP(float4 a, float s) { return raspa_make_float4(a.x OP s, a.y OP s, a.z OP s, a.w OP s); }      \
  __device__ __forceinline__ float4 operator OP(float s, float4 b) { return raspa_make_float4(s OP b.x, s OP b.y, s OP b.z, s OP b.w); }      \
  __device__ __forceinline__ float4& operator OP##=(float4& a, float4 b) { a = a OP b; return a; }                  \
  __device__ __forceinline__ float4& operator OP##=(float4& a, float s) { a = a OP s; return a; }
RASPA_VEC2_OPS(+)
RASPA_VEC2_OPS(-)
RASPA_VEC2_OPS(*)
RASPA_VEC2_OPS(/)
RASPA_VEC3_OPS(+)
RASPA_VEC3_OPS(-)
RASPA_VEC3_OPS(*)
RASPA_VEC3_OPS(/)
RASPA_VEC4_OPS(+)
RASPA_VEC4_OPS(-)
RASPA_VEC4_OPS(*)
RASPA_VEC4_OPS(/)
#undef RASPA_VEC2_OPS
#undef RASPA_VEC3_OPS
#undef RASPA_VEC4_OPS
__device__ __forceinline__ float2 operator-(float2 a) { return raspa_make_float2(-a.x, -a.y); }
__device__ __forceinline__ float3 operator-(float3 a) { return raspa_make_float3(-a.x, -a.y, -a.z); }
__device__ __forceinline__ float4 operator-(float4 a) { return raspa_make_float4(-a.x, -a.y, -a.z, -a.w); }

// geometric functions
__device__ __forceinline__ float dot(float2 a, float2 b) { return a.x * b.x + a.y * b.y; }
__device__ __forceinline__ float dot(float3 a, float3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
__device__ __forceinline__ float dot(float4 a, float4 b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }
__device__ __forceinline__ float length(float2 a) { return sqrtf(dot(a, a)); }
__device__ __forceinline__ float length(float3 a) { return sqrtf(dot(a, a)); }
__device__ __forceinline__ float length(float4 a) { return sqrtf(dot(a, a)); }
__device__ __forceinline__ float3 cross(float3 a, float3 b)
{
  return raspa_make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
__device__ __forceinline__ float3 normalize(float3 a) { return a * rsqrtf(dot(a, a)); }

// common functions and vector-wide math
__device__ __forceinline__ float clamp(float x, float lo, float hi) { return fminf(fmaxf(x, lo), hi); }
__device__ __forceinline__ int clamp(int x, int lo, int hi) { return min(max(x, lo), hi); }
__device__ __forceinline__ uint clamp(uint x, uint lo, uint hi) { return min(max(x, lo), hi); }
__device__ __forceinline__ float2 fmin(float2 a, float2 b) { return raspa_make_float2(fminf(a.x, b.x), fminf(a.y, b.y)); }
__device__ __forceinline__ float2 fmax(float2 a, float2 b) { return raspa_make_float2(fmaxf(a.x, b.x), fmaxf(a.y, b.y)); }
__device__ __forceinline__ float3 fmin(float3 a, float3 b)
{
  return raspa_make_float3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z));
}
__device__ __forceinline__ float3 fmax(float3 a, float3 b)
{
  return raspa_make_float3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z));
}
__device__ __forceinline__ float4 fmin(float4 a, float4 b)
{
  return raspa_make_float4(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z), fminf(a.w, b.w));
}
__device__ __forceinline__ float4 fmax(float4 a, float4 b)
{
  return raspa_make_float4(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z), fmaxf(a.w, b.w));
}
__device__ __forceinline__ float2 rint(float2 a) { return raspa_make_float2(rintf(a.x), rintf(a.y)); }
__device__ __forceinline__ float3 rint(float3 a) { return raspa_make_float3(rintf(a.x), rintf(a.y), rintf(a.z)); }
__device__ __forceinline__ float2 floor(float2 a) { return raspa_make_float2(floorf(a.x), floorf(a.y)); }
__device__ __forceinline__ float3 floor(float3 a) { return raspa_make_float3(floorf(a.x), floorf(a.y), floorf(a.z)); }
__device__ __forceinline__ float2 fabs(float2 a) { return raspa_make_float2(fabsf(a.x), fabsf(a.y)); }
__device__ __forceinline__ float3 fabs(float3 a) { return raspa_make_float3(fabsf(a.x), fabsf(a.y), fabsf(a.z)); }
)CUDA";
