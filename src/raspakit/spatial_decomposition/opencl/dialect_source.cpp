module;

module spatial_decomposition_opencl_handles;

// The OpenCL C (1.2) side of the kernel dialect of spatial_decomposition_device_kernels (device/kernel_sources.ixx):
// prepended to the shared kernel sources by OpenCLDevice::buildDeviceProgram.
const char* const OpenCLDevice::openclKernelDialect = R"CLC(
#define KERNEL __kernel
#define KERNEL_GROUP_SIZE(n) __kernel __attribute__((reqd_work_group_size(n, 1, 1)))
#define KERNEL_INDEX_ARGS
#define VALUE_ARG(type, name) const type name
#define DEVICE_FUNCTION inline
#define GLOBAL __global
#define CONSTANT __constant
#define LOCAL __local
#define PRIVATE __private
#define RESTRICT restrict
#define GLOBAL_ID() ((uint)get_global_id(0))
#define LOCAL_ID() ((uint)get_local_id(0))
#define GROUP_ID() ((uint)get_group_id(0))
#define LOCAL_SIZE() ((uint)get_local_size(0))
#define GLOBAL_SIZE() ((uint)get_global_size(0))
#define LOCAL_BARRIER() barrier(CLK_LOCAL_MEM_FENCE)
#define FLOAT2(...) ((float2)(__VA_ARGS__))
#define FLOAT3(...) ((float3)(__VA_ARGS__))
#define FLOAT4(...) ((float4)(__VA_ARGS__))
#define INT2(...) ((int2)(__VA_ARGS__))
#define AS_UINT(x) as_uint(x)
#define AS_FLOAT(x) as_float(x)
#define ATOMIC_ADD_INT(pointer, value) atomic_add((pointer), (value))
)CLC";
