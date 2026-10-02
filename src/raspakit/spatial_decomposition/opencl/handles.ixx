module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module spatial_decomposition_opencl_handles;

import std;

/// RAII ownership of the OpenCL objects used by the spatial-decomposition device code, and the small helpers the
/// device classes share. Movable, not copyable: the device classes are values.
export namespace OpenCLDevice
{
template <typename T, void (*Release)(T)>
class Handle
{
 public:
  Handle() = default;
  explicit Handle(T h) : handle(h) {}
  ~Handle() { reset(); }
  Handle(const Handle&) = delete;
  Handle& operator=(const Handle&) = delete;
  Handle(Handle&& other) noexcept : handle(std::exchange(other.handle, nullptr)) {}
  Handle& operator=(Handle&& other) noexcept
  {
    if (this != &other)
    {
      reset();
      handle = std::exchange(other.handle, nullptr);
    }
    return *this;
  }
  void reset(T h = nullptr)
  {
    if (handle) Release(handle);
    handle = h;
  }
  T get() const { return handle; }

 private:
  T handle{nullptr};
};

inline void releaseMem(cl_mem m) { clReleaseMemObject(m); }
inline void releaseKernel(cl_kernel k) { clReleaseKernel(k); }
inline void releaseProgram(cl_program p) { clReleaseProgram(p); }
inline void releaseQueue(cl_command_queue q) { clReleaseCommandQueue(q); }
inline void releaseEvent(cl_event e) { clReleaseEvent(e); }

using MemHandle = Handle<cl_mem, &releaseMem>;
using KernelHandle = Handle<cl_kernel, &releaseKernel>;
using ProgramHandle = Handle<cl_program, &releaseProgram>;
using QueueHandle = Handle<cl_command_queue, &releaseQueue>;
using EventHandle = Handle<cl_event, &releaseEvent>;

/// Throws a runtime_error naming the failed call when `error` is not CL_SUCCESS.
inline void check(cl_int error, std::string_view what)
{
  if (error != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("[OpenCL]: {} failed (OpenCL error {})\n", what, error));
  }
}

inline std::size_t roundUp(std::size_t value, std::size_t multiple)
{
  return ((value + multiple - 1) / multiple) * multiple;
}

/// Compiles `source` for `device` with `options`; throws with the build log on failure.
inline ProgramHandle buildProgram(cl_context context, cl_device_id device, const char* source, const char* options,
                                  std::string_view who)
{
  cl_int error = CL_SUCCESS;
  ProgramHandle program(clCreateProgramWithSource(context, 1, &source, nullptr, &error));
  check(error, "clCreateProgramWithSource");
  error = clBuildProgram(program.get(), 1, &device, options, nullptr, nullptr);
  if (error != CL_SUCCESS)
  {
    std::size_t length = 0;
    clGetProgramBuildInfo(program.get(), device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &length);
    std::string log(length, '\0');
    clGetProgramBuildInfo(program.get(), device, CL_PROGRAM_BUILD_LOG, length, log.data(), nullptr);
    throw std::runtime_error(std::format("[{}]: the kernels failed to build:\n{}\n", who, log));
  }
  return program;
}

inline KernelHandle createKernel(cl_program program, const char* name)
{
  cl_int error = CL_SUCCESS;
  KernelHandle kernel(clCreateKernel(program, name, &error));
  check(error, std::format("clCreateKernel ({})", name));
  return kernel;
}

/// Sets consecutive buffer arguments 0..n-1 of a kernel.
inline void setBufferArguments(cl_kernel kernel, std::span<const cl_mem> buffers, std::string_view what)
{
  for (cl_uint k = 0; k < buffers.size(); ++k)
  {
    check(clSetKernelArg(kernel, k, sizeof(cl_mem), &buffers[k]), what);
  }
}

/// Fills a 9-entry column-major (ax ay az bx by bz cx cy cz) float array from the cell matrix columns.
template <typename Matrix>
inline void matrixToFloats(const Matrix& m, float* out)
{
  out[0] = static_cast<float>(m.ax);
  out[1] = static_cast<float>(m.ay);
  out[2] = static_cast<float>(m.az);
  out[3] = static_cast<float>(m.bx);
  out[4] = static_cast<float>(m.by);
  out[5] = static_cast<float>(m.bz);
  out[6] = static_cast<float>(m.cx);
  out[7] = static_cast<float>(m.cy);
  out[8] = static_cast<float>(m.cz);
}
}  // namespace OpenCLDevice
