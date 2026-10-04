module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module spatial_decomposition_opencl_context;

import std;

import opencl;
import spatial_decomposition_device_context;
import spatial_decomposition_opencl_handles;

using OpenCLDevice::check;

namespace
{
/**
 * \brief DeviceContext on OpenCL 1.2: one in-order command queue; buffers are cl_mem (Shared ones allocated with
 * CL_MEM_ALLOC_HOST_PTR and mapped blocking), kernels are cl_kernel objects of cached programs, marks are the
 * events of enqueued markers.
 */
class OpenCLContext final : public DeviceContext
{
 public:
  OpenCLContext();
  ~OpenCLContext() override;

  std::string deviceName() const override;
  std::size_t localMemorySize() const override { return localMemory; }

  DeviceBuffer createBuffer(std::size_t bytes, DeviceMemory memory) override;
  void releaseBuffer(DeviceBuffer& buffer) override;
  void* map(DeviceBuffer buffer, std::size_t bytes, bool forWriting, bool discard) override;
  void unmap(DeviceBuffer buffer) override;
  void readback(DeviceBuffer, std::size_t) override {}  // the mapping of the host-allocated buffer is the download
  void write(DeviceBuffer buffer, std::size_t offset, std::size_t bytes, const void* data, bool blocking) override;
  void read(DeviceBuffer buffer, std::size_t offset, std::size_t bytes, void* data) override;
  void copy(DeviceBuffer source, std::size_t sourceOffset, DeviceBuffer destination, std::size_t destinationOffset,
            std::size_t bytes) override;

  DeviceKernel compileKernel(std::string_view program, const char* source, DeviceMath math,
                             std::string_view kernelName) override;
  std::size_t maxGroupSize(DeviceKernel kernel) const override;
  void launch(DeviceKernel kernel, std::span<const DeviceArg> arguments, std::size_t groups,
              std::size_t groupSize) override;

  DeviceEvent mark() override;
  void wait(DeviceEvent event) override;
  void flush() override;
  void finish() override;

 private:
  struct Buffer
  {
    OpenCLDevice::MemHandle memory{};
    std::size_t bytes{0};
    void* mapped{nullptr};
  };
  struct Kernel
  {
    OpenCLDevice::KernelHandle kernel{};
    std::string name{};
  };

  cl_context context{nullptr};
  cl_device_id device{nullptr};
  OpenCLDevice::QueueHandle queue{};
  std::size_t localMemory{32768};
  std::vector<Buffer> buffers{};  // index = id - 1
  std::vector<Kernel> kernels{};  // index = id - 1
  std::map<std::string, OpenCLDevice::ProgramHandle, std::less<>> programs{};
  std::uint64_t marks{0};
  std::deque<std::pair<std::uint64_t, OpenCLDevice::EventHandle>> pending{};

  Buffer& at(DeviceBuffer buffer)
  {
    if (buffer.id == 0 || buffer.id > buffers.size() || buffers[buffer.id - 1].memory.get() == nullptr)
    {
      throw std::runtime_error("[OpenCL]: invalid device buffer\n");
    }
    return buffers[buffer.id - 1];
  }
  const Kernel& at(DeviceKernel kernel) const
  {
    if (kernel.id == 0 || kernel.id > kernels.size())
    {
      throw std::runtime_error("[OpenCL]: invalid device kernel\n");
    }
    return kernels[kernel.id - 1];
  }
};

OpenCLContext::OpenCLContext()
{
  if (!openclAvailable())
  {
    throw std::runtime_error("[OpenCL pair kernel]: no OpenCL device is available\n");
  }
  cl_int error = CL_SUCCESS;
  context = OpenCL::clContext.value();
  device = OpenCL::clDeviceId.value();
  queue.reset(clCreateCommandQueue(context, device, 0, &error));
  check(error, "clCreateCommandQueue");
  cl_ulong local = 0;
  if (clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local), &local, nullptr) == CL_SUCCESS && local > 0)
  {
    localMemory = static_cast<std::size_t>(local);
  }
}

OpenCLContext::~OpenCLContext()
{
  if (queue.get() == nullptr) return;
  for (Buffer& buffer : buffers)
  {
    if (buffer.mapped != nullptr && buffer.memory.get() != nullptr)
    {
      clEnqueueUnmapMemObject(queue.get(), buffer.memory.get(), buffer.mapped, 0, nullptr, nullptr);
      buffer.mapped = nullptr;
    }
  }
  clFinish(queue.get());
}

std::string OpenCLContext::deviceName() const { return openclDeviceName(); }

DeviceBuffer OpenCLContext::createBuffer(std::size_t bytes, DeviceMemory memory)
{
  const cl_mem_flags flags = memory == DeviceMemory::Shared ? (CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE)
                                                            : CL_MEM_READ_WRITE;
  Buffer buffer{};
  buffer.memory.reset(OpenCL::createBuffer(flags, std::max<std::size_t>(bytes, 1)));
  buffer.bytes = bytes;
  // reuse a released slot
  for (std::size_t k = 0; k < buffers.size(); ++k)
  {
    if (buffers[k].memory.get() == nullptr)
    {
      buffers[k] = std::move(buffer);
      return DeviceBuffer{static_cast<std::uint32_t>(k + 1)};
    }
  }
  buffers.push_back(std::move(buffer));
  return DeviceBuffer{static_cast<std::uint32_t>(buffers.size())};
}

void OpenCLContext::releaseBuffer(DeviceBuffer& handle)
{
  if (handle.id == 0) return;
  Buffer& buffer = at(handle);
  if (buffer.mapped != nullptr)
  {
    check(clEnqueueUnmapMemObject(queue.get(), buffer.memory.get(), buffer.mapped, 0, nullptr, nullptr),
          "clEnqueueUnmapMemObject");
    buffer.mapped = nullptr;
  }
  // the queue may still reference the buffer: the release waits for the enqueued work on it
  buffer.memory.reset();
  buffer.bytes = 0;
  handle = DeviceBuffer{};
}

void* OpenCLContext::map(DeviceBuffer handle, std::size_t bytes, bool forWriting, bool /*discard*/)
{
  // no CL_MAP_WRITE_INVALIDATE_REGION for `discard`: that leaves the mapped contents undefined, while the
  // callers rely on what they wrote earlier (the dummy slots) staying in place
  Buffer& buffer = at(handle);
  if (buffer.mapped != nullptr) unmap(handle);
  cl_int error = CL_SUCCESS;
  buffer.mapped = clEnqueueMapBuffer(queue.get(), buffer.memory.get(), CL_TRUE, forWriting ? CL_MAP_WRITE : CL_MAP_READ,
                                     0, std::max<std::size_t>(bytes, 1), 0, nullptr, nullptr, &error);
  check(error, "clEnqueueMapBuffer");
  return buffer.mapped;
}

void OpenCLContext::unmap(DeviceBuffer handle)
{
  Buffer& buffer = at(handle);
  if (buffer.mapped == nullptr) return;
  check(clEnqueueUnmapMemObject(queue.get(), buffer.memory.get(), buffer.mapped, 0, nullptr, nullptr),
        "clEnqueueUnmapMemObject");
  buffer.mapped = nullptr;
}

void OpenCLContext::write(DeviceBuffer handle, std::size_t offset, std::size_t bytes, const void* data, bool blocking)
{
  if (bytes == 0) return;
  check(clEnqueueWriteBuffer(queue.get(), at(handle).memory.get(), blocking ? CL_TRUE : CL_FALSE, offset, bytes,
                             data, 0, nullptr, nullptr),
        "clEnqueueWriteBuffer");
}

void OpenCLContext::read(DeviceBuffer handle, std::size_t offset, std::size_t bytes, void* data)
{
  if (bytes == 0) return;
  check(clEnqueueReadBuffer(queue.get(), at(handle).memory.get(), CL_FALSE, offset, bytes, data, 0, nullptr, nullptr),
        "clEnqueueReadBuffer");
}

void OpenCLContext::copy(DeviceBuffer source, std::size_t sourceOffset, DeviceBuffer destination,
                         std::size_t destinationOffset, std::size_t bytes)
{
  if (bytes == 0) return;
  check(clEnqueueCopyBuffer(queue.get(), at(source).memory.get(), at(destination).memory.get(), sourceOffset,
                            destinationOffset, bytes, 0, nullptr, nullptr),
        "clEnqueueCopyBuffer");
}

DeviceKernel OpenCLContext::compileKernel(std::string_view programName, const char* source, DeviceMath math,
                                          std::string_view kernelName)
{
  auto found = programs.find(programName);
  if (found == programs.end())
  {
    const char* options = math == DeviceMath::Fast      ? "-cl-mad-enable -cl-no-signed-zeros -cl-fast-relaxed-math"
                          : math == DeviceMath::Relaxed ? "-cl-mad-enable -cl-no-signed-zeros"
                                                        : "";
    OpenCLDevice::ProgramHandle program = OpenCLDevice::buildDeviceProgram(
        context, device, source, options, std::format("OpenCL {}", programName));
    found = programs.emplace(std::string(programName), std::move(program)).first;
  }
  Kernel kernel{};
  kernel.name = std::string(kernelName);
  kernel.kernel = OpenCLDevice::createKernel(found->second.get(), kernel.name.c_str());
  kernels.push_back(std::move(kernel));
  return DeviceKernel{static_cast<std::uint32_t>(kernels.size())};
}

std::size_t OpenCLContext::maxGroupSize(DeviceKernel handle) const
{
  std::size_t limit = 0;
  if (clGetKernelWorkGroupInfo(at(handle).kernel.get(), device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(limit), &limit,
                               nullptr) != CL_SUCCESS)
  {
    return 0;
  }
  return limit;
}

void OpenCLContext::launch(DeviceKernel handle, std::span<const DeviceArg> arguments, std::size_t groups,
                           std::size_t groupSize)
{
  const Kernel& kernel = at(handle);
  for (cl_uint k = 0; k < arguments.size(); ++k)
  {
    const DeviceArg& argument = arguments[k];
    switch (argument.kind)
    {
      case DeviceArg::Kind::Buffer:
      {
        const cl_mem memory = at(argument.buffer).memory.get();
        check(clSetKernelArg(kernel.kernel.get(), k, sizeof(cl_mem), &memory),
              std::format("clSetKernelArg ({}, {})", kernel.name, k));
        break;
      }
      case DeviceArg::Kind::Value:
        check(clSetKernelArg(kernel.kernel.get(), k, argument.bytes, argument.data),
              std::format("clSetKernelArg ({}, {})", kernel.name, k));
        break;
      case DeviceArg::Kind::Local:
        check(clSetKernelArg(kernel.kernel.get(), k, argument.bytes, nullptr),
              std::format("clSetKernelArg ({}, {} local)", kernel.name, k));
        break;
    }
  }
  const std::size_t global = std::max<std::size_t>(groups, 1) * groupSize;
  check(clEnqueueNDRangeKernel(queue.get(), kernel.kernel.get(), 1, nullptr, &global, &groupSize, 0, nullptr, nullptr),
        std::format("clEnqueueNDRangeKernel ({})", kernel.name));
}

DeviceEvent OpenCLContext::mark()
{
  cl_event event = nullptr;
  check(clEnqueueMarker(queue.get(), &event), "clEnqueueMarker");
  ++marks;
  pending.emplace_back(marks, OpenCLDevice::EventHandle(event));
  return DeviceEvent{marks};
}

void OpenCLContext::wait(DeviceEvent event)
{
  if (event.sequence == 0) return;
  // the in-order queue: waiting for the mark completes everything enqueued before it
  while (!pending.empty() && pending.front().first <= event.sequence)
  {
    if (pending.front().first == event.sequence)
    {
      cl_event e = pending.front().second.get();
      check(clWaitForEvents(1, &e), "clWaitForEvents");
    }
    pending.pop_front();
  }
}

void OpenCLContext::flush() { check(clFlush(queue.get()), "clFlush"); }

void OpenCLContext::finish()
{
  check(clFinish(queue.get()), "clFinish");
  pending.clear();
}
}  // namespace

bool openclAvailable()
{
  OpenCL::initialize();
  return OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value();
}

std::string openclDeviceName()
{
  if (!openclAvailable()) return {};
  char name[256] = {};
  std::size_t length = 0;
  if (clGetDeviceInfo(OpenCL::clDeviceId.value(), CL_DEVICE_NAME, sizeof(name) - 1, name, &length) != CL_SUCCESS)
  {
    return "unknown device";
  }
  return std::string(name);
}

std::unique_ptr<DeviceContext> createOpenCLContext() { return std::make_unique<OpenCLContext>(); }
