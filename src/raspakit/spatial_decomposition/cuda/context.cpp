module;

#include "api.h"

module spatial_decomposition_cuda_context;

import std;

import spatial_decomposition_device_context;
import spatial_decomposition_cuda_api;

namespace
{
constexpr std::size_t defaultDynamicShared = 48 * 1024;  // dynamic shared memory allowed without the opt-in
constexpr std::size_t sharedAlignment = 16;

std::size_t roundUp(std::size_t value, std::size_t multiple) { return ((value + multiple - 1) / multiple) * multiple; }

void check(const CUDAApi& api, CUresult result, std::string_view what)
{
  if (result != RASPA_CUDA_SUCCESS)
  {
    throw std::runtime_error(std::format("[CUDA]: {} failed ({})\n", what, api.errorString(result)));
  }
}

/// The device the backend runs on: the one with the most multiprocessors (ties: the lowest ordinal), so that
/// CUDA_VISIBLE_DEVICES still selects among several.
struct SelectedDevice
{
  CUdevice device{0};
  std::string name{};
  int major{0};
  int minor{0};
  std::size_t sharedPerBlock{defaultDynamicShared};
  std::size_t sharedOptIn{defaultDynamicShared};
  bool memoryPools{false};
};

std::optional<SelectedDevice> selectDevice(const CUDAApi& api, std::string& reason)
{
  int count = 0;
  CUresult result = api.cuDeviceGetCount(&count);
  if (result != RASPA_CUDA_SUCCESS)
  {
    reason = std::format("cuDeviceGetCount failed ({})", api.errorString(result));
    return std::nullopt;
  }
  if (count <= 0)
  {
    reason = "no CUDA device is present";
    return std::nullopt;
  }
  std::optional<SelectedDevice> best;
  int bestMultiprocessors = -1;
  for (int ordinal = 0; ordinal < count; ++ordinal)
  {
    CUdevice device = 0;
    if (api.cuDeviceGet(&device, ordinal) != RASPA_CUDA_SUCCESS) continue;
    int multiprocessors = 0;
    api.cuDeviceGetAttribute(&multiprocessors, RASPA_CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, device);
    if (multiprocessors <= bestMultiprocessors) continue;
    SelectedDevice selected{};
    selected.device = device;
    char name[256] = {};
    api.cuDeviceGetName(name, sizeof(name) - 1, device);
    selected.name = name;
    api.cuDeviceGetAttribute(&selected.major, RASPA_CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, device);
    api.cuDeviceGetAttribute(&selected.minor, RASPA_CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, device);
    int value = 0;
    if (api.cuDeviceGetAttribute(&value, RASPA_CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, device) ==
            RASPA_CUDA_SUCCESS &&
        value > 0)
    {
      selected.sharedPerBlock = static_cast<std::size_t>(value);
    }
    value = 0;
    if (api.cuDeviceGetAttribute(&value, RASPA_CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK_OPTIN, device) ==
            RASPA_CUDA_SUCCESS &&
        value > 0)
    {
      selected.sharedOptIn = std::max(selected.sharedPerBlock, static_cast<std::size_t>(value));
    }
    else
    {
      selected.sharedOptIn = selected.sharedPerBlock;
    }
    value = 0;
    if (api.cuDeviceGetAttribute(&value, RASPA_CU_DEVICE_ATTRIBUTE_MEMORY_POOLS_SUPPORTED, device) ==
        RASPA_CUDA_SUCCESS)
    {
      selected.memoryPools = value != 0;
    }
    best = selected;
    bestMultiprocessors = multiprocessors;
  }
  if (!best.has_value()) reason = "no CUDA device could be queried";
  return best;
}

struct Selection
{
  const CUDAApi* api{nullptr};
  std::optional<SelectedDevice> device{};
  std::string reason{};
};

const Selection& selection()
{
  static const Selection instance = []
  {
    Selection s{};
    s.api = cudaApi();
    if (s.api == nullptr)
    {
      s.reason = cudaApiError();
      return s;
    }
    s.device = selectDevice(*s.api, s.reason);
    return s;
  }();
  return instance;
}

/// Compiled programs shared by the contexts of the process (the tests create many contexts): keyed by the
/// program name, the source hash, the options and the architecture.
struct CompiledImage
{
  std::vector<char> image{};
  bool cubin{false};
};

std::mutex& imageCacheMutex()
{
  static std::mutex mutex;
  return mutex;
}
std::map<std::string, std::shared_ptr<const CompiledImage>>& imageCache()
{
  static std::map<std::string, std::shared_ptr<const CompiledImage>> cache;
  return cache;
}

/**
 * \brief DeviceContext on CUDA through the driver API and NVRTC (api.ixx): one non-blocking stream on the primary
 * context of the device, kernels from modules compiled at run time, marks as recorded events.
 *
 * Device buffers are stream-ordered allocations (cuMemAllocAsync where the driver has memory pools, so that
 * neither allocation nor release synchronizes). A Shared buffer pairs the device allocation with a pinned host
 * mirror: map completes the enqueued work and downloads the device copy into the mirror, unmap after a writing
 * map uploads the mirror in stream order. The two synchronizations of a step that this would cost are avoided
 * where the caller allows it: a map with `discard` (the host rewrites the range, or the device never wrote the
 * buffer) skips the download and only waits for the upload of the previous unmap when a mark taken after it has
 * not been waited for yet; a `readback` enqueued at the end of a chain downloads a buffer into its mirror in
 * stream order, so the map for reading after the wait for the chain's mark finds the data without a copy or a
 * synchronization (a launch with the buffer among its arguments, a write or a copy into it discards the
 * readback). Reads and non-blocking writes go through a pool of pinned staging buffers
 * (a copy from or to pageable memory would make the driver synchronize the stream): a write copies the data into
 * a staging buffer and uploads it asynchronously, a read downloads into one and is completed (copied to its
 * destination) when the mark taken after it is waited for, so that the host keeps running while the device works.
 *
 * NVRTC compiles for the exact architecture of the device when it knows it (CUBIN), else for the newest
 * architecture it supports below the device's (PTX, finalized by the driver): an NVRTC older than the GPU still
 * works.
 */
class CUDAContext final : public DeviceContext
{
 public:
  CUDAContext();
  ~CUDAContext() override;

  std::string deviceName() const override { return info.name; }
  std::size_t localMemorySize() const override { return info.sharedOptIn; }

  DeviceBuffer createBuffer(std::size_t bytes, DeviceMemory memory) override;
  void releaseBuffer(DeviceBuffer& buffer) override;
  void* map(DeviceBuffer buffer, std::size_t bytes, bool forWriting, bool discard) override;
  void unmap(DeviceBuffer buffer) override;
  void readback(DeviceBuffer buffer, std::size_t bytes) override;
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
    CUdeviceptr device{0};
    void* host{nullptr};  ///< pinned mirror of a Shared buffer
    std::size_t bytes{0};
    std::size_t mappedBytes{0};
    bool mappedForWriting{false};
    /// The mark counter when the upload of the mirror (unmap after a writing map) was enqueued: the host may
    /// write the mirror again once a later mark has been waited for.
    std::uint64_t uploadAtMark{0};
    /// A download into the mirror enqueued by `readback` (bytes, and the mark counter at the time); 0 when none
    /// is pending or when device work enqueued after it may have changed the buffer.
    std::size_t readbackBytes{0};
    std::uint64_t readbackAtMark{0};
  };
  struct Staging
  {
    void* host{nullptr};
    std::size_t bytes{0};
  };
  /// A staged transfer: a read to complete into `destination`, or a write (destination null) whose staging buffer
  /// returns to the pool; done once a mark taken after it has been waited for (the stream passed the transfer).
  struct PendingTransfer
  {
    Staging staging{};
    std::size_t bytes{0};
    void* destination{nullptr};
    std::uint64_t enqueuedAtMark{0};  ///< the mark counter when enqueued
  };
  struct Kernel
  {
    CUfunction function{nullptr};
    std::string name{};
    std::size_t maxThreads{1024};
    std::size_t staticShared{0};
    std::size_t dynamicLimit{defaultDynamicShared};  ///< dynamic shared memory the function may currently use
  };

  const CUDAApi& api;
  SelectedDevice info{};
  CUcontext context{nullptr};
  CUstream stream{nullptr};
  std::vector<Buffer> buffers{};  // index = id - 1
  std::vector<Kernel> kernels{};  // index = id - 1
  std::map<std::string, CUmodule, std::less<>> modules{};
  std::vector<Staging> stagingPool{};
  std::deque<PendingTransfer> transfers{};
  std::uint64_t marks{0};
  std::uint64_t completedMark{0};  ///< the highest mark sequence waited for: the stream has passed everything before it
  std::deque<std::pair<std::uint64_t, CUevent>> pending{};
  std::vector<CUevent> eventPool{};
  std::string architectureOption{};
  bool emitCubin{false};

  void bind() const { api.cuCtxSetCurrent(context); }
  Buffer& at(DeviceBuffer buffer)
  {
    if (buffer.id == 0 || buffer.id > buffers.size() || buffers[buffer.id - 1].device == 0)
    {
      throw std::runtime_error("[CUDA]: invalid device buffer\n");
    }
    return buffers[buffer.id - 1];
  }
  const Kernel& at(DeviceKernel kernel) const
  {
    if (kernel.id == 0 || kernel.id > kernels.size()) throw std::runtime_error("[CUDA]: invalid device kernel\n");
    return kernels[kernel.id - 1];
  }
  CUdeviceptr allocateDevice(std::size_t bytes);
  void freeDevice(CUdeviceptr pointer);
  Staging takeStaging(std::size_t bytes);
  void completeTransfers(std::uint64_t beforeMark);
  void chooseArchitecture();
  std::shared_ptr<const CompiledImage> compileProgram(std::string_view programName, const std::string& source,
                                                      DeviceMath math);
  CUmodule loadModule(std::string_view programName, const CompiledImage& image);
};

const CUDAApi& requireApi()
{
  const Selection& s = selection();
  if (s.api == nullptr || !s.device.has_value())
  {
    throw std::runtime_error(std::format("[CUDA pair kernel]: {}\n", s.reason));
  }
  return *s.api;
}

CUDAContext::CUDAContext() : api(requireApi())
{
  info = selection().device.value();
  check(api, api.cuDevicePrimaryCtxRetain(&context, info.device), "cuDevicePrimaryCtxRetain");
  bind();
  const CUresult result = api.cuStreamCreate(&stream, RASPA_CU_STREAM_NON_BLOCKING);
  if (result != RASPA_CUDA_SUCCESS)
  {
    api.cuDevicePrimaryCtxRelease(info.device);
    context = nullptr;
    check(api, result, "cuStreamCreate");
  }
  chooseArchitecture();
}

CUDAContext::~CUDAContext()
{
  if (context == nullptr) return;
  bind();
  try
  {
    finish();
  }
  catch (...)
  {
  }
  for (Buffer& buffer : buffers)
  {
    if (buffer.device != 0) freeDevice(buffer.device);
    if (buffer.host != nullptr) api.cuMemFreeHost(buffer.host);
  }
  api.cuStreamSynchronize(stream);
  for (PendingTransfer& transfer : transfers) api.cuMemFreeHost(transfer.staging.host);
  for (Staging& staging : stagingPool) api.cuMemFreeHost(staging.host);
  for (auto& [sequence, event] : pending) api.cuEventDestroy(event);
  for (CUevent event : eventPool) api.cuEventDestroy(event);
  for (auto& [name, module] : modules) api.cuModuleUnload(module);
  if (stream != nullptr) api.cuStreamDestroy(stream);
  api.cuDevicePrimaryCtxRelease(info.device);
}

void CUDAContext::chooseArchitecture()
{
  const int deviceArch = info.major * 10 + info.minor;
  int chosen = deviceArch;
  emitCubin = true;
  if (api.nvrtcGetNumSupportedArchs != nullptr && api.nvrtcGetSupportedArchs != nullptr)
  {
    int count = 0;
    if (api.nvrtcGetNumSupportedArchs(&count) == RASPA_NVRTC_SUCCESS && count > 0)
    {
      std::vector<int> archs(static_cast<std::size_t>(count));
      if (api.nvrtcGetSupportedArchs(archs.data()) == RASPA_NVRTC_SUCCESS)
      {
        if (std::find(archs.begin(), archs.end(), deviceArch) == archs.end())
        {
          // NVRTC predates the GPU: PTX of the newest architecture it knows below the device's, finalized by the
          // driver
          emitCubin = false;
          chosen = 0;
          for (int arch : archs)
          {
            if (arch < deviceArch && arch > chosen) chosen = arch;
          }
          if (chosen == 0)
          {
            throw std::runtime_error(
                std::format("[CUDA]: NVRTC supports no architecture at or below the device's sm_{}\n", deviceArch));
          }
        }
      }
    }
  }
  architectureOption = std::format("--gpu-architecture={}_{}", emitCubin ? "sm" : "compute", chosen);
}

CUdeviceptr CUDAContext::allocateDevice(std::size_t bytes)
{
  CUdeviceptr pointer = 0;
  if (info.memoryPools && api.cuMemAllocAsync != nullptr)
  {
    check(api, api.cuMemAllocAsync(&pointer, bytes, stream), std::format("cuMemAllocAsync ({} bytes)", bytes));
  }
  else
  {
    check(api, api.cuMemAlloc(&pointer, bytes), std::format("cuMemAlloc ({} bytes)", bytes));
  }
  return pointer;
}

void CUDAContext::freeDevice(CUdeviceptr pointer)
{
  if (info.memoryPools && api.cuMemFreeAsync != nullptr)
  {
    api.cuMemFreeAsync(pointer, stream);
  }
  else
  {
    api.cuMemFree(pointer);  // synchronizes the device implicitly
  }
}

DeviceBuffer CUDAContext::createBuffer(std::size_t bytes, DeviceMemory memory)
{
  bind();
  const std::size_t allocated = std::max<std::size_t>(bytes, 1);
  Buffer buffer{};
  buffer.bytes = bytes;
  buffer.device = allocateDevice(allocated);
  if (memory == DeviceMemory::Shared)
  {
    const CUresult result = api.cuMemHostAlloc(&buffer.host, allocated, RASPA_CU_MEMHOSTALLOC_PORTABLE);
    if (result != RASPA_CUDA_SUCCESS)
    {
      freeDevice(buffer.device);
      check(api, result, std::format("cuMemHostAlloc ({} bytes)", bytes));
    }
  }
  for (std::size_t k = 0; k < buffers.size(); ++k)
  {
    if (buffers[k].device == 0)
    {
      buffers[k] = buffer;
      return DeviceBuffer{static_cast<std::uint32_t>(k + 1)};
    }
  }
  buffers.push_back(buffer);
  return DeviceBuffer{static_cast<std::uint32_t>(buffers.size())};
}

void CUDAContext::releaseBuffer(DeviceBuffer& handle)
{
  if (handle.id == 0) return;
  bind();
  Buffer& buffer = at(handle);
  freeDevice(buffer.device);  // stream-ordered: after the enqueued work on it
  if (buffer.host != nullptr)
  {
    // an upload from the mirror may still be in flight
    check(api, api.cuStreamSynchronize(stream), "cuStreamSynchronize");
    api.cuMemFreeHost(buffer.host);
  }
  buffer = Buffer{};
  handle = DeviceBuffer{};
}

void* CUDAContext::map(DeviceBuffer handle, std::size_t bytes, bool forWriting, bool discard)
{
  bind();
  Buffer& buffer = at(handle);
  if (buffer.host == nullptr) throw std::runtime_error("[CUDA]: only Shared buffers can be mapped\n");
  if (buffer.mappedBytes != 0) unmap(handle);
  const std::size_t length = std::min(std::max<std::size_t>(bytes, 1), std::max<std::size_t>(buffer.bytes, 1));
  if (discard)
  {
    // the mirror keeps what the host wrote last; it may be written again once the upload of the previous unmap
    // has passed (the kernels read the device copy, not the mirror)
    if (buffer.uploadAtMark >= completedMark) finish();
  }
  else if (!forWriting && buffer.readbackBytes >= length && buffer.readbackAtMark < completedMark)
  {
    // the readback enqueued in the chain has landed in the mirror: nothing to download, nothing to wait for
  }
  else
  {
    // the mapped memory shows the current contents also when mapped for writing (CL_MAP_WRITE, Metal shared
    // storage): the host may rewrite only part of it, the rest must stay what the device last wrote
    check(api, api.cuMemcpyDtoHAsync(buffer.host, buffer.device, length, stream), "cuMemcpyDtoHAsync (map)");
    finish();  // the host owns the memory from here: the enqueued work on the buffer has completed
  }
  buffer.readbackBytes = 0;
  buffer.mappedBytes = length;
  buffer.mappedForWriting = forWriting;
  return buffer.host;
}

void CUDAContext::unmap(DeviceBuffer handle)
{
  bind();
  Buffer& buffer = at(handle);
  if (buffer.mappedBytes == 0) return;
  if (buffer.mappedForWriting)
  {
    check(api, api.cuMemcpyHtoDAsync(buffer.device, buffer.host, buffer.mappedBytes, stream),
          "cuMemcpyHtoDAsync (unmap)");
    buffer.uploadAtMark = marks;
  }
  buffer.mappedBytes = 0;
  buffer.mappedForWriting = false;
}

void CUDAContext::readback(DeviceBuffer handle, std::size_t bytes)
{
  if (bytes == 0) return;
  bind();
  Buffer& buffer = at(handle);
  if (buffer.host == nullptr) throw std::runtime_error("[CUDA]: only Shared buffers can be read back\n");
  if (buffer.mappedBytes != 0) throw std::runtime_error("[CUDA]: a mapped buffer cannot be read back\n");
  const std::size_t length = std::min(bytes, std::max<std::size_t>(buffer.bytes, 1));
  check(api, api.cuMemcpyDtoHAsync(buffer.host, buffer.device, length, stream), "cuMemcpyDtoHAsync (readback)");
  buffer.readbackBytes = length;
  buffer.readbackAtMark = marks;
}

void CUDAContext::write(DeviceBuffer handle, std::size_t offset, std::size_t bytes, const void* data, bool blocking)
{
  if (bytes == 0) return;
  bind();
  Buffer& buffer = at(handle);
  if (offset + bytes > std::max<std::size_t>(buffer.bytes, 1))
  {
    throw std::runtime_error("[CUDA]: write beyond the end of a buffer\n");
  }
  buffer.readbackBytes = 0;
  if (blocking)
  {
    check(api, api.cuMemcpyHtoDAsync(buffer.device + offset, data, bytes, stream), "cuMemcpyHtoDAsync");
    check(api, api.cuStreamSynchronize(stream), "cuStreamSynchronize (write)");
    return;
  }
  PendingTransfer transfer{};
  transfer.staging = takeStaging(bytes);
  transfer.bytes = bytes;
  transfer.enqueuedAtMark = marks;
  std::memcpy(transfer.staging.host, data, bytes);
  const CUresult result = api.cuMemcpyHtoDAsync(buffer.device + offset, transfer.staging.host, bytes, stream);
  if (result != RASPA_CUDA_SUCCESS)
  {
    stagingPool.push_back(transfer.staging);
    check(api, result, "cuMemcpyHtoDAsync");
  }
  transfers.push_back(transfer);
}

CUDAContext::Staging CUDAContext::takeStaging(std::size_t bytes)
{
  // the smallest pooled buffer that fits
  std::size_t best = stagingPool.size();
  for (std::size_t k = 0; k < stagingPool.size(); ++k)
  {
    if (stagingPool[k].bytes >= bytes && (best == stagingPool.size() || stagingPool[k].bytes < stagingPool[best].bytes))
    {
      best = k;
    }
  }
  if (best < stagingPool.size())
  {
    Staging staging = stagingPool[best];
    stagingPool.erase(stagingPool.begin() + static_cast<std::ptrdiff_t>(best));
    return staging;
  }
  Staging staging{};
  staging.bytes = std::max<std::size_t>(bytes + bytes / 2, 4096);
  check(api, api.cuMemHostAlloc(&staging.host, staging.bytes, RASPA_CU_MEMHOSTALLOC_PORTABLE),
        "cuMemHostAlloc (staging)");
  return staging;
}

void CUDAContext::read(DeviceBuffer handle, std::size_t offset, std::size_t bytes, void* data)
{
  if (bytes == 0) return;
  bind();
  Buffer& buffer = at(handle);
  if (offset + bytes > std::max<std::size_t>(buffer.bytes, 1))
  {
    throw std::runtime_error("[CUDA]: read beyond the end of a buffer\n");
  }
  PendingTransfer transfer{};
  transfer.staging = takeStaging(bytes);
  transfer.bytes = bytes;
  transfer.destination = data;
  transfer.enqueuedAtMark = marks;
  const CUresult result = api.cuMemcpyDtoHAsync(transfer.staging.host, buffer.device + offset, bytes, stream);
  if (result != RASPA_CUDA_SUCCESS)
  {
    stagingPool.push_back(transfer.staging);
    check(api, result, "cuMemcpyDtoHAsync");
  }
  transfers.push_back(transfer);
}

void CUDAContext::completeTransfers(std::uint64_t beforeMark)
{
  while (!transfers.empty() && transfers.front().enqueuedAtMark < beforeMark)
  {
    PendingTransfer& transfer = transfers.front();
    if (transfer.destination != nullptr) std::memcpy(transfer.destination, transfer.staging.host, transfer.bytes);
    stagingPool.push_back(transfer.staging);
    transfers.pop_front();
  }
}

void CUDAContext::copy(DeviceBuffer source, std::size_t sourceOffset, DeviceBuffer destination,
                       std::size_t destinationOffset, std::size_t bytes)
{
  if (bytes == 0) return;
  bind();
  Buffer& from = at(source);
  Buffer& to = at(destination);
  if (sourceOffset + bytes > std::max<std::size_t>(from.bytes, 1) ||
      destinationOffset + bytes > std::max<std::size_t>(to.bytes, 1))
  {
    throw std::runtime_error("[CUDA]: copy beyond the end of a buffer\n");
  }
  to.readbackBytes = 0;
  check(api, api.cuMemcpyDtoDAsync(to.device + destinationOffset, from.device + sourceOffset, bytes, stream),
        "cuMemcpyDtoDAsync");
}

std::shared_ptr<const CompiledImage> CUDAContext::compileProgram(std::string_view programName,
                                                                 const std::string& source, DeviceMath math)
{
  // Strict: no contractions (the double-float arithmetic). Relaxed: contractions (the NVRTC default). Fast: also
  // denormals flushed to zero; not --use_fast_math, whose approximate division, square root and intrinsic
  // transcendentals lose more than the OpenCL backend does with -cl-fast-relaxed-math (the r^-12 terms of
  // overlapping pairs amplify the division error beyond the tolerance of the device tests).
  std::vector<const char*> options{architectureOption.c_str(), "-std=c++17", "-default-device"};
  if (math == DeviceMath::Strict) options.push_back("--fmad=false");
  if (math == DeviceMath::Fast) options.push_back("--ftz=true");

  std::string optionText;
  for (const char* option : options) optionText += std::string(option) + " ";
  const std::string key =
      std::format("{}|{:016x}|{}|{}", programName, std::hash<std::string>{}(source), optionText, emitCubin);
  {
    std::lock_guard<std::mutex> lock(imageCacheMutex());
    auto found = imageCache().find(key);
    if (found != imageCache().end()) return found->second;
  }

  nvrtcProgram program = nullptr;
  const std::string fileName = std::format("{}.cu", programName);
  nvrtcResult result = api.nvrtcCreateProgram(&program, source.c_str(), fileName.c_str(), 0, nullptr, nullptr);
  if (result != RASPA_NVRTC_SUCCESS)
  {
    throw std::runtime_error(
        std::format("[CUDA {}]: nvrtcCreateProgram failed ({})\n", programName, api.nvrtcGetErrorString(result)));
  }
  result = api.nvrtcCompileProgram(program, static_cast<int>(options.size()), options.data());
  std::string log;
  std::size_t logSize = 0;
  if (api.nvrtcGetProgramLogSize(program, &logSize) == RASPA_NVRTC_SUCCESS && logSize > 1)
  {
    log.resize(logSize);
    api.nvrtcGetProgramLog(program, log.data());
    log.resize(std::strlen(log.c_str()));
  }
  if (result != RASPA_NVRTC_SUCCESS)
  {
    api.nvrtcDestroyProgram(&program);
    throw std::runtime_error(
        std::format("[CUDA {}]: the kernels failed to build ({}, {}):\n{}\n", programName, architectureOption,
                    api.nvrtcGetErrorString(result), log));
  }
  auto compiled = std::make_shared<CompiledImage>();
  compiled->cubin = emitCubin;
  std::size_t size = 0;
  if (emitCubin)
  {
    api.nvrtcGetCUBINSize(program, &size);
    compiled->image.resize(size);
    api.nvrtcGetCUBIN(program, compiled->image.data());
  }
  else
  {
    api.nvrtcGetPTXSize(program, &size);
    compiled->image.resize(size);
    api.nvrtcGetPTX(program, compiled->image.data());
  }
  api.nvrtcDestroyProgram(&program);
  if (compiled->image.empty())
  {
    throw std::runtime_error(std::format("[CUDA {}]: NVRTC produced no code\n", programName));
  }

  std::lock_guard<std::mutex> lock(imageCacheMutex());
  return imageCache().emplace(key, compiled).first->second;
}

CUmodule CUDAContext::loadModule(std::string_view programName, const CompiledImage& image)
{
  char errorLog[4096] = {};
  char infoLog[4096] = {};
  CUjit_option optionNames[] = {RASPA_CU_JIT_ERROR_LOG_BUFFER, RASPA_CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES,
                                RASPA_CU_JIT_INFO_LOG_BUFFER, RASPA_CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES};
  void* optionValues[] = {errorLog, reinterpret_cast<void*>(static_cast<std::uintptr_t>(sizeof(errorLog))), infoLog,
                          reinterpret_cast<void*>(static_cast<std::uintptr_t>(sizeof(infoLog)))};
  CUmodule module = nullptr;
  const CUresult result = api.cuModuleLoadDataEx(&module, image.image.data(), 4, optionNames, optionValues);
  if (result != RASPA_CUDA_SUCCESS)
  {
    throw std::runtime_error(std::format("[CUDA {}]: the module could not be loaded ({}):\n{}\n", programName,
                                         api.errorString(result), errorLog));
  }
  return module;
}

DeviceKernel CUDAContext::compileKernel(std::string_view programName, const char* source, DeviceMath math,
                                        std::string_view kernelName)
{
  bind();
  auto found = modules.find(programName);
  if (found == modules.end())
  {
    const std::string text = std::string(cudaKernelDialect) + source;
    std::shared_ptr<const CompiledImage> image = compileProgram(programName, text, math);
    CUmodule module = loadModule(programName, *image);
    found = modules.emplace(std::string(programName), module).first;
  }
  Kernel kernel{};
  kernel.name = std::string(kernelName);
  check(api, api.cuModuleGetFunction(&kernel.function, found->second, kernel.name.c_str()),
        std::format("cuModuleGetFunction ({})", kernelName));
  int value = 0;
  if (api.cuFuncGetAttribute(&value, RASPA_CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK, kernel.function) ==
          RASPA_CUDA_SUCCESS &&
      value > 0)
  {
    kernel.maxThreads = static_cast<std::size_t>(value);
  }
  value = 0;
  if (api.cuFuncGetAttribute(&value, RASPA_CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES, kernel.function) ==
          RASPA_CUDA_SUCCESS &&
      value > 0)
  {
    kernel.staticShared = static_cast<std::size_t>(value);
  }
  const std::size_t base = std::min(info.sharedPerBlock, defaultDynamicShared);
  kernel.dynamicLimit = base > kernel.staticShared ? base - kernel.staticShared : 0;
  kernels.push_back(std::move(kernel));
  return DeviceKernel{static_cast<std::uint32_t>(kernels.size())};
}

std::size_t CUDAContext::maxGroupSize(DeviceKernel handle) const { return at(handle).maxThreads; }

void CUDAContext::launch(DeviceKernel handle, std::span<const DeviceArg> arguments, std::size_t groups,
                         std::size_t groupSize)
{
  bind();
  at(handle);  // validates the handle
  Kernel& kernel = kernels[handle.id - 1];
  if (groupSize == 0 || groupSize > kernel.maxThreads)
  {
    throw std::runtime_error(std::format("[CUDA]: kernel {} launched with a group of {} threads (at most {})\n",
                                         kernel.name, groupSize, kernel.maxThreads));
  }
  // the arguments in order; a local-memory argument becomes the byte offset of its slice of the dynamic shared
  // memory (LOCAL_ARG / LOCAL_ARG_BIND of the dialect)
  std::vector<void*> parameters;
  parameters.reserve(arguments.size());
  unsigned int localOffsets[16];
  std::size_t localCount = 0;
  std::size_t sharedBytes = 0;
  for (const DeviceArg& argument : arguments)
  {
    switch (argument.kind)
    {
      case DeviceArg::Kind::Buffer:
      {
        Buffer& buffer = at(argument.buffer);
        buffer.readbackBytes = 0;  // the kernel may write it
        parameters.push_back(&buffer.device);
        break;
      }
      case DeviceArg::Kind::Value:
        parameters.push_back(const_cast<void*>(argument.data));
        break;
      case DeviceArg::Kind::Local:
        if (localCount >= 16) throw std::runtime_error("[CUDA]: too many local-memory arguments\n");
        sharedBytes = roundUp(sharedBytes, sharedAlignment);
        localOffsets[localCount] = static_cast<unsigned int>(sharedBytes);
        parameters.push_back(&localOffsets[localCount]);
        ++localCount;
        sharedBytes += std::max<std::size_t>(argument.bytes, 1);
        break;
    }
  }
  sharedBytes = roundUp(sharedBytes, sharedAlignment);
  if (sharedBytes > kernel.dynamicLimit)
  {
    const std::size_t limit = info.sharedOptIn - std::min(kernel.staticShared, info.sharedOptIn);
    if (sharedBytes > limit)
    {
      throw std::runtime_error(std::format(
          "[CUDA]: kernel {} needs {} bytes of local memory, the device allows {}\n", kernel.name, sharedBytes, limit));
    }
    check(api,
          api.cuFuncSetAttribute(kernel.function, RASPA_CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES,
                                 static_cast<int>(sharedBytes)),
          std::format("cuFuncSetAttribute ({})", kernel.name));
    kernel.dynamicLimit = sharedBytes;
  }
  check(api,
        api.cuLaunchKernel(kernel.function, static_cast<unsigned int>(std::max<std::size_t>(groups, 1)), 1, 1,
                           static_cast<unsigned int>(groupSize), 1, 1, static_cast<unsigned int>(sharedBytes), stream,
                           parameters.data(), nullptr),
        std::format("cuLaunchKernel ({})", kernel.name));
}

DeviceEvent CUDAContext::mark()
{
  bind();
  CUevent event = nullptr;
  if (!eventPool.empty())
  {
    event = eventPool.back();
    eventPool.pop_back();
  }
  else
  {
    check(api, api.cuEventCreate(&event, RASPA_CU_EVENT_DISABLE_TIMING), "cuEventCreate");
  }
  const CUresult result = api.cuEventRecord(event, stream);
  if (result != RASPA_CUDA_SUCCESS)
  {
    eventPool.push_back(event);
    check(api, result, "cuEventRecord");
  }
  ++marks;
  pending.emplace_back(marks, event);
  return DeviceEvent{marks};
}

void CUDAContext::wait(DeviceEvent event)
{
  if (event.sequence == 0) return;
  bind();
  // the in-order stream: waiting for the mark completes everything enqueued before it
  while (!pending.empty() && pending.front().first <= event.sequence)
  {
    if (pending.front().first == event.sequence)
    {
      check(api, api.cuEventSynchronize(pending.front().second), "cuEventSynchronize");
    }
    eventPool.push_back(pending.front().second);
    pending.pop_front();
  }
  completedMark = std::max(completedMark, event.sequence);
  completeTransfers(event.sequence);
}

void CUDAContext::flush()
{
  // the driver submits the enqueued work as it comes
}

void CUDAContext::finish()
{
  bind();
  check(api, api.cuStreamSynchronize(stream), "cuStreamSynchronize");
  for (auto& [sequence, event] : pending) eventPool.push_back(event);
  pending.clear();
  completedMark = marks + 1;  // everything enqueued so far (at mark counters <= marks) has passed
  completeTransfers(marks + 1);
}
}  // namespace

bool cudaAvailable()
{
  const Selection& s = selection();
  return s.api != nullptr && s.device.has_value();
}

std::string cudaDeviceName()
{
  const Selection& s = selection();
  return s.device.has_value() ? s.device->name : std::string{};
}

std::string cudaUnavailableReason() { return selection().reason; }

std::unique_ptr<DeviceContext> createCUDAContext() { return std::make_unique<CUDAContext>(); }
