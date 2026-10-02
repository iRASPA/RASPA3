module;

#include "bridge.h"

module spatial_decomposition_metal_context;

import std;

import spatial_decomposition_device_context;

namespace
{
constexpr std::size_t errorCapacity = 4096;

/**
 * \brief DeviceContext on Metal through the C bridge (bridge.h / bridge.mm).
 *
 * All buffers live in shared storage (unified memory): a Shared buffer is mapped by returning its contents (the
 * map completes the in-flight work first; unmap is a no-op, the memory is coherent for command buffers committed
 * afterwards). Writes and reads are ordered in the stream through blit copies with staging buffers of a pool: a
 * write copies the host data into a staging buffer and blits it into the destination; a read blits into a
 * staging buffer and the host copy is made when the batch completed (wait / finish). A write to a buffer with
 * no work in flight copies directly.
 *
 * The stream is a sequence of batches (command buffers): commands are recorded into the open batch, flush()
 * commits it, mark() names the point in the stream, wait(mark) commits the open batch when it holds work from
 * before the mark and waits for every committed batch up to there, performing their reads and returning their
 * staging buffers. Released buffers are kept until no batch in flight can reference them (the command buffers
 * do not retain their resources).
 */
class MetalContext final : public DeviceContext
{
 public:
  MetalContext();
  ~MetalContext() override;

  std::string deviceName() const override { return name; }
  std::size_t localMemorySize() const override { return localMemory; }

  DeviceBuffer createBuffer(std::size_t bytes, DeviceMemory memory) override;
  void releaseBuffer(DeviceBuffer& buffer) override;
  void* map(DeviceBuffer buffer, std::size_t bytes, bool forWriting) override;
  void unmap(DeviceBuffer buffer) override;
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
    RaspaMetalBuffer* buffer{nullptr};
    std::size_t bytes{0};
  };
  struct Staging
  {
    RaspaMetalBuffer* buffer{nullptr};
    std::size_t bytes{0};
  };
  struct PendingRead
  {
    Staging staging{};
    std::size_t bytes{0};
    void* destination{nullptr};
  };
  struct Batch
  {
    RaspaMetalCommands* commands{nullptr};
    std::uint64_t startMark{0};  ///< the mark counter when the batch was opened
    bool hasWork{false};
    std::vector<Staging> stagings{};  ///< staging buffers of the writes (returned to the pool on completion)
    std::vector<PendingRead> reads{};
    std::vector<RaspaMetalBuffer*> graveyard{};  ///< buffers released while this batch may still use them
  };
  struct Kernel
  {
    RaspaMetalPipeline* pipeline{nullptr};
    std::string name{};
    std::size_t maxThreads{0};
  };

  RaspaMetalDevice* device{nullptr};
  std::string name{};
  std::size_t localMemory{32768};
  std::vector<Buffer> buffers{};
  std::vector<Kernel> kernels{};
  std::map<std::string, RaspaMetalLibrary*, std::less<>> libraries{};
  std::vector<Staging> stagingPool{};
  std::uint64_t marks{0};
  std::optional<Batch> open{};
  std::deque<Batch> committed{};

  Buffer& at(DeviceBuffer buffer)
  {
    if (buffer.id == 0 || buffer.id > buffers.size() || buffers[buffer.id - 1].buffer == nullptr)
    {
      throw std::runtime_error("[Metal]: invalid device buffer\n");
    }
    return buffers[buffer.id - 1];
  }
  const Kernel& at(DeviceKernel kernel) const
  {
    if (kernel.id == 0 || kernel.id > kernels.size()) throw std::runtime_error("[Metal]: invalid device kernel\n");
    return kernels[kernel.id - 1];
  }
  Batch& openBatch();
  void commitOpen();
  void completeFront();
  bool inFlight() const { return !committed.empty() || (open.has_value() && open->hasWork); }
  Staging takeStaging(std::size_t bytes);
  void destroyBatch(Batch& batch);
};

MetalContext::MetalContext()
{
  char error[errorCapacity] = {};
  device = raspaMetalCreateDevice(error, sizeof(error));
  if (device == nullptr)
  {
    throw std::runtime_error(std::format("[Metal pair kernel]: {}\n", error));
  }
  name = metalDeviceName();
  const std::size_t local = raspaMetalLocalMemorySize(device);
  if (local > 0) localMemory = local;
}

MetalContext::~MetalContext()
{
  if (device == nullptr) return;
  try
  {
    finish();
  }
  catch (...)
  {
  }
  if (open.has_value()) destroyBatch(*open);
  for (Kernel& kernel : kernels) raspaMetalDestroyPipeline(kernel.pipeline);
  for (auto& [key, library] : libraries) raspaMetalDestroyLibrary(library);
  for (Buffer& buffer : buffers) raspaMetalDestroyBuffer(buffer.buffer);
  for (Staging& staging : stagingPool) raspaMetalDestroyBuffer(staging.buffer);
  raspaMetalDestroyDevice(device);
}

void MetalContext::destroyBatch(Batch& batch)
{
  for (Staging& staging : batch.stagings) stagingPool.push_back(staging);
  for (PendingRead& read : batch.reads) stagingPool.push_back(read.staging);
  for (RaspaMetalBuffer* buffer : batch.graveyard) raspaMetalDestroyBuffer(buffer);
  if (batch.commands != nullptr) raspaMetalDestroyCommands(batch.commands);
  batch = Batch{};
}

MetalContext::Batch& MetalContext::openBatch()
{
  if (!open.has_value())
  {
    Batch batch{};
    batch.commands = raspaMetalBeginCommands(device);
    if (batch.commands == nullptr) throw std::runtime_error("[Metal]: no command buffer could be created\n");
    batch.startMark = marks;
    open = std::move(batch);
  }
  return *open;
}

void MetalContext::commitOpen()
{
  if (!open.has_value()) return;
  if (!open->hasWork && open->graveyard.empty())
  {
    destroyBatch(*open);
    open.reset();
    return;
  }
  raspaMetalCommit(open->commands);
  committed.push_back(std::move(*open));
  open.reset();
}

void MetalContext::completeFront()
{
  Batch& batch = committed.front();
  char error[errorCapacity] = {};
  const int status = raspaMetalWait(batch.commands, error, sizeof(error));
  if (status != 0)
  {
    const std::string message = error;
    destroyBatch(batch);
    committed.pop_front();
    throw std::runtime_error(std::format("[Metal]: {}\n", message));
  }
  for (PendingRead& read : batch.reads)
  {
    std::memcpy(read.destination, raspaMetalBufferContents(read.staging.buffer), read.bytes);
  }
  destroyBatch(batch);
  committed.pop_front();
}

MetalContext::Staging MetalContext::takeStaging(std::size_t bytes)
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
  staging.buffer = raspaMetalCreateBuffer(device, staging.bytes);
  if (staging.buffer == nullptr) throw std::runtime_error("[Metal]: a staging buffer could not be allocated\n");
  return staging;
}

DeviceBuffer MetalContext::createBuffer(std::size_t bytes, DeviceMemory)
{
  Buffer buffer{};
  buffer.buffer = raspaMetalCreateBuffer(device, std::max<std::size_t>(bytes, 1));
  if (buffer.buffer == nullptr)
  {
    throw std::runtime_error(std::format("[Metal]: a buffer of {} bytes could not be allocated\n", bytes));
  }
  buffer.bytes = bytes;
  for (std::size_t k = 0; k < buffers.size(); ++k)
  {
    if (buffers[k].buffer == nullptr)
    {
      buffers[k] = buffer;
      return DeviceBuffer{static_cast<std::uint32_t>(k + 1)};
    }
  }
  buffers.push_back(buffer);
  return DeviceBuffer{static_cast<std::uint32_t>(buffers.size())};
}

void MetalContext::releaseBuffer(DeviceBuffer& handle)
{
  if (handle.id == 0) return;
  Buffer& buffer = at(handle);
  if (inFlight())
  {
    // destroyed when the batches that may use it completed: appended to the open batch (completed last)
    openBatch().graveyard.push_back(buffer.buffer);
  }
  else
  {
    raspaMetalDestroyBuffer(buffer.buffer);
  }
  buffer = Buffer{};
  handle = DeviceBuffer{};
}

void* MetalContext::map(DeviceBuffer handle, std::size_t, bool)
{
  Buffer& buffer = at(handle);
  finish();
  return raspaMetalBufferContents(buffer.buffer);
}

void MetalContext::unmap(DeviceBuffer) {}

void MetalContext::write(DeviceBuffer handle, std::size_t offset, std::size_t bytes, const void* data, bool blocking)
{
  if (bytes == 0) return;
  Buffer& buffer = at(handle);
  if (offset + bytes > std::max<std::size_t>(buffer.bytes, 1))
  {
    throw std::runtime_error("[Metal]: write beyond the end of a buffer\n");
  }
  if (blocking || !inFlight())
  {
    if (blocking) finish();
    std::memcpy(static_cast<std::byte*>(raspaMetalBufferContents(buffer.buffer)) + offset, data, bytes);
    return;
  }
  Batch& batch = openBatch();
  Staging staging = takeStaging(bytes);
  std::memcpy(raspaMetalBufferContents(staging.buffer), data, bytes);
  raspaMetalCopy(batch.commands, staging.buffer, 0, buffer.buffer, offset, bytes);
  batch.stagings.push_back(staging);
  batch.hasWork = true;
}

void MetalContext::read(DeviceBuffer handle, std::size_t offset, std::size_t bytes, void* data)
{
  if (bytes == 0) return;
  Buffer& buffer = at(handle);
  Batch& batch = openBatch();
  PendingRead pending{};
  pending.staging = takeStaging(bytes);
  pending.bytes = bytes;
  pending.destination = data;
  raspaMetalCopy(batch.commands, buffer.buffer, offset, pending.staging.buffer, 0, bytes);
  batch.reads.push_back(pending);
  batch.hasWork = true;
}

void MetalContext::copy(DeviceBuffer source, std::size_t sourceOffset, DeviceBuffer destination,
                        std::size_t destinationOffset, std::size_t bytes)
{
  if (bytes == 0) return;
  Buffer& from = at(source);
  Buffer& to = at(destination);
  if (sourceOffset + bytes > std::max<std::size_t>(from.bytes, 1) ||
      destinationOffset + bytes > std::max<std::size_t>(to.bytes, 1))
  {
    throw std::runtime_error("[Metal]: copy beyond the end of a buffer\n");
  }
  Batch& batch = openBatch();
  raspaMetalCopy(batch.commands, from.buffer, sourceOffset, to.buffer, destinationOffset, bytes);
  batch.hasWork = true;
}

DeviceKernel MetalContext::compileKernel(std::string_view programName, const char* source, DeviceMath math,
                                         std::string_view kernelName)
{
  auto found = libraries.find(programName);
  if (found == libraries.end())
  {
    const std::string text = std::string(metalKernelDialect) + source;
    char error[errorCapacity] = {};
    RaspaMetalLibrary* library = raspaMetalCompile(device, text.c_str(), static_cast<int>(math), error, sizeof(error));
    if (library == nullptr)
    {
      throw std::runtime_error(std::format("[Metal {}]: the kernels failed to build:\n{}\n", programName, error));
    }
    found = libraries.emplace(std::string(programName), library).first;
  }
  Kernel kernel{};
  kernel.name = std::string(kernelName);
  char error[errorCapacity] = {};
  kernel.pipeline = raspaMetalCreatePipeline(device, found->second, kernel.name.c_str(), error, sizeof(error));
  if (kernel.pipeline == nullptr)
  {
    throw std::runtime_error(std::format("[Metal {}]: {}\n", programName, error));
  }
  kernel.maxThreads = raspaMetalMaxThreadsPerThreadgroup(kernel.pipeline);
  kernels.push_back(std::move(kernel));
  return DeviceKernel{static_cast<std::uint32_t>(kernels.size())};
}

std::size_t MetalContext::maxGroupSize(DeviceKernel handle) const { return at(handle).maxThreads; }

void MetalContext::launch(DeviceKernel handle, std::span<const DeviceArg> arguments, std::size_t groups,
                          std::size_t groupSize)
{
  const Kernel& kernel = at(handle);
  if (groupSize > kernel.maxThreads)
  {
    throw std::runtime_error(std::format("[Metal]: kernel {} launched with a group of {} threads (at most {})\n",
                                         kernel.name, groupSize, kernel.maxThreads));
  }
  RaspaMetalBuffer* bufferHandles[32];
  std::uint32_t bufferIndices[32];
  const void* values[32];
  std::size_t valueBytes[32];
  std::uint32_t valueIndices[32];
  std::size_t localBytes[8];
  std::size_t bufferCount = 0, valueCount = 0, localCount = 0;
  std::uint32_t slot = 0;
  for (const DeviceArg& argument : arguments)
  {
    switch (argument.kind)
    {
      case DeviceArg::Kind::Buffer:
        if (bufferCount >= 32) throw std::runtime_error("[Metal]: too many buffer arguments\n");
        bufferHandles[bufferCount] = at(argument.buffer).buffer;
        bufferIndices[bufferCount] = slot++;
        ++bufferCount;
        break;
      case DeviceArg::Kind::Value:
        if (valueCount >= 32) throw std::runtime_error("[Metal]: too many value arguments\n");
        values[valueCount] = argument.data;
        valueBytes[valueCount] = argument.bytes;
        valueIndices[valueCount] = slot++;
        ++valueCount;
        break;
      case DeviceArg::Kind::Local:
        if (localCount >= 8) throw std::runtime_error("[Metal]: too many local-memory arguments\n");
        localBytes[localCount++] = argument.bytes;
        break;
    }
  }
  Batch& batch = openBatch();
  raspaMetalDispatch(batch.commands, kernel.pipeline, bufferHandles, bufferIndices, bufferCount, values, valueBytes,
                     valueIndices, valueCount, localBytes, localCount, groups, groupSize);
  batch.hasWork = true;
}

DeviceEvent MetalContext::mark()
{
  ++marks;
  return DeviceEvent{marks};
}

void MetalContext::wait(DeviceEvent event)
{
  if (event.sequence == 0) return;
  // the open batch holds work from before the mark when it was opened before the mark was taken
  if (open.has_value() && open->startMark < event.sequence) commitOpen();
  while (!committed.empty() && committed.front().startMark < event.sequence) completeFront();
}

void MetalContext::flush() { commitOpen(); }

void MetalContext::finish()
{
  commitOpen();
  while (!committed.empty()) completeFront();
}
}  // namespace

bool metalAvailable() { return raspaMetalAvailable() != 0; }

std::string metalDeviceName()
{
  char name[256] = {};
  const std::size_t length = raspaMetalDeviceName(name, sizeof(name));
  return std::string(name, length);
}

std::unique_ptr<DeviceContext> createMetalContext() { return std::make_unique<MetalContext>(); }
