module;

export module spatial_decomposition_device_context;

import std;

/**
 * \file context.ixx
 * \brief The device API behind the spatial-decomposition device code: buffers, kernels, launches, transfers and
 * completion, in the common denominator of OpenCL, Metal and CUDA.
 *
 * Everything above this interface (DeviceStep, DeviceMesh, DeviceBonded, the resident integrator) is written
 * once against it and runs on every backend; the kernels are the shared sources of
 * spatial_decomposition_device_kernels, compiled by the context with its dialect header.
 *
 * Conventions:
 *  - One in-order stream: writes, launches and reads complete in the order they were enqueued. All calls are made
 *    from one thread.
 *  - `read` and non-blocking `write` are asynchronous; the host memory stays valid until the enqueued work
 *    completed (wait on a mark taken after the call, or finish).
 *  - A Shared buffer is host-visible through map/unmap (zero-copy on unified memory): the host owns the memory
 *    between map and unmap, the device must not touch it then. A Device buffer is reached by write/read only.
 *  - launch binds the arguments in order: buffer and value arguments by their position among the buffer and
 *    value arguments, local-memory arguments by their position among the local-memory arguments.
 *  - mark() names the completion of everything enqueued so far; wait(mark) blocks until then (and makes the
 *    reads enqueued before it visible), flush() starts the device on the enqueued work without blocking, finish()
 *    blocks until the stream is idle.
 */

/// Where a buffer lives: on the device only, or host-visible (map/unmap).
export enum class DeviceMemory : std::uint8_t
{
  Device = 0,
  Shared = 1
};

/// Floating-point contract of a program: Strict keeps every operation as written (needed by the double-float
/// arithmetic), Relaxed allows contractions and ignores signed zeros, Fast allows the approximate math too.
export enum class DeviceMath : std::uint8_t
{
  Strict = 0,
  Relaxed = 1,
  Fast = 2
};

export struct DeviceBuffer
{
  std::uint32_t id{0};
  explicit operator bool() const { return id != 0; }
};

export struct DeviceKernel
{
  std::uint32_t id{0};
  explicit operator bool() const { return id != 0; }
};

/// A point in the stream (see DeviceContext::mark).
export struct DeviceEvent
{
  std::uint64_t sequence{0};
};

/// One kernel argument of a launch.
export struct DeviceArg
{
  enum class Kind : std::uint8_t
  {
    Buffer,
    Value,
    Local
  };
  Kind kind{Kind::Buffer};
  DeviceBuffer buffer{};
  const void* data{nullptr};
  std::size_t bytes{0};

  static DeviceArg of(DeviceBuffer b) { return DeviceArg{Kind::Buffer, b, nullptr, 0}; }
  /// A VALUE_ARG of the kernel; `v` must be an lvalue that outlives the launch call (the bytes are copied by it).
  template <typename T>
    requires std::is_trivially_copyable_v<T>
  static DeviceArg value(const T& v)
  {
    return DeviceArg{Kind::Value, DeviceBuffer{}, &v, sizeof(T)};
  }
  /// A LOCAL pointer argument of `bytes` of local memory.
  static DeviceArg local(std::size_t bytes) { return DeviceArg{Kind::Local, DeviceBuffer{}, nullptr, bytes}; }
};

export class DeviceContext
{
 public:
  virtual ~DeviceContext() = default;

  virtual std::string deviceName() const = 0;
  /// Local (threadgroup) memory per work-group, bytes.
  virtual std::size_t localMemorySize() const = 0;

  virtual DeviceBuffer createBuffer(std::size_t bytes, DeviceMemory memory) = 0;
  virtual void releaseBuffer(DeviceBuffer& buffer) = 0;
  /// Host pointer to the first `bytes` of a Shared buffer (blocking: completes the enqueued device work on it).
  /// With `discard` the host does not need the device's contents of the range: either it rewrites everything it
  /// will read, or the device has not written the buffer since the host last unmapped it. A backend with separate
  /// host and device memories then skips the download; the mapped memory still holds what the host last wrote
  /// into it (it is not invalidated).
  virtual void* map(DeviceBuffer buffer, std::size_t bytes, bool forWriting, bool discard) = 0;
  virtual void unmap(DeviceBuffer buffer) = 0;
  /// Enqueues, in stream order, the download of the first `bytes` of a Shared buffer into its host copy, so that a
  /// map for reading after a wait for a later mark finds the contents there without synchronizing. A no-op with
  /// unified memory. A kernel launch, copy or write after it invalidates the download (the buffer may have changed).
  virtual void readback(DeviceBuffer buffer, std::size_t bytes) = 0;
  virtual void write(DeviceBuffer buffer, std::size_t offset, std::size_t bytes, const void* data, bool blocking) = 0;
  virtual void read(DeviceBuffer buffer, std::size_t offset, std::size_t bytes, void* data) = 0;
  /// Device-to-device copy, in stream order (asynchronous).
  virtual void copy(DeviceBuffer source, std::size_t sourceOffset, DeviceBuffer destination,
                    std::size_t destinationOffset, std::size_t bytes) = 0;

  /// A kernel of a program; programs are compiled once per `program` name (the dialect header prepended to
  /// `source`) and cached, so all kernels of a source share one compilation.
  virtual DeviceKernel compileKernel(std::string_view program, const char* source, DeviceMath math,
                                     std::string_view kernelName) = 0;
  virtual std::size_t maxGroupSize(DeviceKernel kernel) const = 0;
  virtual void launch(DeviceKernel kernel, std::span<const DeviceArg> arguments, std::size_t groups,
                      std::size_t groupSize) = 0;

  virtual DeviceEvent mark() = 0;
  virtual void wait(DeviceEvent event) = 0;
  virtual void flush() = 0;
  virtual void finish() = 0;
};

/// Owns a DeviceBuffer of a context (released on destruction / reset); movable.
export class DeviceBufferOwner
{
 public:
  DeviceBufferOwner() = default;
  explicit DeviceBufferOwner(DeviceContext* c) : context(c) {}
  ~DeviceBufferOwner() { reset(); }
  DeviceBufferOwner(const DeviceBufferOwner&) = delete;
  DeviceBufferOwner& operator=(const DeviceBufferOwner&) = delete;
  DeviceBufferOwner(DeviceBufferOwner&& other) noexcept
      : context(std::exchange(other.context, nullptr)), buffer(std::exchange(other.buffer, DeviceBuffer{}))
  {
  }
  DeviceBufferOwner& operator=(DeviceBufferOwner&& other) noexcept
  {
    if (this != &other)
    {
      reset();
      context = std::exchange(other.context, nullptr);
      buffer = std::exchange(other.buffer, DeviceBuffer{});
    }
    return *this;
  }
  /// Replaces the buffer by a new one of `bytes` (at least 1 byte is allocated).
  void allocate(DeviceContext& c, std::size_t bytes, DeviceMemory memory)
  {
    reset();
    context = &c;
    buffer = c.createBuffer(std::max<std::size_t>(bytes, 1), memory);
  }
  void reset()
  {
    if (context != nullptr && buffer) context->releaseBuffer(buffer);
    buffer = DeviceBuffer{};
  }
  DeviceBuffer get() const { return buffer; }
  explicit operator bool() const { return static_cast<bool>(buffer); }

 private:
  DeviceContext* context{nullptr};
  DeviceBuffer buffer{};
};
