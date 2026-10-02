#pragma once

// The C interface of the Metal bridge (bridge.mm, Objective-C++) used by the C++ MetalContext
// (spatial_decomposition_metal_context). Plain C only: opaque handles, C strings and byte counts, so that nothing
// of the C++ runtime crosses between the Objective-C++ unit and the C++ modules. All handles belong to one
// RaspaMetalDevice and must be destroyed before it. Errors are reported through a caller-provided character
// buffer; functions returning a pointer return NULL on failure.
//
// Stream model: one MTLCommandQueue per device; the caller records commands into a command buffer (one at a
// time), commits it, and later waits for its completion. Command buffers of one queue execute in commit order.

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct RaspaMetalDevice RaspaMetalDevice;
  typedef struct RaspaMetalBuffer RaspaMetalBuffer;
  typedef struct RaspaMetalLibrary RaspaMetalLibrary;
  typedef struct RaspaMetalPipeline RaspaMetalPipeline;
  typedef struct RaspaMetalCommands RaspaMetalCommands;

  /// 1 when a Metal device exists on this machine.
  int raspaMetalAvailable(void);
  /// Name of the default device (empty string when none); returns the length written.
  size_t raspaMetalDeviceName(char* name, size_t capacity);

  RaspaMetalDevice* raspaMetalCreateDevice(char* error, size_t errorCapacity);
  void raspaMetalDestroyDevice(RaspaMetalDevice* device);
  /// Threadgroup memory per threadgroup, bytes.
  size_t raspaMetalLocalMemorySize(const RaspaMetalDevice* device);

  /// A buffer in shared storage (unified memory: host-visible through raspaMetalBufferContents).
  RaspaMetalBuffer* raspaMetalCreateBuffer(RaspaMetalDevice* device, size_t bytes);
  void raspaMetalDestroyBuffer(RaspaMetalBuffer* buffer);
  void* raspaMetalBufferContents(RaspaMetalBuffer* buffer);

  /// Compiles MSL source; mathMode 0 safe (no reassociation, no contraction), 1 relaxed, 2 fast.
  RaspaMetalLibrary* raspaMetalCompile(RaspaMetalDevice* device, const char* source, int mathMode, char* error,
                                       size_t errorCapacity);
  void raspaMetalDestroyLibrary(RaspaMetalLibrary* library);
  RaspaMetalPipeline* raspaMetalCreatePipeline(RaspaMetalDevice* device, RaspaMetalLibrary* library,
                                               const char* kernelName, char* error, size_t errorCapacity);
  void raspaMetalDestroyPipeline(RaspaMetalPipeline* pipeline);
  size_t raspaMetalMaxThreadsPerThreadgroup(const RaspaMetalPipeline* pipeline);

  /// A new command buffer to record into.
  RaspaMetalCommands* raspaMetalBeginCommands(RaspaMetalDevice* device);
  /// A compute dispatch: `bufferCount` buffers at buffer indices 0.., `valueCount` by-value arguments
  /// interleaved by `bufferIndex` (each value occupies one buffer index), `localCount` threadgroup memory sizes at
  /// threadgroup indices 0...
  void raspaMetalDispatch(RaspaMetalCommands* commands, RaspaMetalPipeline* pipeline, RaspaMetalBuffer* const* buffers,
                          const uint32_t* bufferIndices, size_t bufferCount, const void* const* values,
                          const size_t* valueBytes, const uint32_t* valueIndices, size_t valueCount,
                          const size_t* localBytes, size_t localCount, size_t groups, size_t groupSize);
  /// A copy between buffers (blit).
  void raspaMetalCopy(RaspaMetalCommands* commands, RaspaMetalBuffer* source, size_t sourceOffset,
                      RaspaMetalBuffer* destination, size_t destinationOffset, size_t bytes);
  /// Ends the encoders and commits; the handle stays valid for raspaMetalWait / raspaMetalDestroyCommands.
  void raspaMetalCommit(RaspaMetalCommands* commands);
  /// Blocks until the committed command buffer completed; returns 0 on success, else writes the error.
  int raspaMetalWait(RaspaMetalCommands* commands, char* error, size_t errorCapacity);
  void raspaMetalDestroyCommands(RaspaMetalCommands* commands);

#ifdef __cplusplus
}
#endif
