// The Objective-C++ side of the Metal bridge (bridge.h). Compiled with ARC; the opaque handles own strong
// references to the Metal objects.

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <string.h>

#include "bridge.h"

struct RaspaMetalDevice
{
  id<MTLDevice> device;
  id<MTLCommandQueue> queue;
};

struct RaspaMetalBuffer
{
  id<MTLBuffer> buffer;
};

struct RaspaMetalLibrary
{
  id<MTLLibrary> library;
};

struct RaspaMetalPipeline
{
  id<MTLComputePipelineState> pipeline;
};

struct RaspaMetalCommands
{
  id<MTLCommandBuffer> commandBuffer;
  id<MTLComputeCommandEncoder> compute;
  id<MTLBlitCommandEncoder> blit;
  bool committed;
};

namespace
{
void writeError(char* error, size_t capacity, NSString* text)
{
  if (error == NULL || capacity == 0) return;
  const char* utf8 = text == nil ? "" : [text UTF8String];
  strncpy(error, utf8, capacity - 1);
  error[capacity - 1] = '\0';
}

void endEncoders(RaspaMetalCommands* commands)
{
  if (commands->compute != nil)
  {
    [commands->compute endEncoding];
    commands->compute = nil;
  }
  if (commands->blit != nil)
  {
    [commands->blit endEncoding];
    commands->blit = nil;
  }
}
}  // namespace

extern "C"
{
  int raspaMetalAvailable(void)
  {
    @autoreleasepool
    {
      id<MTLDevice> device = MTLCreateSystemDefaultDevice();
      return device != nil ? 1 : 0;
    }
  }

  size_t raspaMetalDeviceName(char* name, size_t capacity)
  {
    @autoreleasepool
    {
      if (name == NULL || capacity == 0) return 0;
      id<MTLDevice> device = MTLCreateSystemDefaultDevice();
      if (device == nil)
      {
        name[0] = '\0';
        return 0;
      }
      writeError(name, capacity, [device name]);
      return strlen(name);
    }
  }

  RaspaMetalDevice* raspaMetalCreateDevice(char* error, size_t errorCapacity)
  {
    @autoreleasepool
    {
      id<MTLDevice> device = MTLCreateSystemDefaultDevice();
      if (device == nil)
      {
        writeError(error, errorCapacity, @"no Metal device is available");
        return NULL;
      }
      id<MTLCommandQueue> queue = [device newCommandQueue];
      if (queue == nil)
      {
        writeError(error, errorCapacity, @"the Metal command queue could not be created");
        return NULL;
      }
      RaspaMetalDevice* handle = new RaspaMetalDevice{};
      handle->device = device;
      handle->queue = queue;
      return handle;
    }
  }

  void raspaMetalDestroyDevice(RaspaMetalDevice* device)
  {
    if (device == NULL) return;
    @autoreleasepool
    {
      device->queue = nil;
      device->device = nil;
    }
    delete device;
  }

  size_t raspaMetalLocalMemorySize(const RaspaMetalDevice* device)
  {
    if (device == NULL) return 0;
    return static_cast<size_t>([device->device maxThreadgroupMemoryLength]);
  }

  RaspaMetalBuffer* raspaMetalCreateBuffer(RaspaMetalDevice* device, size_t bytes)
  {
    @autoreleasepool
    {
      id<MTLBuffer> buffer = [device->device newBufferWithLength:(bytes == 0 ? 16 : bytes)
                                                         options:MTLResourceStorageModeShared];
      if (buffer == nil) return NULL;
      RaspaMetalBuffer* handle = new RaspaMetalBuffer{};
      handle->buffer = buffer;
      return handle;
    }
  }

  void raspaMetalDestroyBuffer(RaspaMetalBuffer* buffer)
  {
    if (buffer == NULL) return;
    @autoreleasepool
    {
      buffer->buffer = nil;
    }
    delete buffer;
  }

  void* raspaMetalBufferContents(RaspaMetalBuffer* buffer) { return [buffer->buffer contents]; }

  RaspaMetalLibrary* raspaMetalCompile(RaspaMetalDevice* device, const char* source, int mathMode, char* error,
                                       size_t errorCapacity)
  {
    @autoreleasepool
    {
      MTLCompileOptions* options = [[MTLCompileOptions alloc] init];
      if (@available(macOS 15.0, *))
      {
        options.mathMode = mathMode == 0 ? MTLMathModeSafe : mathMode == 1 ? MTLMathModeRelaxed : MTLMathModeFast;
      }
      else
      {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
        options.fastMathEnabled = mathMode == 2;
#pragma clang diagnostic pop
      }
      NSError* nsError = nil;
      id<MTLLibrary> library = [device->device newLibraryWithSource:[NSString stringWithUTF8String:source]
                                                            options:options
                                                              error:&nsError];
      if (library == nil)
      {
        writeError(error, errorCapacity, nsError == nil ? @"unknown Metal compile error" : [nsError localizedDescription]);
        return NULL;
      }
      RaspaMetalLibrary* handle = new RaspaMetalLibrary{};
      handle->library = library;
      return handle;
    }
  }

  void raspaMetalDestroyLibrary(RaspaMetalLibrary* library)
  {
    if (library == NULL) return;
    @autoreleasepool
    {
      library->library = nil;
    }
    delete library;
  }

  RaspaMetalPipeline* raspaMetalCreatePipeline(RaspaMetalDevice* device, RaspaMetalLibrary* library,
                                               const char* kernelName, char* error, size_t errorCapacity)
  {
    @autoreleasepool
    {
      id<MTLFunction> function = [library->library newFunctionWithName:[NSString stringWithUTF8String:kernelName]];
      if (function == nil)
      {
        writeError(error, errorCapacity,
                   [NSString stringWithFormat:@"kernel '%s' not found in the Metal library", kernelName]);
        return NULL;
      }
      NSError* nsError = nil;
      id<MTLComputePipelineState> pipeline = [device->device newComputePipelineStateWithFunction:function
                                                                                           error:&nsError];
      if (pipeline == nil)
      {
        writeError(error, errorCapacity,
                   nsError == nil ? @"the compute pipeline could not be created" : [nsError localizedDescription]);
        return NULL;
      }
      RaspaMetalPipeline* handle = new RaspaMetalPipeline{};
      handle->pipeline = pipeline;
      return handle;
    }
  }

  void raspaMetalDestroyPipeline(RaspaMetalPipeline* pipeline)
  {
    if (pipeline == NULL) return;
    @autoreleasepool
    {
      pipeline->pipeline = nil;
    }
    delete pipeline;
  }

  size_t raspaMetalMaxThreadsPerThreadgroup(const RaspaMetalPipeline* pipeline)
  {
    return static_cast<size_t>([pipeline->pipeline maxTotalThreadsPerThreadgroup]);
  }

  RaspaMetalCommands* raspaMetalBeginCommands(RaspaMetalDevice* device)
  {
    @autoreleasepool
    {
      id<MTLCommandBuffer> commandBuffer = [device->queue commandBufferWithUnretainedReferences];
      if (commandBuffer == nil) return NULL;
      RaspaMetalCommands* handle = new RaspaMetalCommands{};
      handle->commandBuffer = commandBuffer;
      handle->compute = nil;
      handle->blit = nil;
      handle->committed = false;
      return handle;
    }
  }

  void raspaMetalDispatch(RaspaMetalCommands* commands, RaspaMetalPipeline* pipeline, RaspaMetalBuffer* const* buffers,
                          const uint32_t* bufferIndices, size_t bufferCount, const void* const* values,
                          const size_t* valueBytes, const uint32_t* valueIndices, size_t valueCount,
                          const size_t* localBytes, size_t localCount, size_t groups, size_t groupSize)
  {
    @autoreleasepool
    {
      if (commands->blit != nil)
      {
        [commands->blit endEncoding];
        commands->blit = nil;
      }
      if (commands->compute == nil)
      {
        commands->compute = [commands->commandBuffer computeCommandEncoder];
      }
      id<MTLComputeCommandEncoder> encoder = commands->compute;
      [encoder setComputePipelineState:pipeline->pipeline];
      for (size_t k = 0; k < bufferCount; ++k)
      {
        [encoder setBuffer:buffers[k]->buffer offset:0 atIndex:bufferIndices[k]];
      }
      for (size_t k = 0; k < valueCount; ++k)
      {
        [encoder setBytes:values[k] length:valueBytes[k] atIndex:valueIndices[k]];
      }
      for (size_t k = 0; k < localCount; ++k)
      {
        // threadgroup memory lengths are multiples of 16 bytes
        const size_t rounded = ((localBytes[k] + 15) / 16) * 16;
        [encoder setThreadgroupMemoryLength:(rounded == 0 ? 16 : rounded) atIndex:k];
      }
      [encoder dispatchThreadgroups:MTLSizeMake(groups == 0 ? 1 : groups, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(groupSize, 1, 1)];
    }
  }

  void raspaMetalCopy(RaspaMetalCommands* commands, RaspaMetalBuffer* source, size_t sourceOffset,
                      RaspaMetalBuffer* destination, size_t destinationOffset, size_t bytes)
  {
    @autoreleasepool
    {
      if (bytes == 0) return;
      if (commands->compute != nil)
      {
        [commands->compute endEncoding];
        commands->compute = nil;
      }
      if (commands->blit == nil)
      {
        commands->blit = [commands->commandBuffer blitCommandEncoder];
      }
      [commands->blit copyFromBuffer:source->buffer
                        sourceOffset:sourceOffset
                            toBuffer:destination->buffer
                   destinationOffset:destinationOffset
                                size:bytes];
    }
  }

  void raspaMetalCommit(RaspaMetalCommands* commands)
  {
    @autoreleasepool
    {
      if (commands->committed) return;
      endEncoders(commands);
      [commands->commandBuffer commit];
      commands->committed = true;
    }
  }

  int raspaMetalWait(RaspaMetalCommands* commands, char* error, size_t errorCapacity)
  {
    @autoreleasepool
    {
      if (!commands->committed) raspaMetalCommit(commands);
      [commands->commandBuffer waitUntilCompleted];
      if (commands->commandBuffer.status == MTLCommandBufferStatusError)
      {
        NSError* nsError = commands->commandBuffer.error;
        writeError(error, errorCapacity,
                   nsError == nil ? @"the Metal command buffer failed" : [nsError localizedDescription]);
        return 1;
      }
      return 0;
    }
  }

  void raspaMetalDestroyCommands(RaspaMetalCommands* commands)
  {
    if (commands == NULL) return;
    @autoreleasepool
    {
      endEncoders(commands);
      commands->commandBuffer = nil;
    }
    delete commands;
  }
}
