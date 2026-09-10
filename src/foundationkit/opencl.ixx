module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module opencl;

import std;

export namespace OpenCL
{
extern std::optional<cl_context> clContext;
extern std::optional<cl_device_id> clDeviceId;
extern std::optional<cl_command_queue> clCommandQueue;

void initialize(void);
std::optional<cl_device_id> bestOpenCLDevice(cl_device_type device_type);
std::string printBestOpenCLDevice();
bool supportsImageFormatCapabilities(cl_context &trial_clContext, cl_device_id &trial_clDeviceId);

// Apple's OpenCL 1.2 treats a device allocation or copy of exactly 4 GiB as size 0.
// createBuffer pads that size; readBuffer/writeBuffer copy in 256 MiB chunks.
cl_mem createBuffer(cl_mem_flags flags, std::size_t bytes);
void writeBuffer(cl_mem mem, std::size_t bytes, const void *host);
void readBuffer(cl_mem mem, std::size_t bytes, void *host);
}  // namespace OpenCL
