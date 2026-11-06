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

// What a routine made on the device, let go of when the routine is left, by whichever way it is left.
//
// A routine that builds a program, makes a kernel of it and a few buffers, and releases them on its last lines
// releases them only if it gets to its last lines. Every call in between can fail, and a failure is thrown:
// whatever had been made by then stayed on the device, a program and its kernel and each buffer, for as long
// as the context lived. Held here instead, they are released when this goes out of scope, thrown past or not.
// A handle that was never made is left as nothing and nothing is done about it.
struct Handles
{
  cl_program program{nullptr};
  cl_kernel kernel{nullptr};
  std::vector<cl_mem> buffers;

  Handles() = default;
  Handles(const Handles &) = delete;
  Handles &operator=(const Handles &) = delete;
  ~Handles();

  // Takes a buffer to be released with the rest and hands it back, so that it can be made and kept in a line.
  cl_mem keep(cl_mem buffer);
};
}  // namespace OpenCL
