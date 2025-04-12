module;

#include <cstddef>
#include <optional>
#include <string>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

export module opencl;

export namespace OpenCL
{
  std::optional<cl_context> clContext{std::nullopt};
  std::optional<cl_device_id> clDeviceId{std::nullopt};
  std::optional<cl_command_queue> clCommandQueue{std::nullopt};

  void initialize(void);
  std::optional<cl_device_id> bestOpenCLDevice(cl_device_type device_type);
  std::string printBestOpenCLDevice();
  bool supportsImageFormatCapabilities(cl_context &trial_clContext, cl_device_id &trial_clDeviceId);
}

