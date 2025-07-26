module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <format>
#include <optional>
#include <ostream>
#include <print>
#include <sstream>
#include <string>
#include <vector>
#endif

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module opencl;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

std::optional<cl_context> OpenCL::clContext = std::nullopt;
std::optional<cl_device_id> OpenCL::clDeviceId = std::nullopt;
std::optional<cl_command_queue> OpenCL::clCommandQueue = std::nullopt;

void OpenCL::initialize()
{
  cl_int err;

  // Check first if there is a suitable GPU device for OpenCL
  std::optional<cl_device_id> bestGPUDevice = bestOpenCLDevice(CL_DEVICE_TYPE_GPU);

  if (bestGPUDevice)
  {
    std::optional<cl_device_id> device_id = bestGPUDevice;

    if (device_id.has_value())
    {
      cl_context context = clCreateContext(nullptr, 1, &device_id.value(), nullptr, nullptr, &err);
      if (err == CL_SUCCESS)
      {
        cl_command_queue command_queue = clCreateCommandQueue(context, device_id.value(), 0, &err);
        if (err == CL_SUCCESS)
        {
          clContext = context;
          clDeviceId = device_id;
          clCommandQueue = command_queue;
          return;
        }
        clReleaseContext(context);
      }
    }
  }

  // Check next if there is a suitable CPU device for OpenCL if no GPU is available
  std::optional<cl_device_id> bestCPUDevice = bestOpenCLDevice(CL_DEVICE_TYPE_CPU);

  if (bestCPUDevice)
  {
    std::optional<cl_device_id> device_id = bestCPUDevice;

    if (device_id.has_value())
    {
      cl_context context = clCreateContext(nullptr, 1, &device_id.value(), nullptr, nullptr, &err);
      if (err == CL_SUCCESS)
      {
        cl_command_queue command_queue = clCreateCommandQueue(context, device_id.value(), 0, &err);
        if (err == CL_SUCCESS)
        {
          clContext = context;
          clDeviceId = device_id;
          clCommandQueue = command_queue;
          return;
        }
        clReleaseContext(context);
      }
    }
  }
}

std::string OpenCL::printBestOpenCLDevice()
{
  std::ostringstream stream;

  if (clDeviceId.has_value())
  {
    cl_int err;

    char device_string[1024];
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_NAME, sizeof(device_string), device_string, NULL);
    if (err == CL_SUCCESS) stream << "OpenCL device: " << device_string << "\n";

    char vendor_string[1024];
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_VENDOR, sizeof(vendor_string), vendor_string, NULL);
    if (err == CL_SUCCESS) stream << "Vendor: " << vendor_string << "\n";

    char profile_string[1024];
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_PROFILE, sizeof(profile_string), profile_string, NULL);
    if (err == CL_SUCCESS) stream << "Profile: " << profile_string << "\n";

    char driver_string[1024];
    err = clGetDeviceInfo(clDeviceId.value(), CL_DRIVER_VERSION, sizeof(driver_string), driver_string, NULL);
    if (err == CL_SUCCESS) stream << "Driver version: " << driver_string << "\n";

    // CL_DEVICE_MAX_COMPUTE_UNITS
    cl_uint compute_units;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(compute_units), &compute_units, NULL);
    if (err == CL_SUCCESS) stream << "CL_DEVICE_MAX_COMPUTE_UNITS: " << compute_units << "\n";

    // CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
    size_t workitem_dims;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(workitem_dims), &workitem_dims,
                          NULL);
    if (err == CL_SUCCESS) stream << "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS:" << workitem_dims << "\n";

    // CL_DEVICE_MAX_WORK_ITEM_SIZES
    size_t workitem_size[3];
    err =
        clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(workitem_size), &workitem_size, NULL);
    if (err == CL_SUCCESS)
      stream << "CL_DEVICE_MAX_WORK_ITEM_SIZES: " << workitem_size[0] << "," << workitem_size[1] << ","
             << workitem_size[2] << "\n";

    // CL_DEVICE_MAX_WORK_GROUP_SIZE
    size_t workgroup_size;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(workgroup_size), &workgroup_size,
                          NULL);
    if (err == CL_SUCCESS) stream << "CL_DEVICE_MAX_WORK_GROUP_SIZE: " << workgroup_size << "\n";

    // CL_DEVICE_MAX_CLOCK_FREQUENCY
    cl_uint clock_frequency;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(clock_frequency), &clock_frequency,
                          NULL);
    if (err == CL_SUCCESS) stream << "CL_DEVICE_MAX_CLOCK_FREQUENCY:" << clock_frequency << " MHz\n";

    // CL_DEVICE_ADDRESS_BITS
    cl_uint addr_bits;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_ADDRESS_BITS, sizeof(addr_bits), &addr_bits, NULL);
    if (err == CL_SUCCESS) stream << "CL_DEVICE_ADDRESS_BITS: " << addr_bits << "\n";

    // CL_DEVICE_MAX_MEM_ALLOC_SIZE
    cl_ulong max_mem_alloc_size;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(max_mem_alloc_size),
                          &max_mem_alloc_size, NULL);
    if (err == CL_SUCCESS)
      stream << std::format("CL_DEVICE_MAX_MEM_ALLOC_SIZE: {} MByte\n",
                            static_cast<unsigned int>(max_mem_alloc_size / (1024 * 1024)));

    // CL_DEVICE_GLOBAL_MEM_SIZE
    cl_ulong mem_size;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
    if (err == CL_SUCCESS)
      stream << std::format("CL_DEVICE_GLOBAL_MEM_SIZE: {} MByte\n",
                            static_cast<unsigned int>(mem_size / (1024 * 1024)));

    // CL_DEVICE_ERROR_CORRECTION_SUPPORT
    cl_bool error_correction_support;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(error_correction_support),
                          &error_correction_support, NULL);
    if (err == CL_SUCCESS)
      stream << std::format("CL_DEVICE_ERROR_CORRECTION_SUPPORT: {}\n",
                            error_correction_support == CL_TRUE ? "yes" : "no");

    // CL_DEVICE_LOCAL_MEM_TYPE
    cl_device_local_mem_type local_mem_type;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_LOCAL_MEM_TYPE, sizeof(local_mem_type), &local_mem_type, NULL);
    if (err == CL_SUCCESS)
      stream << std::format("CL_DEVICE_LOCAL_MEM_TYPE: {}\n", local_mem_type == 1 ? "local" : "global");

    // CL_DEVICE_LOCAL_MEM_SIZE
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
    if (err == CL_SUCCESS)
      stream << std::format("CL_DEVICE_LOCAL_MEM_SIZE: {} KByte\n", static_cast<unsigned int>(mem_size / 1024));

    // CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(mem_size), &mem_size, NULL);
    if (err == CL_SUCCESS)
      stream << std::format("CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE: {} KByte\n",
                            static_cast<unsigned int>(mem_size / 1024));

    // CL_DEVICE_IMAGE_SUPPORT
    cl_bool image_support;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_IMAGE_SUPPORT, sizeof(image_support), &image_support, NULL);
    if (err == CL_SUCCESS) stream << std::format("CL_DEVICE_IMAGE_SUPPORT: {}\n", image_support);

    // CL_DEVICE_MAX_READ_IMAGE_ARGS
    cl_uint max_read_image_args;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_MAX_READ_IMAGE_ARGS, sizeof(max_read_image_args),
                          &max_read_image_args, NULL);
    if (err == CL_SUCCESS) stream << std::format("CL_DEVICE_MAX_READ_IMAGE_ARGS: {}\n", max_read_image_args);

    // CL_DEVICE_MAX_WRITE_IMAGE_ARGS
    cl_uint max_write_image_args;
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_MAX_WRITE_IMAGE_ARGS, sizeof(max_write_image_args),
                          &max_write_image_args, NULL);
    if (err == CL_SUCCESS) stream << std::format("CL_DEVICE_MAX_WRITE_IMAGE_ARGS: {}\n", max_write_image_args);

    // CL_DEVICE_IMAGE2D_MAX_WIDTH, CL_DEVICE_IMAGE2D_MAX_HEIGHT, CL_DEVICE_IMAGE3D_MAX_WIDTH,
    // CL_DEVICE_IMAGE3D_MAX_HEIGHT, CL_DEVICE_IMAGE3D_MAX_DEPTH
    size_t szMaxDims[5];
    if (err == CL_SUCCESS) stream << "CL_DEVICE_IMAGE <dim>\n";
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_IMAGE2D_MAX_WIDTH, sizeof(size_t), &szMaxDims[0], NULL);
    if (err == CL_SUCCESS) stream << std::format("2D_MAX_WIDTH {}\n", szMaxDims[0]);
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_IMAGE2D_MAX_HEIGHT, sizeof(size_t), &szMaxDims[1], NULL);
    if (err == CL_SUCCESS) stream << std::format("2D_MAX_HEIGHT {}\n", szMaxDims[1]);
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_IMAGE3D_MAX_WIDTH, sizeof(size_t), &szMaxDims[2], NULL);
    if (err == CL_SUCCESS) stream << std::format("3D_MAX_WIDTH {}\n", szMaxDims[2]);
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_IMAGE3D_MAX_HEIGHT, sizeof(size_t), &szMaxDims[3], NULL);
    if (err == CL_SUCCESS) stream << std::format("3D_MAX_HEIGHT {}\n", szMaxDims[3]);
    err = clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_IMAGE3D_MAX_DEPTH, sizeof(size_t), &szMaxDims[4], NULL);
    if (err == CL_SUCCESS) stream << std::format("3D_MAX_DEPTH {}\n", szMaxDims[4]);

    // CL_DEVICE_PREFERRED_VECTOR_WIDTH_<type>
    cl_uint vec_width[6];
    clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, sizeof(cl_uint), &vec_width[0], NULL);
    clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, sizeof(cl_uint), &vec_width[1], NULL);
    clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, sizeof(cl_uint), &vec_width[2], NULL);
    clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, sizeof(cl_uint), &vec_width[3], NULL);
    clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, sizeof(cl_uint), &vec_width[4], NULL);
    clGetDeviceInfo(clDeviceId.value(), CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(cl_uint), &vec_width[5], NULL);
    if (err == CL_SUCCESS)
    {
      stream << "CL_DEVICE_PREFERRED_VECTOR_WIDTH\n";
      stream << std::format("CHAR {}, SHORT {}, INT {}, FLOAT {}, DOUBLE {}\n", vec_width[0], vec_width[1],
                            vec_width[2], vec_width[3], vec_width[4]);
    }

    return stream.str();
  }

  return std::string("No OpenCL devise available\n");
}

std::optional<cl_device_id> OpenCL::bestOpenCLDevice(cl_device_type device_type)
{
  cl_int err;

  // get the number of platforms
  cl_uint platformCount;
  err = clGetPlatformIDs(0, nullptr, &platformCount);

  if (platformCount <= 0)
  {
    return std::nullopt;
  }

  // get platform ids
  std::vector<cl_platform_id> platforms;
  platforms.resize(platformCount);
  err = clGetPlatformIDs(static_cast<cl_uint>(platforms.size()), platforms.data(), nullptr);
  if (err != CL_SUCCESS)
  {
    return std::nullopt;
  }

  std::optional<std::pair<cl_device_id, cl_uint>> bestDevice = std::nullopt;

  // loop over all platforms and devices of type 'device_type'
  for (cl_uint i = 0; i < platformCount; ++i)
  {
    // get the number of devices of type 'device_type' in the platform
    cl_uint deviceCount;
    err = clGetDeviceIDs(platforms[i], device_type, 0, nullptr, &deviceCount);
    if (err != CL_SUCCESS)
    {
      continue;
    }

    // get the IDs of GPU devices
    std::vector<cl_device_id> devices;
    devices.resize(deviceCount);
    err = clGetDeviceIDs(platforms[i], device_type, static_cast<cl_uint>(devices.size()), devices.data(), nullptr);
    if (err != CL_SUCCESS)
    {
      continue;
    }

    // loop over all devices
    for (cl_uint j = 0; j < deviceCount; j++)
    {
      cl_uint maxComputeUnits;
      err = clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(maxComputeUnits), &maxComputeUnits, NULL);
      if (err != CL_SUCCESS)
      {
        continue;
      }

      cl_context context = clCreateContext(nullptr, 1, &devices[j], nullptr, nullptr, &err);

      if (err != CL_SUCCESS)
      {
        continue;
      }

      if (supportsImageFormatCapabilities(context, devices[j]))
      {
        cl_command_queue command_queue = clCreateCommandQueue(context, devices[j], 0, &err);
        if (err != CL_SUCCESS)
        {
          clReleaseContext(context);
          continue;
        }
        clReleaseCommandQueue(command_queue);

        if (bestDevice)
        {
          // use the one with the highest compute units
          if (maxComputeUnits > std::get<1>(*bestDevice))
          {
            bestDevice = std::make_pair(devices[j], maxComputeUnits);
          }
        }
        else
        {
          bestDevice = std::make_pair(devices[j], maxComputeUnits);
        }
      }

      clReleaseContext(context);
    }
  }

  if (bestDevice)
  {
    // we have a suitable device that can be used
    return std::get<0>(*bestDevice);
  }

  return std::nullopt;
}

bool OpenCL::supportsImageFormatCapabilities(cl_context &trial_clContext, cl_device_id &trialclDeviceId)
{
  cl_int err;

  // check image support
  cl_bool image_support = false;
  clGetDeviceInfo(trialclDeviceId, CL_DEVICE_IMAGE_SUPPORT, sizeof(image_support), &image_support, nullptr);

  if (!image_support)
  {
    return false;
  }

  // check the needed image formats
  cl_image_format imageFormat_RGBA_INT8{CL_RGBA, CL_UNSIGNED_INT8};
  cl_image_desc imageDescriptor_RGBA_INT8{CL_MEM_OBJECT_IMAGE3D, 256, 256, 256, 0, 0, 0, 0, 0, nullptr};
  cl_mem image_RGBA_INT8 = clCreateImage(trial_clContext, CL_MEM_READ_WRITE, &imageFormat_RGBA_INT8,
                                         &imageDescriptor_RGBA_INT8, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    return false;
  }
  clReleaseMemObject(image_RGBA_INT8);

  cl_image_format imageFormat_R_INT8{CL_R, CL_UNSIGNED_INT8};
  cl_image_desc imageDescriptor_R_INT8{CL_MEM_OBJECT_IMAGE3D, 256, 256, 256, 0, 0, 0, 0, 0, nullptr};
  cl_mem image_R_INT8 =
      clCreateImage(trial_clContext, CL_MEM_READ_WRITE, &imageFormat_R_INT8, &imageDescriptor_R_INT8, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    return false;
  }
  clReleaseMemObject(image_R_INT8);

  cl_image_format imageFormat_R_INT16{CL_R, CL_UNSIGNED_INT16};
  cl_image_desc imageDescriptor_R_INT16{CL_MEM_OBJECT_IMAGE3D, 256, 256, 256, 0, 0, 0, 0, 0, nullptr};
  cl_mem image_R_INT16 =
      clCreateImage(trial_clContext, CL_MEM_READ_WRITE, &imageFormat_R_INT16, &imageDescriptor_R_INT16, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    return false;
  }
  clReleaseMemObject(image_R_INT16);

  cl_image_format imageFormat_R_INT32{CL_R, CL_UNSIGNED_INT16};
  cl_image_desc imageDescriptor_R_INT32{CL_MEM_OBJECT_IMAGE3D, 128, 128, 128, 0, 0, 0, 0, 0, nullptr};
  cl_mem image_R_INT32 =
      clCreateImage(trial_clContext, CL_MEM_READ_WRITE, &imageFormat_R_INT32, &imageDescriptor_R_INT32, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    return false;
  }
  clReleaseMemObject(image_R_INT32);

  cl_image_format imageFormat_R_FLOAT{CL_R, CL_FLOAT};
  cl_image_desc imageDescriptor_R_FLOAT{CL_MEM_OBJECT_IMAGE3D, 128, 128, 128, 0, 0, 0, 0, 0, nullptr};
  cl_mem image_R_FLOAT =
      clCreateImage(trial_clContext, CL_MEM_READ_WRITE, &imageFormat_R_FLOAT, &imageDescriptor_R_FLOAT, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    return false;
  }
  clReleaseMemObject(image_R_FLOAT);

  return true;
}
