module;
     
#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>  
#include <chrono> 
#include <complex>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <exception>
#include <format>
#include <fstream>
#include <istream>
#include <map>                                      
#include <ostream>                                  
#include <print>                                    
#include <source_location>
#include <sstream>
#include <tuple>          
#include <utility>
#include <vector>
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
  #include <OpenCL/cl.h>
#elif _WIN32
  #include <CL/cl.h>
#else
  #include <CL/opencl.h>
#endif
#endif
  
module energy_void_fraction;

import opencl;
import float4;
import double4;

EnergyVoidFraction::EnergyVoidFraction()
{
  if(OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;
    const char* shaderSourceCode = EnergyVoidFraction::energyVoidFractionKernel;
    program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &shaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error("clCreateProgramWithSource failed");}

    err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error("clBuildProgram failed");}

    kernel = clCreateKernel(program, "ComputeVoidFraction", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error("clCreateKernel failed");}

    err = clGetKernelWorkGroupInfo(kernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &workGroupSize, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error("clGetKernelWorkGroupInfo failed");}
  }
};

EnergyVoidFraction::~EnergyVoidFraction()
{
  if(OpenCL::clContext.has_value())
  {
    clReleaseKernel(kernel);
    clReleaseProgram(program);
  }
}


double EnergyVoidFraction::computeVoidFraction(std::vector<cl_float> *voxels)
{
  if(OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value() && OpenCL::clCommandQueue.has_value())
  {
    cl_int err;

    // make sure the the global work size is an multiple of the work group size
    size_t temp = voxels->size();
    size_t numberOfGridPoints = (temp  + workGroupSize-1) & ~(workGroupSize-1);
    size_t global_work_size = numberOfGridPoints;


    // Transfer dataset to device
    cl_mem rawData = clCreateBuffer(OpenCL::clContext.value(),  CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,  sizeof(cl_float) * voxels->size(), voxels->data(), &err);
    if (err != CL_SUCCESS) {throw std::runtime_error("SKComputeVoidFraction: error in clCreateBuffer");}

    size_t nWorkGroups = numberOfGridPoints/workGroupSize;

    // Allocate cumulative error array
    float* sumReduction = new float[nWorkGroups];
    cl_mem reductionBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, nWorkGroups*sizeof(float), nullptr, &err);

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), reductionBuffer, CL_TRUE, 0, nWorkGroups*sizeof(float), sumReduction, 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error("SKComputeVoidFraction: error in clEnqueueWriteBuffer");}

    // Set the arguments of the kernel
    clSetKernelArg(kernel, 0, sizeof(cl_mem), &rawData);
    clSetKernelArg(kernel, 1, sizeof(cl_mem), &reductionBuffer);
    clSetKernelArg(kernel, 2, workGroupSize*sizeof(cl_float),nullptr);

    // Execute kernel code
    err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 1, nullptr, &global_work_size, &workGroupSize, 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error("SKComputeVoidFraction: error in clEnqueueNDRangeKernel");}

    // Read the buffer back to the array
    err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), reductionBuffer, CL_TRUE, 0, nWorkGroups*sizeof(float), sumReduction, 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error("SKComputeVoidFraction: error in clEnqueueReadBuffer");}

    clFinish(OpenCL::clCommandQueue.value());

    // Final summation with CPU
    double fraction = 0.0;
    for (size_t i=0; i<nWorkGroups; i++)
      fraction += double(sumReduction[i]);

    clReleaseMemObject(rawData);
    clReleaseMemObject(reductionBuffer);

    return fraction/(128.0*128.0*128.0);
  }

  return 0.0;
}

const char* EnergyVoidFraction::energyVoidFractionKernel = R"foo(
__kernel void ComputeVoidFraction (__global const float *input,
                           __global float *partialSums,
                           __local float *localSums)
{
   size_t local_id = get_local_id(0);
   size_t global_id = get_global_id(0);
   size_t group_size = get_local_size(0);
   size_t group_id = get_group_id(0);

   // Copy from global memory to local memory
   localSums[local_id] = exp(-(1.0f/298.0f)*input[global_id]);

   // Loop for computing localSums
   for (uint stride = group_size/2; stride>0; stride/=2)
   {
      // Waiting for each 2x2 addition into given workgroup
      barrier(CLK_LOCAL_MEM_FENCE);

      // Divide WorkGroup into 2 parts and add elements 2 by 2
      // between local_id and local_id + stride
      if (local_id < stride)
        localSums[local_id] += localSums[local_id + stride];
   }

   // Write result into partialSums[nWorkGroups]
   if (local_id == 0)
     partialSums[group_id] = localSums[0];
}
)foo";

