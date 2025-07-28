module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <print>
#include <string>
#include <tuple>
#include <vector>
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

#ifndef USE_LEGACY_HEADERS
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

module tessellation;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import opencl;
import float4;
import double2;
import double3;
import double4;
import double3x3;
import randomnumbers;
import atom;
import framework;
import forcefield;
import component;
import system;
import units;
#if !(defined(__has_include) && __has_include(<mdspan>))
import mdspan;
#endif

Tessellation::Tessellation(int3 grid_size) : grid_size(grid_size)
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* tessellationShaderSourceCode = Tessellation::tessellationKernelSource;
    tessellationProgram =
        clCreateProgramWithSource(OpenCL::clContext.value(), 1, &tessellationShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
    }

    err = clBuildProgram(tessellationProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      std::size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(tessellationProgram, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer),
                            buffer, &len);
      std::string message = std::format("Tessellation: OpenCL Failed to build program at {} (line {} error: {})\n",
                                        __FILE__, __LINE__, std::string(buffer));
      throw std::runtime_error(message);
    }

    tessellationKernel = clCreateKernel(tessellationProgram, "Tessellation", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateKernel failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clGetKernelWorkGroupInfo(tessellationKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(std::size_t), &tessellationWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));
    }
  }
}

Tessellation::~Tessellation()
{
  if (OpenCL::clContext.has_value())
  {
    clReleaseKernel(tessellationKernel);
    clReleaseProgram(tessellationProgram);
  }
}

void Tessellation::run(const ForceField& forceField, const Framework& framework)
{
  double2 probeParameter = double2(10.9 * Units::KelvinToEnergy, 2.64);
  double cutoff = forceField.cutOffFrameworkVDW;
  double3x3 unitCell = framework.simulationBox.cell;
  int3 numberOfReplicas = framework.simulationBox.smallestNumberOfUnitCellsForMinimumImagesConvention(cutoff);
  std::vector<double3> positions = framework.fractionalAtomPositionsUnitCell();
  std::vector<double2> potentialParameters = framework.atomUnitCellLennardJonesPotentialParameters(forceField);
  std::chrono::system_clock::time_point time_begin, time_end;
  cl_int err;

  if (!OpenCL::clContext.has_value()) return;

  time_begin = std::chrono::system_clock::now();

  std::vector<cl_float4> pos(positions.size());
  std::vector<cl_float> sigma(positions.size());

  cl_mem image;
  cl_image_format cl_image_format{CL_R, CL_UNSIGNED_INT32};
  cl_image_desc imageDescriptor = cl_image_desc{CL_MEM_OBJECT_IMAGE3D,
                                                static_cast<std::size_t>(grid_size.x),
                                                static_cast<std::size_t>(grid_size.y),
                                                static_cast<std::size_t>(grid_size.z),
                                                0,
                                                0,
                                                0,
                                                0,
                                                0,
                                                nullptr};
  image =
      clCreateImage(OpenCL::clContext.value(), CL_MEM_READ_WRITE, &cl_image_format, &imageDescriptor, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clCreateImage failed at {} line {}\n", __FILE__, __LINE__));
  }

  std::vector<int32_t> output_data(static_cast<std::size_t>(grid_size.x * grid_size.y * grid_size.x));

  if (positions.size() > 0)
  {
    double3 correction =
        double3(1.0 / double(numberOfReplicas.x), 1.0 / double(numberOfReplicas.y), 1.0 / double(numberOfReplicas.z));
    for (std::size_t i = 0; i < positions.size(); ++i)
    {
      double3 position = correction * positions[i];
      double2 currentPotentialParameters = potentialParameters[i];

      // fill in the fractional position
      pos[i] = {{cl_float(position.x), cl_float(position.y), cl_float(position.z), 0.0f}};

      // mixing rule for the atom and the probe
      sigma[i] = cl_float(currentPotentialParameters.y);
    }

    cl_mem inputPos =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(float4) * pos.size(), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputPos, CL_TRUE, 0, sizeof(float4) * pos.size(),
                               pos.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
    }

    cl_mem inputSigma =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * sigma.size(), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputSigma, CL_TRUE, 0, sizeof(cl_float) * sigma.size(),
                               sigma.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
    }

    double3x3 replicaCell =
        double3x3(double(numberOfReplicas.x) * unitCell[0], double(numberOfReplicas.y) * unitCell[1],
                  double(numberOfReplicas.z) * unitCell[2]);

    cl_float4 clCella = {
        {cl_float(replicaCell[0][0]), cl_float(replicaCell[1][0]), cl_float(replicaCell[2][0]), cl_float(0.0)}};
    cl_float4 clCellb = {
        {cl_float(replicaCell[0][1]), cl_float(replicaCell[1][1]), cl_float(replicaCell[2][1]), cl_float(0.0)}};
    cl_float4 clCellc = {
        {cl_float(replicaCell[0][2]), cl_float(replicaCell[1][2]), cl_float(replicaCell[2][2]), cl_float(0.0)}};

    cl_int numberOfAtoms = static_cast<cl_int>(positions.size());
    err = clSetKernelArg(tessellationKernel, 0, sizeof(cl_mem), &inputPos);
    err |= clSetKernelArg(tessellationKernel, 1, sizeof(cl_mem), &inputSigma);
    err |= clSetKernelArg(tessellationKernel, 2, sizeof(cl_mem), &image);
    err |= clSetKernelArg(tessellationKernel, 3, sizeof(cl_float4), &clCella);
    err |= clSetKernelArg(tessellationKernel, 4, sizeof(cl_float4), &clCellb);
    err |= clSetKernelArg(tessellationKernel, 5, sizeof(cl_float4), &clCellc);
    err |= clSetKernelArg(tessellationKernel, 6, sizeof(cl_int), &numberOfAtoms);
    err |= clSetKernelArg(tessellationKernel, 7, sizeof(cl_int3), &numberOfReplicas);
    err |= clSetKernelArg(tessellationKernel, 8, sizeof(cl_int3), &grid_size);

    std::size_t global_work_size[3] = {static_cast<std::size_t>(grid_size.x), static_cast<std::size_t>(grid_size.y),
                                       static_cast<std::size_t>(grid_size.z)};
    err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), tessellationKernel, 3, nullptr, global_work_size,
                                 nullptr, 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("OpenCL clEnqueueNDRangeKernel tessellationKernel failed at {} line {}\n", __FILE__, __LINE__));
    }

    clFinish(OpenCL::clCommandQueue.value());

    std::size_t originRawData[3] = {0, 0, 0};
    std::size_t regionRawData[3] = {static_cast<std::size_t>(grid_size.x), static_cast<std::size_t>(grid_size.y),
                                    static_cast<std::size_t>(grid_size.z)};
    err = clEnqueueReadImage(OpenCL::clCommandQueue.value(), image, CL_TRUE, originRawData, regionRawData, 0, 0,
                             output_data.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
    }

    clFinish(OpenCL::clCommandQueue.value());

    clReleaseMemObject(inputPos);
    clReleaseMemObject(inputSigma);
    clReleaseMemObject(image);
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::mdspan<int32_t, std::dextents<std::size_t, 3>, std::layout_left> data_view(output_data.data(), grid_size.x,
                                                                                  grid_size.y, grid_size.z);

  std::ofstream myfile;
  myfile.open(framework.name + ".tessellation.gpu.txt");
  std::print(myfile, "# Tesselation using grid-based method\n");
  std::print(myfile, "# GPU Timing: {} [s]\n", timing.count());
  for (std::size_t k = 0; k < static_cast<std::size_t>(grid_size.z); ++k)
  {
    for (std::size_t j = 0; j < static_cast<std::size_t>(grid_size.y); ++j)
    {
      for (std::size_t i = 0; i < static_cast<std::size_t>(grid_size.x); ++i)
      {
        std::print(myfile, "{} {} {} {}\n", i, j, k, data_view[i, j, k]);
      }
    }
  }

  myfile.close();
}

const char* Tessellation::tessellationKernelSource = R"foo(
__kernel void Tessellation(__global float4 *position,
                           __global float *sigma,
                           __write_only image3d_t image,
                           const float4 cella,
                           const float4 cellb,
                           const float4 cellc,
                           const int numberOfAtoms,
                           const int3 numberOfReplicas,
                           const int3 grid_size)
{
  #pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

  int4 ipos = {get_global_id(0), get_global_id(1), get_global_id(2), 0};

  int iatom;
  float value = 0.0f;
  float4 s,t,dr,pos;

  float4 correction = (float4)(1.0/(float)(numberOfReplicas.x), 1.0/(float)(numberOfReplicas.y), 1.0/(float)(numberOfReplicas.z), 0.0f);

  float4 grid_position = correction * (float4)((float)(ipos.x) / (float)(grid_size.x-1), (float)(ipos.y) / (float)(grid_size.y-1), (float)(ipos.z) / (float)(grid_size.z-1), 0.0f);

  int closest_atom = -1;
  float closest_distance = 1e5;


  for(int i=0; i<numberOfReplicas.x; i++)
  {
    for(int j=0; j<numberOfReplicas.y; j++)
    {
      for(int k=0; k<numberOfReplicas.z; k++)
      {
        float4 replicaVector = (float4)(
                    (float)(i)/(float)(numberOfReplicas.x),
                    (float)(j)/(float)(numberOfReplicas.y),
                    (float)(k)/(float)(numberOfReplicas.z),
                    0.0);
        for( iatom = 0; iatom < numberOfAtoms; iatom++)
        {
          dr = grid_position - position[iatom] - replicaVector;

          t = dr - rint(dr);

          dr.x = dot(cella,t);
          dr.y = dot(cellb,t);
          dr.z = dot(cellc,t);
          dr.w = 0.0f;

          float size = sigma[iatom];

          float r = sqrt(dot(dr,dr)) - size;
          if(r < closest_distance)
          {
            closest_distance = r;
            closest_atom = iatom;
          }
        }
      }
    }
  }

  write_imageui(image, ipos, closest_atom);
}
)foo";
