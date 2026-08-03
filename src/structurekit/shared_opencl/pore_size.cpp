module;

#define CL_TARGET_OPENCL_VERSION 120
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module grid_pore_size_opencl;

import std;

import opencl;
import uint3;
import double3;
import double3x3;
import unit_cell;


std::vector<float> poreRadiusFieldOpenCL(uint3 gridSize, const UnitCell &unitCell,
                                         std::span<const float> distance, double slack)
{
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value() || !OpenCL::clCommandQueue.has_value())
  {
    throw std::runtime_error("poreRadiusFieldOpenCL: no OpenCL device found\n");
  }

  const std::int32_t nx = static_cast<std::int32_t>(gridSize.x);
  const std::int32_t ny = static_cast<std::int32_t>(gridSize.y);
  const std::int32_t nz = static_cast<std::int32_t>(gridSize.z);
  const std::size_t numberOfVoxels = distance.size();

  std::vector<float> poreRadius(numberOfVoxels, 0.0f);
  if (numberOfVoxels == 0) return poreRadius;

  std::vector<cl_int> centres;
  centres.reserve(numberOfVoxels / 2);
  for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
  {
    const float radius = distance[voxel];
    if (radius > 0.0f && radius != std::numeric_limits<float>::max())
    {
      centres.push_back(static_cast<cl_int>(voxel));
    }
  }
  if (centres.empty()) return poreRadius;

  cl_int err;
  const char *source = poreSizeKernelSource;
  cl_program program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &source, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("poreRadiusFieldOpenCL: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
  }

  err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    std::size_t length;
    char buffer[2048];
    clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
    clReleaseProgram(program);
    throw std::runtime_error(
        std::format("poreRadiusFieldOpenCL: OpenCL failed to build program (error: {})\n", std::string(buffer)));
  }

  cl_kernel kernel = clCreateKernel(program, "PoreSize", &err);
  if (err != CL_SUCCESS)
  {
    clReleaseProgram(program);
    throw std::runtime_error(std::format("poreRadiusFieldOpenCL: OpenCL clCreateKernel failed at {}\n", __LINE__));
  }

  std::vector<cl_float> field(numberOfVoxels);
  for (std::size_t v = 0; v < numberOfVoxels; ++v) field[v] = distance[v];

  cl_int centreError, distanceError, poreError;
  cl_mem centreBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_int) * centres.size(),
                                       nullptr, &centreError);
  cl_mem distanceBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * numberOfVoxels,
                                         nullptr, &distanceError);
  cl_mem poreBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, sizeof(cl_float) * numberOfVoxels,
                                     nullptr, &poreError);
  if (centreError != CL_SUCCESS || distanceError != CL_SUCCESS || poreError != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("poreRadiusFieldOpenCL: OpenCL clCreateBuffer failed at {}\n", __LINE__));
  }

  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), centreBuffer, CL_TRUE, 0, sizeof(cl_int) * centres.size(),
                             centres.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), distanceBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * numberOfVoxels, field.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), poreBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * numberOfVoxels, poreRadius.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("poreRadiusFieldOpenCL: OpenCL clEnqueueWriteBuffer failed at {}\n", __LINE__));
  }

  const double3x3 cell = unitCell.cell;
  cl_float4 cellRow0 = {{cl_float(cell[0][0]), cl_float(cell[1][0]), cl_float(cell[2][0]), 0.0f}};
  cl_float4 cellRow1 = {{cl_float(cell[0][1]), cl_float(cell[1][1]), cl_float(cell[2][1]), 0.0f}};
  cl_float4 cellRow2 = {{cl_float(cell[0][2]), cl_float(cell[1][2]), cl_float(cell[2][2]), 0.0f}};

  // How many grid steps along each axis a ball of a given radius can possibly reach. The perpendicular width
  // is the right measure for an oblique cell: a step along one axis moves that far towards the opposite face
  // at least, so the box built from it holds the ball.
  double3 widths = unitCell.perpendicularWidths();
  cl_float4 stepsPerLength = {{cl_float(static_cast<double>(nx) / widths.x),
                               cl_float(static_cast<double>(ny) / widths.y),
                               cl_float(static_cast<double>(nz) / widths.z), 0.0f}};

  cl_int numberOfCentresArgument = static_cast<cl_int>(centres.size());
  cl_int3 clGridSize = {{nx, ny, nz, 0}};
  cl_float clSlack = static_cast<cl_float>(slack);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &centreBuffer);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &distanceBuffer);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &poreBuffer);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_float4), &cellRow0);
  err |= clSetKernelArg(kernel, 4, sizeof(cl_float4), &cellRow1);
  err |= clSetKernelArg(kernel, 5, sizeof(cl_float4), &cellRow2);
  err |= clSetKernelArg(kernel, 6, sizeof(cl_float4), &stepsPerLength);
  err |= clSetKernelArg(kernel, 7, sizeof(cl_float), &clSlack);
  err |= clSetKernelArg(kernel, 8, sizeof(cl_int), &numberOfCentresArgument);
  err |= clSetKernelArg(kernel, 9, sizeof(cl_int3), &clGridSize);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("poreRadiusFieldOpenCL: OpenCL clSetKernelArg failed at {}\n", __LINE__));
  }

  std::size_t globalWorkSize = centres.size();
  err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 1, nullptr, &globalWorkSize, nullptr, 0, nullptr,
                               nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("poreRadiusFieldOpenCL: OpenCL clEnqueueNDRangeKernel failed at {}\n", __LINE__));
  }

  err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), poreBuffer, CL_TRUE, 0, sizeof(cl_float) * numberOfVoxels,
                            poreRadius.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("poreRadiusFieldOpenCL: OpenCL clEnqueueReadBuffer failed at {}\n", __LINE__));
  }

  clFinish(OpenCL::clCommandQueue.value());

  clReleaseMemObject(centreBuffer);
  clReleaseMemObject(distanceBuffer);
  clReleaseMemObject(poreBuffer);
  clReleaseKernel(kernel);
  clReleaseProgram(program);

  return poreRadius;
}


const char *poreSizeKernelSource = R"foo(
__kernel void PoreSize(__global const int *centre,
                       __global const float *distance,
                       __global float *poreRadius,
                       const float4 cellRow0,
                       const float4 cellRow1,
                       const float4 cellRow2,
                       const float4 stepsPerLength,
                       const float slack,
                       const int numberOfCentres,
                       const int3 gridSize)
{
  int index = get_global_id(0);
  if(index >= numberOfCentres) return;

  int voxel = centre[index];
  float radius = distance[voxel];
  if(radius <= 0.0f) return;

  int cz = voxel / (gridSize.x * gridSize.y);
  int cy = (voxel / gridSize.x) % gridSize.y;
  int cx = voxel % gridSize.x;

  // A grid point stands for its voxel, so the sphere is allowed to reach a little past its radius before the
  // point is said to be outside it.
  float reach = radius + slack;

  int reachX = (int)(ceil(reach * stepsPerLength.x));
  int reachY = (int)(ceil(reach * stepsPerLength.y));
  int reachZ = (int)(ceil(reach * stepsPerLength.z));

  // Every point the sphere covers is told how big the sphere was, and keeps the largest it hears. A
  // non-negative float compares as its own bit pattern does, so the largest of them can be kept with an
  // integer maximum rather than with a compare-and-swap that has to be retried.
  int bits = as_int(radius);

  for(int dk = -reachZ; dk <= reachZ; dk++)
  {
    for(int dj = -reachY; dj <= reachY; dj++)
    {
      for(int di = -reachX; di <= reachX; di++)
      {
        float4 ds = (float4)((float)(di) / (float)(gridSize.x),
                             (float)(dj) / (float)(gridSize.y),
                             (float)(dk) / (float)(gridSize.z),
                             0.0f);
        float4 dr;
        dr.x = dot(cellRow0, ds);
        dr.y = dot(cellRow1, ds);
        dr.z = dot(cellRow2, ds);
        dr.w = 0.0f;
        if(dot(dr, dr) > reach * reach) continue;

        int ix = ((cx + di) % gridSize.x + gridSize.x) % gridSize.x;
        int iy = ((cy + dj) % gridSize.y + gridSize.y) % gridSize.y;
        int iz = ((cz + dk) % gridSize.z + gridSize.z) % gridSize.z;

        atomic_max((__global int *)(poreRadius + ((iz * gridSize.y + iy) * gridSize.x + ix)), bits);
      }
    }
  }
}
)foo";
