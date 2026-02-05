module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <print>
#include <string>
#include <vector>
#include <numeric>
#include <chrono>
#include <random>
#endif

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module mc_opencl_surface_area;

#ifdef USE_STD_IMPORT
import std;
#endif

import opencl;
import double2;
import double3;
import float4;
import double3x3;
import skspacegroupdatabase;
import atom;
import randomnumbers;
import framework;
import forcefield;
import units;

MC_OpenCL_SurfaceArea::MC_OpenCL_SurfaceArea()
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* surfaceAreaShaderSourceCode = MC_OpenCL_SurfaceArea::surfaceAreaKernelSource;
    surfaceAreaProgram =
        clCreateProgramWithSource(OpenCL::clContext.value(), 1, &surfaceAreaShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
    }

    err = clBuildProgram(surfaceAreaProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(surfaceAreaProgram, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer,
                            &len);
      std::string message =
          std::format("MC_OpenCL_SurfaceArea: OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__,
                      __LINE__, std::string(buffer));
      throw std::runtime_error(message);
    }

    surfaceAreaKernel = clCreateKernel(surfaceAreaProgram, "ComputeSurfaceArea", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateKernel failed {} : {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(surfaceAreaKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &surfaceAreaWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));
    }
  }
}

void MC_OpenCL_SurfaceArea::run(const ForceField &forceField, const Framework &framework, double wellDepthFactor, std::string probePseudoAtom,
                                std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps) const
{
  RandomNumber random{std::nullopt};
  cl_int err;
  std::chrono::system_clock::time_point time_begin, time_end;
  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);

  if (!probeType.has_value())
  {
    throw std::runtime_error(std::format("MC_SurfaceArea: Unknown probe-atom type\n"));
  }

  std::size_t number_of_iterations = numberOfIterations.value_or(100);
  std::size_t number_of_inner_steps = numberOfInnerSteps.value_or(1000);

  double3x3 unit_cell = framework.simulationBox.cell;
  double3x3 inverse_unit_cell = framework.simulationBox.inverseCell;
  std::vector<double3> positions = framework.cartesianAtomPositionsUnitCell();

  std::size_t numberOfAtoms = positions.size();
  size_t global_work_size = (numberOfAtoms + surfaceAreaWorkGroupSize - 1) & ~(surfaceAreaWorkGroupSize - 1);

  std::vector<cl_float4> pos(numberOfAtoms);
  std::vector<cl_float> sigma(numberOfAtoms);

  std::vector<cl_float> output(numberOfAtoms);

  time_begin = std::chrono::system_clock::now();

  for (size_t i = 0; i < numberOfAtoms; i++)
  {
    std::size_t atomType = static_cast<std::size_t>(framework.unitCellAtoms[i].type);
    double size_parameter = forceField(probeType.value(), atomType).sizeParameter();
    double equilibrium_distance = wellDepthFactor * size_parameter;

    // fill in the Cartesian position
    double3 position = positions[i];
    pos[i] = {{cl_float(position.x), cl_float(position.y), cl_float(position.z), 0.0f}};

    // mixing rule for the atom and the probe
    sigma[i] = cl_float(equilibrium_distance);
  }

  // upload position array
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

  // upload equilibrium size array
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

  // setup output array
  cl_mem outputArray =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, sizeof(cl_float) * output.size(), nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
  }

  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), outputArray, CL_TRUE, 0,
                             sizeof(cl_float) * output.size(), output.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
  }


  // setup unit cell matrix
  cl_float4 clCella = {{cl_float(unit_cell[0][0]), cl_float(unit_cell[1][0]), cl_float(unit_cell[2][0]), cl_float(0.0)}};
  cl_float4 clCellb = {{cl_float(unit_cell[0][1]), cl_float(unit_cell[1][1]), cl_float(unit_cell[2][1]), cl_float(0.0)}};
  cl_float4 clCellc = {{cl_float(unit_cell[0][2]), cl_float(unit_cell[1][2]), cl_float(unit_cell[2][2]), cl_float(0.0)}};

  cl_float4 inverse_clCella = {{cl_float(inverse_unit_cell[0][0]), cl_float(inverse_unit_cell[1][0]), cl_float(inverse_unit_cell[2][0]), cl_float(0.0)}};
  cl_float4 inverse_clCellb = {{cl_float(inverse_unit_cell[0][1]), cl_float(inverse_unit_cell[1][1]), cl_float(inverse_unit_cell[2][1]), cl_float(0.0)}};
  cl_float4 inverse_clCellc = {{cl_float(inverse_unit_cell[0][2]), cl_float(inverse_unit_cell[1][2]), cl_float(inverse_unit_cell[2][2]), cl_float(0.0)}};

  cl_int end_index = cl_int(numberOfAtoms);
  cl_int number_of_slices = cl_int(number_of_inner_steps);

  err = clSetKernelArg(surfaceAreaKernel, 0, sizeof(cl_mem), &inputPos);
  err |= clSetKernelArg(surfaceAreaKernel, 1, sizeof(cl_mem), &inputSigma);
  err |= clSetKernelArg(surfaceAreaKernel, 3, sizeof(cl_mem), &outputArray);
  err |= clSetKernelArg(surfaceAreaKernel, 4, sizeof(cl_int), &end_index);
  err |= clSetKernelArg(surfaceAreaKernel, 5, sizeof(cl_int), &number_of_slices);
  err |= clSetKernelArg(surfaceAreaKernel, 6, sizeof(cl_float4), &clCella);
  err |= clSetKernelArg(surfaceAreaKernel, 7, sizeof(cl_float4), &clCellb);
  err |= clSetKernelArg(surfaceAreaKernel, 8, sizeof(cl_float4), &clCellc);
  err |= clSetKernelArg(surfaceAreaKernel, 9, sizeof(cl_float4), &inverse_clCella);
  err |= clSetKernelArg(surfaceAreaKernel, 10, sizeof(cl_float4), &inverse_clCellb);
  err |= clSetKernelArg(surfaceAreaKernel, 11, sizeof(cl_float4), &inverse_clCellc);

  double accumulated_surface_area{};
  for(std::size_t i = 0; i < number_of_iterations; ++i)
  {
    // upload random Cartesian positions
    std::size_t number_of_random_unit_vectors{number_of_inner_steps};
    std::vector<cl_float4> random_unit_vectors(number_of_random_unit_vectors);
    for (size_t j = 0; j < number_of_random_unit_vectors; j++)
    {
      double3 vec = random.randomVectorOnUnitSphere();
      random_unit_vectors[j] = {{cl_float(vec.x), cl_float(vec.y), cl_float(vec.z), 0.0f}};
    }
    cl_mem random_unit_vectors_mem =
          clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(float4) *  random_unit_vectors.size(), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), random_unit_vectors_mem, CL_TRUE, 0, sizeof(float4) * random_unit_vectors.size(),
                                 random_unit_vectors.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
    }

    err |= clSetKernelArg(surfaceAreaKernel, 2, sizeof(cl_mem), &random_unit_vectors_mem);
    err |= clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), surfaceAreaKernel, 1, nullptr, &global_work_size,
                                  &surfaceAreaWorkGroupSize, 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel failed (err: {}) {} : {}\n", err, __FILE__, __LINE__));
    }


    // Read the buffer back to the array
    err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), outputArray, CL_TRUE, 0, output.size() * sizeof(float),
                              output.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("MC_OpenCL_SurfaceArea: error in clEnqueueReadBuffer");
    }

    double surface_area = std::accumulate(output.begin(), output.end(), 0.0);
    accumulated_surface_area += surface_area;

    clReleaseMemObject(random_unit_vectors_mem);
  }

  clFinish(OpenCL::clCommandQueue.value());


  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".mc.sa.gpu.txt");
  std::print(myfile, "# Surface area using Monte Carlo-based method\n");
  std::print(myfile, "# Framework: {}\n", framework.name);
  std::print(myfile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(myfile, "# Space-group Hall-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HallString());
  std::print(myfile, "# Space-group HM-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Space-group IT number: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].number());
  std::print(myfile, "# Probe atom: {} well-depth-factor: {} sigma: {}\n", probePseudoAtom, wellDepthFactor, forceField[probeType.value()].sizeParameter());
  std::print(myfile, "# Number of iterations: {}\n", number_of_iterations);
  std::print(myfile, "# Number of framework atoms: {}\n", framework.unitCellAtoms.size());
  std::print(myfile, "# Number of inner-steps: {}\n", number_of_inner_steps);
  std::print(myfile, "# GPU Timing: {} [s]\n", timing.count());
  myfile << accumulated_surface_area / static_cast<double>(number_of_iterations)  << " [A^2]" << std::endl;
  myfile << (accumulated_surface_area / static_cast<double>(number_of_iterations)) * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant /
                framework.unitCellMass
         << " [m^2/g]" << std::endl;
  myfile << 1.0e4 * (accumulated_surface_area / static_cast<double>(number_of_iterations)) / framework.simulationBox.volume << " [m^2/cm^3]" << std::endl;

  myfile.close();
}

const char* MC_OpenCL_SurfaceArea::surfaceAreaKernelSource = R"foo(
__kernel void ComputeSurfaceArea(__global float4 *position,
                                 __global float *sigma,
                                 __global float4 *randomCartesianPositions,
                                 __global float *output,
                                 const int numberOfAtoms,
                                 const int numberOfSlices,
                                 const float4 cella,
                                 const float4 cellb,
                                 const float4 cellc,
                                 const float4 inverse_cella,
                                 const float4 inverse_cellb,
                                 const float4 inverse_cellc)
{
  float4 dr, ds;
  int iatom = get_global_id(0);
  float counted = 0.0f;
  float total = 0.0f;

  if(iatom < numberOfAtoms)
  {
    float radius_i = sigma[iatom];
    float4 sphere_center = position[iatom];

    for(int slice = 0; slice < numberOfSlices; ++slice)
    {
      float4 unit_vector = randomCartesianPositions[slice];

      // check overlap with other atoms
      bool overlap = false;
      for(int jatom = 0; jatom < numberOfAtoms; ++jatom)
      {
        if(jatom != iatom)
        {
          float4 dr = (sphere_center + radius_i * unit_vector) - position[jatom];

          ds.x = dot(inverse_cella, dr);
          ds.y = dot(inverse_cellb, dr);
          ds.z = dot(inverse_cellc, dr);
          ds.w = 0.0f;

          float4 t = ds - rint(ds);

          dr.x = dot(cella, t);
          dr.y = dot(cellb, t);
          dr.z = dot(cellc, t);
          dr.w = 0.0f;

          float rr = dot(dr, dr);

          if(rr < sigma[jatom] * sigma[jatom])
          {
            overlap = true;
            break;
          }
        }
      }

      if(!overlap)
      {
        counted += 1.0f;
      }

      total += 1.0f;
    }

    output[ iatom ] = (counted / total) * 4.0 * M_PI * radius_i * radius_i;
  }
}
)foo";
