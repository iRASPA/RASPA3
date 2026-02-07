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
#endif

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module integration_opencl_surface_area;

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

Integration_OpenCL_SurfaceArea::Integration_OpenCL_SurfaceArea()
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* surfaceAreaShaderSourceCode = Integration_OpenCL_SurfaceArea::surfaceAreaKernelSource;
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
          std::format("Integration_OpenCL_SurfaceArea: OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__,
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

void Integration_OpenCL_SurfaceArea::run(const ForceField &forceField,
                                const Framework &framework, 
                                double wellDepthFactor,
                                std::string probePseudoAtom,
                                std::optional<std::size_t> numberOfSlices) const
{
  cl_int err;
  std::chrono::system_clock::time_point time_begin, time_end;
  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);

  if (!probeType.has_value())
  {
    throw std::runtime_error(std::format("Integration_SurfaceArea: Unknown probe-atom type\n"));
  }

  std::size_t number_of_slices = numberOfSlices.value_or(1024);

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
  cl_int number_of_slices_parameter = cl_int(number_of_slices);

  err = clSetKernelArg(surfaceAreaKernel, 0, sizeof(cl_mem), &inputPos);
  err |= clSetKernelArg(surfaceAreaKernel, 1, sizeof(cl_mem), &inputSigma);
  err |= clSetKernelArg(surfaceAreaKernel, 2, sizeof(cl_mem), &outputArray);
  err |= clSetKernelArg(surfaceAreaKernel, 3, sizeof(cl_int), &end_index);
  err |= clSetKernelArg(surfaceAreaKernel, 4, sizeof(cl_int), &number_of_slices_parameter);
  err |= clSetKernelArg(surfaceAreaKernel, 5, sizeof(cl_float4), &clCella);
  err |= clSetKernelArg(surfaceAreaKernel, 6, sizeof(cl_float4), &clCellb);
  err |= clSetKernelArg(surfaceAreaKernel, 7, sizeof(cl_float4), &clCellc);
  err |= clSetKernelArg(surfaceAreaKernel, 8, sizeof(cl_float4), &inverse_clCella);
  err |= clSetKernelArg(surfaceAreaKernel, 9, sizeof(cl_float4), &inverse_clCellb);
  err |= clSetKernelArg(surfaceAreaKernel, 10, sizeof(cl_float4), &inverse_clCellc);
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
    throw std::runtime_error("Integration_OpenCL_SurfaceArea: error in clEnqueueReadBuffer");
  }

  clFinish(OpenCL::clCommandQueue.value());

  double accumulated_surface_area = std::accumulate(output.begin(), output.end(), 0.0);

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".integration.sa.gpu.txt");
  std::print(myfile, "# Surface area using Monte Carlo-based method\n");
  std::print(myfile, "# Framework: {}\n", framework.name);
  std::print(myfile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(myfile, "# Space-group Hall-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HallString());
  std::print(myfile, "# Space-group HM-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Space-group IT number: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].number());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.unitCellAtoms.size());
  std::print(myfile, "# Framework volume: {} [Å³]\n", framework.simulationBox.volume);
  std::print(myfile, "# Framework mass: {} [g/mol]\n", framework.unitCellMass);
  std::print(myfile, "# Framework density: {} [kg/m³]\n", 1e-3 * framework.unitCellMass /
      (framework.simulationBox.volume * Units::Angstrom * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant));
  std::print(myfile, "# Probe atom: {} well-depth-factor: {} sigma: {}\n", probePseudoAtom, wellDepthFactor, forceField[probeType.value()].sizeParameter());
  std::print(myfile, "# Number of integration points per atom: {}\n", (number_of_slices + 1) * number_of_slices);
  std::print(myfile, "# GPU Timing: {} [s]\n", timing.count());
  myfile << accumulated_surface_area << " [Å²]" << std::endl;
  myfile << 1.0e4 * accumulated_surface_area / framework.simulationBox.volume << " [m²/cm³]" << std::endl;
  myfile << accumulated_surface_area * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant /
                framework.unitCellMass
         << " [m²/g]" << std::endl;

  myfile.close();
}

const char* Integration_OpenCL_SurfaceArea::surfaceAreaKernelSource = R"foo(
__kernel void ComputeSurfaceArea(__global float4 *position,
                                 __global float *sigma,
                                 __global float *output,
                                 const int numberOfAtoms,
                                 const int number_of_slices,
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

    for(int stackNumber = 0; stackNumber <= number_of_slices; ++stackNumber)
    {
      for(int sliceNumber = 0; sliceNumber < number_of_slices; ++sliceNumber)
      {
        float theta = float(sliceNumber) * 2.0 * M_PI / float(number_of_slices);

        float u = float(stackNumber) / float(number_of_slices);
        float phi = acos(2.0 * u - 1.0);

        float sinTheta = sin(theta);
        float sinPhi = sin(phi);
        float cosTheta = cos(theta);
        float cosPhi = cos(phi);
        float4 unit_vector = (float4)(sinPhi * cosTheta, sinPhi * sinTheta, cosPhi, 0.0f);

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
    }

    output[ iatom ] = (counted / total) * 4.0 * M_PI * radius_i * radius_i;
  }
}
)foo";
