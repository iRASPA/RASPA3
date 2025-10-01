module;
  
#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>              
#include <optional>
#include <limits> 
#include <algorithm>
#include <iostream>
#include <fstream>
#include <tuple> 
#include <chrono>
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

module pore_size_distribution_ban_vlugt;
    
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

BanVlugtPoreSizeDistribution::BanVlugtPoreSizeDistribution(int3 grid_size):
  grid_size(grid_size)
{
  if(OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* PSDShaderSourceCode = BanVlugtPoreSizeDistribution::PSDKernelSource;
    PSDProgram = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &PSDShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));}

    err = clBuildProgram(PSDProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(PSDProgram, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
      std::string message = std::format("Pore size distribution (Ban, Vlugt): OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__, __LINE__, std::string(buffer));
      throw std::runtime_error(message);
    }

    PSDKernel = clCreateKernel(PSDProgram, "PSD", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateKernel failed {} : {}\n", __FILE__ , __LINE__));}

    err = clGetKernelWorkGroupInfo(PSDKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &PSDWorkGroupSize, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));}
  }
}

BanVlugtPoreSizeDistribution::~BanVlugtPoreSizeDistribution()
{
  if(OpenCL::clContext.has_value())
  {
    clReleaseKernel(PSDKernel);
    clReleaseProgram(PSDProgram);
  }
}

void BanVlugtPoreSizeDistribution::run(const ForceField &forceField, const Framework &framework)
{
  double2 probeParameter = double2(10.9 * Units::KelvinToEnergy, 2.64);
  double cutoff = forceField.cutOffFrameworkVDW;
  double3x3 unitCell = framework.simulationBox.cell;
  int3 numberOfReplicas = framework.simulationBox.smallestNumberOfUnitCellsForMinimumImagesConvention(cutoff);
  std::vector<double3> positions = framework.fractionalAtomPositionsUnitCell();
  std::vector<double2> potentialParameters = framework.atomUnitCellLennardJonesPotentialParameters(forceField);
  std::chrono::system_clock::time_point time_begin, time_end;
  cl_int err;

  if(!OpenCL::clContext.has_value()) return;

  time_begin = std::chrono::system_clock::now();

  std::vector<cl_float4> pos(positions.size());
  std::vector<cl_float> sigma(positions.size());

  std::vector<cl_float4> output_data(static_cast<size_t>(grid_size.x * grid_size.y * grid_size.z));

  size_t nGrid = grid_size.x * grid_size.y * grid_size.z;
  std::vector<int> pore_sizes_host(nGrid, std::bit_cast<int>(0.0f));
  cl_mem pore_sizes = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(int) * nGrid, pore_sizes_host.data(), &err);
  if (err != CL_SUCCESS) { throw std::runtime_error("clCreateBuffer for pore_sizes failed"); }

  if (positions.size() > 0)
  { 
    double3 correction = double3(1.0/double(numberOfReplicas.x), 1.0/double(numberOfReplicas.y), 1.0/double(numberOfReplicas.z));
    for(size_t i = 0 ; i < positions.size(); ++i)
    { 
      double3 position = correction * positions[i];
      double2 currentPotentialParameters = potentialParameters[i];
        
      // fill in the fractional position
      pos[i] = {{cl_float(position.x),cl_float(position.y),cl_float(position.z),0.0f}};
          
      // mixing rule for the atom and the probe
      sigma[i] = cl_float(currentPotentialParameters.y);
    }

    cl_mem inputPos = clCreateBuffer(OpenCL::clContext.value(),  CL_MEM_READ_ONLY,  sizeof(float4) * pos.size(), nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputPos, CL_TRUE, 0, sizeof(float4) * pos.size(), pos.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));}

    cl_mem inputSigma = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * sigma.size(), nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputSigma, CL_TRUE, 0, sizeof(cl_float) * sigma.size(), sigma.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));}


    double3x3 replicaCell = double3x3(double(numberOfReplicas.x) * unitCell[0],
                                      double(numberOfReplicas.y) * unitCell[1],
                                      double(numberOfReplicas.z) * unitCell[2]);

    cl_float4 clCella = {{cl_float(replicaCell[0][0]), cl_float(replicaCell[1][0]), cl_float(replicaCell[2][0]), cl_float(0.0)}};
    cl_float4 clCellb = {{cl_float(replicaCell[0][1]), cl_float(replicaCell[1][1]), cl_float(replicaCell[2][1]), cl_float(0.0)}};
    cl_float4 clCellc = {{cl_float(replicaCell[0][2]), cl_float(replicaCell[1][2]), cl_float(replicaCell[2][2]), cl_float(0.0)}};

    cl_int numberOfAtoms = static_cast<cl_int>(positions.size());
    err  = clSetKernelArg(PSDKernel,  0, sizeof(cl_mem), &inputPos);
    err |= clSetKernelArg(PSDKernel,  1, sizeof(cl_mem), &inputSigma);
    err |= clSetKernelArg(PSDKernel,  2, sizeof(cl_mem), &pore_sizes);
    err |= clSetKernelArg(PSDKernel,  3, sizeof(cl_float4), &clCella);
    err |= clSetKernelArg(PSDKernel,  4, sizeof(cl_float4), &clCellb);
    err |= clSetKernelArg(PSDKernel,  5, sizeof(cl_float4), &clCellc);
    err |= clSetKernelArg(PSDKernel,  6, sizeof(cl_int), &numberOfAtoms);
    err |= clSetKernelArg(PSDKernel,  7, sizeof(cl_int3), &numberOfReplicas);
    err |= clSetKernelArg(PSDKernel,  8, sizeof(cl_int3), &grid_size);

    size_t global_work_size[3] = {static_cast<size_t>(grid_size.x), static_cast<size_t>(grid_size.y), static_cast<size_t>(grid_size.z)};
    err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), PSDKernel, 3, nullptr, global_work_size, nullptr, 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel PSDKernel failed at {} line {}\n", __FILE__, __LINE__));}

    clFinish(OpenCL::clCommandQueue.value());

    size_t originRawData[3] = {0,0,0};
    size_t regionRawData[3] = {static_cast<size_t>(grid_size.x), static_cast<size_t>(grid_size.y), static_cast<size_t>(grid_size.z)};

    err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), pore_sizes, CL_TRUE, 0,
                         sizeof(int) * nGrid, pore_sizes_host.data(), 0, nullptr, nullptr);

    std::vector<float> pore_sizes_float(nGrid);
    for (size_t i = 0; i < nGrid; ++i)
      pore_sizes_float[i] = *reinterpret_cast<float*>(&pore_sizes_host[i]);

    clFinish(OpenCL::clCommandQueue.value());

    clReleaseMemObject(inputPos);
    clReleaseMemObject(inputSigma);
    clReleaseMemObject(pore_sizes);
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::mdspan<cl_float4, std::dextents<size_t, 3>, std::layout_left> data_view(
      output_data.data(),  grid_size.x, grid_size.y, grid_size.z);

  std::ofstream myfile;
  myfile.open(framework.name + ".PSD_ban_vlugt.gpu.txt");
  std::print(myfile, "# Pore size distribution using Ban, Vlugt method\n");
  std::print(myfile, "# GPU Timing: {} [s]\n", timing.count());
  for (size_t k = 0; k < static_cast<size_t>(grid_size.z); ++k)
  {
    for (size_t j = 0; j < static_cast<size_t>(grid_size.y); ++j)
    {
      for (size_t i = 0; i < static_cast<size_t>(grid_size.x); ++i)
      {
        size_t idx = i + j * grid_size.x + k * grid_size.x * grid_size.y;
        float pore_size = std::bit_cast<float>(pore_sizes_host[idx]);
        std::print(myfile,"{} {} {} {}\n", i, j, k, pore_size);
      }
    }
  }

  std::cout << "PSD completed in " << timing.count() << " seconds.\n";
  std::cout << "PSD results saved to " << framework.name << ".PSD_ban_vlugt.gpu.txt\n";

  myfile.close();

}

const char* BanVlugtPoreSizeDistribution::PSDKernelSource = R"foo(
__kernel void PSD(
    __global float4 *position,
    __global float *sigma,
    __global float *pore_sizes, // 3D grid flattened: [z*Y*X + y*X + x]
    const float4 cella,
    const float4 cellb,
    const float4 cellc,
    const int numberOfAtoms,
    const int3 numberOfReplicas,
    const int3 grid_size)
{
    int ix = get_global_id(0);
    int iy = get_global_id(1);
    int iz = get_global_id(2);

    int grid_idx = iz * grid_size.y * grid_size.x + iy * grid_size.x + ix;

    float4 correction = (float4)(1.0f/(float)(numberOfReplicas.x), 1.0f/(float)(numberOfReplicas.y), 1.0f/(float)(numberOfReplicas.z), 0.0f);

    float4 grid_position = correction * (float4)(
        (float)(ix) / (float)(grid_size.x-1),
        (float)(iy) / (float)(grid_size.y-1),
        (float)(iz) / (float)(grid_size.z-1),
        0.0f);

    // Step 1: Calculate distance to closest atom
    float closest_distance = 1e5f;
    for(int i=0; i<numberOfReplicas.x; i++)
    for(int j=0; j<numberOfReplicas.y; j++)
    for(int k=0; k<numberOfReplicas.z; k++)
    {
        float4 replicaVector = (float4)(
            (float)(i)/(float)(numberOfReplicas.x),
            (float)(j)/(float)(numberOfReplicas.y),
            (float)(k)/(float)(numberOfReplicas.z),
            0.0f);
        for(int iatom = 0; iatom < numberOfAtoms; iatom++)
        {
            float4 dr = grid_position - position[iatom] - replicaVector;
            float4 t;
            t.x = dr.x - rint(dr.x);
            t.y = dr.y - rint(dr.y);
            t.z = dr.z - rint(dr.z);
            t.w = 0.0f;
            dr.x = dot(cella,t);
            dr.y = dot(cellb,t);
            dr.z = dot(cellc,t);
            dr.w = 0.0f;

          // hotfix to overwrite values of sigma for testing
          // float size = 0.0f;
          // if (sigma[iatom] == 2.30f){
          //   size = 0.0f;
          // }
          // else if (sigma[iatom] == 3.30f){
          //   size = 1.35f;     
          // }

            float size = sigma[iatom];
            float r = sqrt(dot(dr,dr)) - size;
            if(r < closest_distance)
                closest_distance = r;
        }
    }

    // Step 2: For all grid points inside the sphere of radius closest_distance, update their pore size if this value is larger
    float len_a = sqrt(cella.x * cella.x + cella.y * cella.y + cella.z * cella.z);
    float len_b = sqrt(cellb.x * cellb.x + cellb.y * cellb.y + cellb.z * cellb.z);
    float len_c = sqrt(cellc.x * cellc.x + cellc.y * cellc.y + cellc.z * cellc.z);

    float dx = len_a / (float)(grid_size.x - 1);
    float dy = len_b / (float)(grid_size.y - 1);
    float dz = len_c / (float)(grid_size.z - 1);

    int num_points_x = (int)(closest_distance / dx);
    int num_points_y = (int)(closest_distance / dy);
    int num_points_z = (int)(closest_distance / dz);

    for (int dz_ = -num_points_z; dz_ <= num_points_z; dz_++) {
        int zz = iz + dz_;
        if (zz < 0 || zz >= grid_size.z) continue;
        for (int dy_ = -num_points_y; dy_ <= num_points_y; dy_++) {
            int yy = iy + dy_;
            if (yy < 0 || yy >= grid_size.y) continue;
            for (int dx_ = -num_points_x; dx_ <= num_points_x; dx_++) {
                int xx = ix + dx_;
                if (xx < 0 || xx >= grid_size.x) continue;

                // Compute Cartesian offset from center
                float off_x = dx_ * dx;
                float off_y = dy_ * dy;
                float off_z = dz_ * dz;
                float dist = sqrt(off_x*off_x + off_y*off_y + off_z*off_z);
                if (dist > closest_distance) continue;

                int neighbor_idx = zz * grid_size.y * grid_size.x + yy * grid_size.x + xx;

                // Atomically update the pore size if this value is larger
                int new_val = as_int(closest_distance);
                int old_val = pore_sizes[neighbor_idx];
                int assumed;
                do {
                    assumed = old_val;
                    float assumed_float = as_float(assumed);
                    if (closest_distance > assumed_float)
                        old_val = atomic_cmpxchg(&pore_sizes[neighbor_idx], assumed, new_val);
                    else
                        break;
                } while (old_val != assumed);
            }
        }
    }
}
)foo";
