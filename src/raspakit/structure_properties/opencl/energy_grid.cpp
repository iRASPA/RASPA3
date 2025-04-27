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
  
module energy_grid;

import opencl;
import float4;
import double4;

EnergyGrid::EnergyGrid()
{
  if(OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;
    const char* shaderSourceCode = EnergyGrid::energyGridKernel;
    program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &shaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));}

    // Build the program executable
    err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS) 
    {
      size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
      std::string message = std::format("SKComputeIsosurface: OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__, __LINE__, std::string(buffer));
      throw std::runtime_error(message);
    }

    kernel = clCreateKernel(program, "ComputeEnergyGrid", &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateKernel failed {} : {}\n", __FILE__ , __LINE__));}

    err = clGetKernelWorkGroupInfo(kernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &workGroupSize, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));}
  }
};

EnergyGrid::~EnergyGrid()
{
  if(OpenCL::clContext.has_value())
  {
    clReleaseKernel(kernel);
    clReleaseProgram(program);
  }
}


std::vector<cl_float> EnergyGrid::computeEnergyGrid(int3 size, double2 probeParameter,
                                        std::vector<double3> positions, std::vector<double2> potentialParameters,
                                        double3x3 unitCell, int3 numberOfReplicas) noexcept(false)
{

  if(OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    return EnergyGrid::computeEnergyGridGPUImplementation(size, probeParameter, positions, potentialParameters, unitCell, numberOfReplicas);
  }
  return EnergyGrid::computeEnergyGridCPUImplementation(size, probeParameter, positions, potentialParameters, unitCell, numberOfReplicas);
}

std::vector<cl_float> EnergyGrid::computeEnergyGridGPUImplementation(int3 size, double2 probeParameter,
                                                                              std::vector<double3> positions, std::vector<double2> potentialParameters,
                                                                              double3x3 unitCell, int3 numberOfReplicas) noexcept(false)
{
  size_t numberOfAtoms = positions.size();
  size_t temp = static_cast<size_t>(size.x * size.y * size.z);
  cl_int err = 0;

  // make sure the the global work size is an multiple of the work group size
  // (detected on NVIDIA)
  size_t numberOfGridPoints = (temp  + workGroupSize-1) & ~(workGroupSize-1);
  size_t global_work_size = numberOfGridPoints;

  std::vector<cl_float4> pos(numberOfAtoms);
  std::vector<cl_float> epsilon(numberOfAtoms);
  std::vector<cl_float> sigma(numberOfAtoms);

  std::vector<cl_float> outputData = std::vector<cl_float>(numberOfGridPoints);

  std::vector<cl_float4> gridPositions(numberOfGridPoints);
  std::vector<cl_float> output(numberOfGridPoints);

  double3 correction = double3(1.0/double(numberOfReplicas.x), 1.0/double(numberOfReplicas.y), 1.0/double(numberOfReplicas.z));

  if (numberOfAtoms > 0)
  {
    for(size_t i=0 ; i<numberOfAtoms; i++)
    {
      double3 position = correction * positions[i];
      double2 currentPotentialParameters = potentialParameters[i];

      // fill in the Cartesian position
      pos[i] = {{cl_float(position.x),cl_float(position.y),cl_float(position.z),0.0f}};

      // use 4 x epsilon for a probe epsilon of unity
      epsilon[i] = cl_float(4.0*sqrt(currentPotentialParameters.x * probeParameter.x));

      // mixing rule for the atom and the probe
      sigma[i] = cl_float(0.5 * (currentPotentialParameters.y + probeParameter.y));
    }

    size_t index = 0;
    for(int k=0; k<size.z;k++)
    {
      for(int j=0; j<size.y; j++)
      {
        // X various the fastest (contiguous in x)
        for(int i=0 ; i<size.x;i++)
        {
          double3 position = correction * double3(double(i)/double(size.x-1),double(j)/double(size.y-1),double(k)/double(size.z-1));
          gridPositions[index] = {{cl_float(position.x),cl_float(position.y),cl_float(position.z),cl_float(0.0)}};
          index += 1;
        }
      }
    }

    cl_mem inputPos = clCreateBuffer(OpenCL::clContext.value(),  CL_MEM_READ_ONLY,  sizeof(float4) * pos.size(), nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputPos, CL_TRUE, 0, sizeof(float4) * pos.size(), pos.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));}

    cl_mem inputGridPos = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,  sizeof(cl_float4) * gridPositions.size(), nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputGridPos, CL_TRUE, 0, sizeof(cl_float4) * gridPositions.size(), gridPositions.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));}

    cl_mem inputEpsilon = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,  sizeof(cl_float) * epsilon.size(), nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputEpsilon, CL_TRUE, 0, sizeof(cl_float) * epsilon.size(), epsilon.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));}

    cl_mem inputSigma = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * sigma.size(), nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputSigma, CL_TRUE, 0, sizeof(cl_float) * sigma.size(), sigma.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));}

    // set work-item dimensions
    size_t totalNumberOfReplicas = static_cast<size_t>(numberOfReplicas.x * numberOfReplicas.y * numberOfReplicas.z);
    cl_int clNumberOfReplicas =  cl_int(totalNumberOfReplicas);
    std::vector<cl_float4> replicaVector(totalNumberOfReplicas);
    index = 0;
    for(int i=0; i<numberOfReplicas.x; i++)
    {
      for(int j=0; j<numberOfReplicas.y; j++)
      {
        for(int k=0; k<numberOfReplicas.z; k++)
        {
          replicaVector[index] = {
                   {cl_float(double(i)/double(numberOfReplicas.x)),
                   cl_float(double(j)/double(numberOfReplicas.y)),
                   cl_float(double(k)/double(numberOfReplicas.z)),
                   cl_float(0.0)}};
          index += 1;
        }
      }
    }

    // allocate xpos memory and queue it to the device
    cl_mem replicaCellBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,  sizeof(cl_float4) * replicaVector.size(), nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), replicaCellBuffer, CL_TRUE, 0, sizeof(cl_float4) * replicaVector.size(), replicaVector.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));}

    // allocate memory for the output and queue it to the device
    cl_mem outputMemory = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, sizeof(cl_float) * output.size(), nullptr, &err);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));}

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), outputMemory, CL_TRUE, 0, sizeof(cl_float) * output.size(), output.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));}

    double3x3 replicaCell = double3x3(double(numberOfReplicas.x) * unitCell[0],
                                      double(numberOfReplicas.y) * unitCell[1],
                                      double(numberOfReplicas.z) * unitCell[2]);

    cl_float4 clCella = {{cl_float(replicaCell[0][0]), cl_float(replicaCell[1][0]), cl_float(replicaCell[2][0]), cl_float(0.0)}};
    cl_float4 clCellb = {{cl_float(replicaCell[0][1]), cl_float(replicaCell[1][1]), cl_float(replicaCell[2][1]), cl_float(0.0)}};
    cl_float4 clCellc = {{cl_float(replicaCell[0][2]), cl_float(replicaCell[1][2]), cl_float(replicaCell[2][2]), cl_float(0.0)}};

    size_t unitsOfWorkDone = 0;
    size_t sizeOfWorkBatch = 4096;
    while(unitsOfWorkDone < positions.size())
    {
      size_t numberOfAtomsPerThreadgroup = std::min(sizeOfWorkBatch, positions.size()) - unitsOfWorkDone;

      cl_int startIndex = cl_int(unitsOfWorkDone);
      cl_int endIndex = cl_int(unitsOfWorkDone + numberOfAtomsPerThreadgroup);
      err  = clSetKernelArg(kernel,  0, sizeof(cl_mem), &inputPos);
      err |= clSetKernelArg(kernel,  1, sizeof(cl_mem), &inputGridPos);
      err |= clSetKernelArg(kernel,  2, sizeof(cl_mem), &inputEpsilon);
      err |= clSetKernelArg(kernel,  3, sizeof(cl_mem), &inputSigma);
      err |= clSetKernelArg(kernel,  4, sizeof(cl_mem), &replicaCellBuffer);
      err |= clSetKernelArg(kernel,  5, sizeof(cl_mem), &outputMemory);
      err |= clSetKernelArg(kernel,  6, sizeof(cl_int), &clNumberOfReplicas);
      err |= clSetKernelArg(kernel,  7, sizeof(cl_float4), &clCella);
      err |= clSetKernelArg(kernel,  8, sizeof(cl_float4), &clCellb);
      err |= clSetKernelArg(kernel,  9, sizeof(cl_float4), &clCellc);
      err |= clSetKernelArg(kernel,  10, sizeof(cl_int), &startIndex);
      err |= clSetKernelArg(kernel,  11, sizeof(cl_int), &endIndex);
      err |= clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 1, nullptr, &global_work_size, &workGroupSize, 0, nullptr, nullptr);
      if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel failed {} : {}\n", __FILE__, __LINE__));}

      clFinish(OpenCL::clCommandQueue.value());

      unitsOfWorkDone += sizeOfWorkBatch;
    }

    // read output image using SAME size as before
    err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), outputMemory, CL_TRUE, 0, sizeof(cl_float) * outputData.size(), outputData.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));}

    clFinish(OpenCL::clCommandQueue.value());

    outputData.resize(static_cast<size_t>(size.x * size.y * size.z));

    clReleaseMemObject(inputPos);
    clReleaseMemObject(inputGridPos);
    clReleaseMemObject(inputEpsilon);
    clReleaseMemObject(inputSigma);
    clReleaseMemObject(replicaCellBuffer);
    clReleaseMemObject(outputMemory);
  }
  return outputData;
}



// brute-force implementation
std::vector<cl_float> EnergyGrid::computeEnergyGridCPUImplementation(int3 size, double2 probeParameter,
                                                                              std::vector<double3> positions, std::vector<double2> potentialParameters,
                                                                              double3x3 unitCell, int3 numberOfReplicas) noexcept
{

  size_t numberOfAtoms = positions.size();
  size_t temp = static_cast<size_t>(size.x * size.y * size.z);

  std::vector<cl_float> outputData = std::vector<cl_float>(temp);

  double3 correction = double3(1.0/double(numberOfReplicas.x), 1.0/double(numberOfReplicas.y), 1.0/double(numberOfReplicas.z));

  double3x3 replicaCell = double3x3(double(numberOfReplicas.x) * unitCell[0],
                                    double(numberOfReplicas.y) * unitCell[1],
                                    double(numberOfReplicas.z) * unitCell[2]);

  size_t totalNumberOfReplicas = static_cast<size_t>(numberOfReplicas.x * numberOfReplicas.y * numberOfReplicas.z);
  std::vector<double3> replicaVector(totalNumberOfReplicas);
  size_t index = 0;
  for(size_t i=0; i < static_cast<size_t>(numberOfReplicas.x); i++)
  {
    for(size_t j=0; j < static_cast<size_t>(numberOfReplicas.y); j++)
    {
      for(size_t k=0; k < static_cast<size_t>(numberOfReplicas.z); k++)
      {
        replicaVector[index] = double3((double(i)/double(numberOfReplicas.x)),
                                       (double(j)/double(numberOfReplicas.y)),
                                       (double(k)/double(numberOfReplicas.z)));
        index += 1;
      }
    }
  }

  for(size_t z=0; z < static_cast<size_t>(size.z); z++)
  {
    for(size_t y=0; y < static_cast<size_t>(size.y); y++)
    {
      for(size_t x=0 ; x < static_cast<size_t>(size.x); x++)
      {
        double3 gridPosition = correction * double3(double(x)/double(size.x-1),double(y)/double(size.y-1),double(z)/double(size.z-1));

        double value = 0.0;
        for(size_t i=0 ; i<numberOfAtoms; i++)
        {
          double3 position = correction * positions[i];
          double2 currentPotentialParameters = potentialParameters[i];

          // use 4 x epsilon for a probe epsilon of unity
          double epsilon = 4.0*sqrt(currentPotentialParameters.x * probeParameter.x);

          // mixing rule for the atom and the probe
          double sigma = 0.5 * (currentPotentialParameters.y + probeParameter.y);

          for(size_t j=0;j<totalNumberOfReplicas;j++)
          {
            double3 ds = gridPosition - position - replicaVector[j];
            ds.x -= std::rint(ds.x);
            ds.y -= std::rint(ds.y);
            ds.z -= std::rint(ds.z);
            double3 dr = replicaCell * ds;

            double rr = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
            if (rr<12.0*12.0)
            {
              double sigma2rr = sigma*sigma/rr;
              double rri3 = sigma2rr * sigma2rr * sigma2rr;
              value += epsilon*(rri3*(rri3-1.0));
            }
          }
        }

        outputData[x + y * static_cast<size_t>(size.x) + z * static_cast<size_t>(size.x * size.y)] += static_cast<cl_float>(std::min(value,10000000.0));
      }
    }
  }

  return outputData;
}

const char* EnergyGrid::energyGridKernel = R"foo(
__kernel void ComputeEnergyGrid(__global float4 *position,
                                __global float4 *gridposition,
                                __global float *epsilon,
                                __global float *sigma,
                                __global float4 *replicaCell,
                                __global float *output,
                                const int numberOfReplicas,
                                const float4 cella,
                                const float4 cellb,
                                const float4 cellc,
                                const int startIndexAtoms,
                                const int endIndexAtoms)
{
  int igrid = get_global_id(0);
  int lsize = get_local_size(0);
  int lid = get_local_id(0);

  int iatom;
  float value = 0.0f;
  float4 s,t,dr,pos;

  float4 gridpos = gridposition[igrid];

  for(int j=0;j<numberOfReplicas;j++)
  {
    for( iatom = startIndexAtoms; iatom < endIndexAtoms; iatom++)
    {
      pos = position[iatom];

      dr = gridpos - pos - replicaCell[j];

      t = dr - rint(dr);

      dr.x = dot(cella,t);
      dr.y = dot(cellb,t);
      dr.z = dot(cellc,t);
      dr.w = 0.0f;

      float size = sigma[iatom];

      float rr = dot(dr,dr);
      if (rr<12.0f*12.0f)
      {
        float temp = size*size/rr;
        float rri3 = temp * temp * temp;
        value += epsilon[iatom]*(rri3*(rri3-1.0f));
      }
    }
  }

  output[ igrid ] += min(value,10000000.0f);
}
)foo";
