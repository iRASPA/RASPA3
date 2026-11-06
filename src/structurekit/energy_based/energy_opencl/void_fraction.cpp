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

module energy_opencl_void_fraction;

import std;

import opencl;
import uint3;
import float4;
import double4;
import skspacegroupdatabase;
import pair_interactions;
import crystal;
import units;
import blocking_spheres;
import energy_shared_blocking_mask;
#if !(defined(__has_include) && __has_include(<mdspan>))
//import mdspan;
#endif

EnergyOpenCLVoidFraction::EnergyOpenCLVoidFraction()
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* energyGridShaderSourceCode = EnergyOpenCLVoidFraction::energyGridKernelSource;
    energyGridProgram =
        clCreateProgramWithSource(OpenCL::clContext.value(), 1, &energyGridShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
    }

    err = clBuildProgram(energyGridProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(energyGridProgram, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer,
                            &len);
      std::string message =
          std::format("SKComputeIsosurface: OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__,
                      __LINE__, std::string(buffer));
      throw std::runtime_error(message);
    }

    energyGridKernel = clCreateKernel(energyGridProgram, "ComputeEnergyGrid", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateKernel failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clGetKernelWorkGroupInfo(energyGridKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &energyGridWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));
    }

    const char* energyVoidFractionShaderSourceCode = EnergyOpenCLVoidFraction::energyVoidFractionKernelSource;
    energyVoidFractionProgram =
        clCreateProgramWithSource(OpenCL::clContext.value(), 1, &energyVoidFractionShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("clCreateProgramWithSource failed");
    }

    err = clBuildProgram(energyVoidFractionProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("clBuildProgram failed");
    }

    energyVoidFractionKernel = clCreateKernel(energyVoidFractionProgram, "ComputeVoidFraction", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("clCreateKernel failed");
    }

    err = clGetKernelWorkGroupInfo(energyVoidFractionKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &energyVoidFractionWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("clGetKernelWorkGroupInfo failed");
    }
  }
};

EnergyOpenCLVoidFraction::~EnergyOpenCLVoidFraction()
{
  if (OpenCL::clContext.has_value())
  {
    clReleaseKernel(energyGridKernel);
    clReleaseProgram(energyGridProgram);

    clReleaseKernel(energyVoidFractionKernel);
    clReleaseProgram(energyVoidFractionProgram);
  }
}

void EnergyOpenCLVoidFraction::run(const PairInteractions& interactions, const Crystal& framework,
                                   std::string probePseudoAtom, uint3 gridSize, double temperature,
                                   std::span<const BlockingSphere> blockingSpheres, double blockedEnergyPerAngstrom,
                                   double ceiling)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(std::format("EnergyOpenCLVoidFraction: unknown probe-atom type '{}'\n", probePseudoAtom));
  }
  if (temperature <= 0.0)
  {
    throw std::runtime_error("EnergyOpenCLVoidFraction: a Boltzmann average needs a temperature above zero\n");
  }

  double cutoff = interactions.cutOffVDW;
  double3x3 unitCell = framework.unitCell.cell;
  int3 numberOfReplicas = framework.unitCell.smallestNumberOfUnitCellsForMinimumImagesConvention(cutoff);
  std::vector<double3> positions = framework.fractionalPositions;
  std::chrono::steady_clock::time_point time_begin, time_end;

  time_begin = std::chrono::steady_clock::now();

  // Energy-grid computation step
  // ==================================================================================================================================================================

  size_t numberOfAtoms = positions.size();
  size_t temp = static_cast<size_t>(gridSize.x * gridSize.y * gridSize.z);
  cl_int err = 0;

  // The global work size has to be a multiple of the work-group size (detected on NVIDIA). Two kernels run
  // over the same range, so it is padded to a multiple of both their group sizes, which a device may report
  // as any number and not only a power of two: to the least number both divide, by division and not by a
  // mask, which is only a multiple of anything when the unit is a power of two. The padding points are
  // never counted: the energy kernel evaluates them, harmlessly, and the reduction is told how many points
  // are real.
  if (energyGridWorkGroupSize == 0 || energyVoidFractionWorkGroupSize == 0)
  {
    throw std::runtime_error("EnergyOpenCLVoidFraction: a work group of no points\n");
  }
  const size_t paddingUnit = std::lcm(energyGridWorkGroupSize, energyVoidFractionWorkGroupSize);
  size_t numberOfGridPoints = ((temp + paddingUnit - 1) / paddingUnit) * paddingUnit;
  size_t global_work_size = numberOfGridPoints;

  std::vector<cl_float4> pos(numberOfAtoms);
  std::vector<cl_float> epsilon(numberOfAtoms);
  std::vector<cl_float> sigma(numberOfAtoms);
  std::vector<cl_float> shift(numberOfAtoms);

  // std::vector<cl_float> outputData = std::vector<cl_float>(numberOfGridPoints);

  std::vector<cl_float4> gridPositions(numberOfGridPoints);
  std::vector<cl_float> output(numberOfGridPoints);

  double3 correction =
      double3(1.0 / double(numberOfReplicas.x), 1.0 / double(numberOfReplicas.y), 1.0 / double(numberOfReplicas.z));

  if (numberOfAtoms > 0)
  {
    for (size_t i = 0; i < numberOfAtoms; i++)
    {
      double3 position = correction * positions[i];

      // Asking the force field for the pair rather than mixing here means its own rule is the one applied,
      // whichever that is, and the probe named on the command line is the one actually used.
      std::size_t atomType = framework.atoms[i].type;
      double sizeParameter = interactions(probeType.value(), atomType).sizeParameter;
      double strengthParameter = interactions(probeType.value(), atomType).strengthParameter;

      pos[i] = {{cl_float(position.x), cl_float(position.y), cl_float(position.z), 0.0f}};
      epsilon[i] = cl_float(4.0 * strengthParameter);
      sigma[i] = cl_float(sizeParameter);

      // A truncated potential is usually shifted so that it reaches the cutoff at zero rather than stepping
      // there. The force field says whether it is, and leaving the step in puts the energy out by however
      // many neighbours happen to be in range, which is not the same everywhere.
      shift[i] = cl_float(interactions(probeType.value(), atomType).shift);
    }

    size_t index = 0;
    for (std::size_t k = 0; k < gridSize.z; ++k)
    {
      for (std::size_t j = 0; j < gridSize.y; ++j)
      {
        // X various the fastest (contiguous in x)
        for (std::size_t i = 0; i < gridSize.x; ++i)
        {
          // Dividing by the count rather than one less leaves the far face out, which is the near face of the
          // next cell along. Including both counts that plane twice and thins every other one, which biases
          // the average towards whatever happens to lie on an arbitrary plane of the structure.
          double3 position = correction * double3(double(i) / double(gridSize.x), double(j) / double(gridSize.y),
                                                  double(k) / double(gridSize.z));
          gridPositions[index] = {{cl_float(position.x), cl_float(position.y), cl_float(position.z), cl_float(0.0)}};
          ++index;
        }
      }
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

    cl_mem inputGridPos = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                         sizeof(cl_float4) * gridPositions.size(), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputGridPos, CL_TRUE, 0,
                               sizeof(cl_float4) * gridPositions.size(), gridPositions.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    cl_mem inputEpsilon =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * epsilon.size(), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputEpsilon, CL_TRUE, 0,
                               sizeof(cl_float) * epsilon.size(), epsilon.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    cl_mem inputShift =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * shift.size(), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputShift, CL_TRUE, 0, sizeof(cl_float) * shift.size(),
                               shift.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
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

    size_t totalNumberOfReplicas = static_cast<size_t>(numberOfReplicas.x * numberOfReplicas.y * numberOfReplicas.z);
    cl_int clNumberOfReplicas = cl_int(totalNumberOfReplicas);
    std::vector<cl_float4> replicaVector(totalNumberOfReplicas);
    index = 0;
    for (int i = 0; i < numberOfReplicas.x; i++)
    {
      for (int j = 0; j < numberOfReplicas.y; j++)
      {
        for (int k = 0; k < numberOfReplicas.z; k++)
        {
          replicaVector[index] = {{cl_float(double(i) / double(numberOfReplicas.x)),
                                   cl_float(double(j) / double(numberOfReplicas.y)),
                                   cl_float(double(k) / double(numberOfReplicas.z)), cl_float(0.0)}};
          index += 1;
        }
      }
    }

    cl_mem replicaCellBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                              sizeof(cl_float4) * replicaVector.size(), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), replicaCellBuffer, CL_TRUE, 0,
                               sizeof(cl_float4) * replicaVector.size(), replicaVector.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    cl_mem outputMemory =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, sizeof(cl_float) * output.size(), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), outputMemory, CL_TRUE, 0,
                               sizeof(cl_float) * output.size(), output.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
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

    cl_float clCutOffSquared = cl_float(cutoff * cutoff);

    size_t unitsOfWorkDone = 0;
    size_t sizeOfWorkBatch = 4096;
    while (unitsOfWorkDone < positions.size())
    {
      cl_int startIndex = cl_int(unitsOfWorkDone);
      cl_int endIndex = cl_int(std::min(unitsOfWorkDone + sizeOfWorkBatch, positions.size()));
      err = clSetKernelArg(energyGridKernel, 0, sizeof(cl_mem), &inputPos);
      err |= clSetKernelArg(energyGridKernel, 1, sizeof(cl_mem), &inputGridPos);
      err |= clSetKernelArg(energyGridKernel, 2, sizeof(cl_mem), &inputEpsilon);
      err |= clSetKernelArg(energyGridKernel, 3, sizeof(cl_mem), &inputSigma);
      err |= clSetKernelArg(energyGridKernel, 4, sizeof(cl_mem), &replicaCellBuffer);
      err |= clSetKernelArg(energyGridKernel, 5, sizeof(cl_mem), &outputMemory);
      err |= clSetKernelArg(energyGridKernel, 6, sizeof(cl_int), &clNumberOfReplicas);
      err |= clSetKernelArg(energyGridKernel, 7, sizeof(cl_float4), &clCella);
      err |= clSetKernelArg(energyGridKernel, 8, sizeof(cl_float4), &clCellb);
      err |= clSetKernelArg(energyGridKernel, 9, sizeof(cl_float4), &clCellc);
      err |= clSetKernelArg(energyGridKernel, 10, sizeof(cl_int), &startIndex);
      err |= clSetKernelArg(energyGridKernel, 11, sizeof(cl_int), &endIndex);
      err |= clSetKernelArg(energyGridKernel, 12, sizeof(cl_float), &clCutOffSquared);
      err |= clSetKernelArg(energyGridKernel, 13, sizeof(cl_mem), &inputShift);
      err |= clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), energyGridKernel, 1, nullptr, &global_work_size,
                                    &energyGridWorkGroupSize, 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel failed {} : {}\n", __FILE__, __LINE__));
      }

      clFinish(OpenCL::clCommandQueue.value());

      unitsOfWorkDone += sizeOfWorkBatch;
    }

    // The pockets the probe cannot reach, written onto the field the way the backend writes them onto its own.
    // The field is on the device, so it is fetched, marked and sent back; with no spheres it stays where it is.
    if (!blockingSpheres.empty())
    {
      err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), outputMemory, CL_TRUE, 0,
                                sizeof(cl_float) * output.size(), output.data(), 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error(std::format("OpenCL clEnqueueReadBuffer failed {} : {}\n", __FILE__, __LINE__));
      }
      applyBlockingSpheresToEnergy(std::span(output).first(temp), gridSize, framework.unitCell, blockingSpheres,
                                   blockedEnergyPerAngstrom, ceiling);
      err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), outputMemory, CL_TRUE, 0,
                                 sizeof(cl_float) * output.size(), output.data(), 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
      }
    }

    // Accumulation step
    // ==================================================================================================================================================================

    size_t nWorkGroups = numberOfGridPoints / energyVoidFractionWorkGroupSize;

    // Every group writes its partial sum here and nothing is sent down beforehand, so the buffer is the
    // kernel's to write and starts empty.
    std::vector<float> sumReduction(nWorkGroups, 0.0f);
    cl_mem reductionBuffer =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_WRITE_ONLY, nWorkGroups * sizeof(float), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    cl_int clNumberOfPoints = cl_int(temp);
    cl_float clBeta = cl_float(1.0 / (Units::KB * temperature));
    err = clSetKernelArg(energyVoidFractionKernel, 0, sizeof(cl_mem), &outputMemory);
    err |= clSetKernelArg(energyVoidFractionKernel, 1, sizeof(cl_mem), &reductionBuffer);
    err |= clSetKernelArg(energyVoidFractionKernel, 2, energyVoidFractionWorkGroupSize * sizeof(cl_float), nullptr);
    err |= clSetKernelArg(energyVoidFractionKernel, 3, sizeof(cl_float), &clBeta);
    err |= clSetKernelArg(energyVoidFractionKernel, 4, sizeof(cl_int), &clNumberOfPoints);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("SKComputeVoidFraction: error in clSetKernelArg");
    }

    err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), energyVoidFractionKernel, 1, nullptr,
                                 &global_work_size, &energyVoidFractionWorkGroupSize, 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("SKComputeVoidFraction: error in clEnqueueNDRangeKernel");
    }

    // Read the buffer back to the array
    err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), reductionBuffer, CL_TRUE, 0, nWorkGroups * sizeof(float),
                              sumReduction.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("SKComputeVoidFraction: error in clEnqueueReadBuffer");
    }

    clFinish(OpenCL::clCommandQueue.value());

    // Final summation with CPU
    double fraction = 0.0;
    for (size_t i = 0; i < nWorkGroups; i++) fraction += double(sumReduction[i]);
    this->voidFraction = fraction / static_cast<double>(temp);

    time_end = std::chrono::steady_clock::now();

    std::chrono::duration<double> timing = time_end - time_begin;

    std::ofstream myfile;
    myfile.open(framework.name + ".energy.vf.gpu.txt");
    std::print(myfile, "# Void-fraction using energy-based method\n");
    std::print(myfile, "# Crystal: {}\n", framework.name);
    std::print(myfile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
    std::print(myfile, "# Space-group Hall-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HallString());
    std::print(myfile, "# Space-group HM-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
    std::print(myfile, "# Space-group IT number: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].number());
    std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
    std::print(myfile, "# Crystal volume: {} [Å³]\n", framework.unitCell.volume);
    std::print(myfile, "# Crystal mass: {} [g/mol]\n", framework.mass);
    std::print(myfile, "# Crystal density: {} [kg/m³]\n", 1e-3 * framework.mass /
        (framework.unitCell.volume * Units::Angstrom * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant));
    std::print(myfile, "# Probe atom: {} strength/kʙ: {} [K] size: {} [Å]\n", probePseudoAtom,
               interactions[probeType.value()].strengthParameter * Units::EnergyToKelvin,
               interactions[probeType.value()].sizeParameter);
    std::print(myfile, "# Temperature: {} [K]\n", temperature);
    std::print(myfile, "# Blocking spheres: {}\n", blockingSpheres.size());
    std::print(myfile, "# Cutoff: {} [Å], replicas: {}x{}x{}\n", cutoff, numberOfReplicas.x, numberOfReplicas.y,
               numberOfReplicas.z);
    std::print(myfile, "# Grid: {} x {} x {} points, sampled without the far face\n", gridSize.x, gridSize.y,
               gridSize.z);
    std::print(myfile, "# GPU Timing: {} [s]\n", timing.count());
    myfile << fraction / (static_cast<double>(gridSize.x) * static_cast<double>(gridSize.y) *
                          static_cast<double>(gridSize.z))
           << std::endl;
    myfile.close();

    clReleaseMemObject(inputPos);
    clReleaseMemObject(inputGridPos);
    clReleaseMemObject(inputEpsilon);
    clReleaseMemObject(inputSigma);
    clReleaseMemObject(inputShift);
    clReleaseMemObject(replicaCellBuffer);
    clReleaseMemObject(outputMemory);
    clReleaseMemObject(reductionBuffer);
  }
}

const char* EnergyOpenCLVoidFraction::energyGridKernelSource = R"foo(
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
                                const int endIndexAtoms,
                                const float cutOffSquared,
                                __global float *shift)
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

      // A grid point can lie on an atom's centre, where the pair energy has no value. Held off zero the
      // separation gives an energy far above the ceiling, which is what the point is worth. A pair with no
      // strength is skipped outright: its bracket overflows there and nought times that is not a number.
      float rr = fmax(dot(dr,dr), 1.0e-6f);
      if (rr<cutOffSquared && epsilon[iatom] != 0.0f)
      {
        float temp = size*size/rr;
        float rri3 = temp * temp * temp;
        value += epsilon[iatom]*(rri3*(rri3-1.0f)) - shift[iatom];
      }
    }
  }

  output[ igrid ] += min(value,10000000.0f);
}
)foo";

const char* EnergyOpenCLVoidFraction::energyVoidFractionKernelSource = R"foo(
__kernel void ComputeVoidFraction (__global const float *input,
                           __global float *partialSums,
                           __local float *localSums,
                           const float beta,
                           const int numberOfPoints)
{
   size_t local_id = get_local_id(0);
   size_t global_id = get_global_id(0);
   size_t group_size = get_local_size(0);
   size_t group_id = get_group_id(0);

   // Copy from global memory to local memory. The range is padded out to whole groups; a point past the
   // end of the grid is padding and weighs nothing.
   localSums[local_id] = global_id < (size_t)numberOfPoints ? exp(-beta * input[global_id]) : 0.0f;

   // Fold the upper part of what is left onto the lower, until one sum is left. A group may be of any size,
   // so the part kept is the larger half, and the one in the middle of an odd count stays where it is to be
   // folded the next time round. Halving the size and adding element for element only reaches every element
   // when the size is a power of two: a group of twelve comes down to three, and the third is never added.
   for (size_t left = group_size; left > 1; )
   {
      // Waiting for each addition into the given workgroup
      barrier(CLK_LOCAL_MEM_FENCE);

      size_t kept = (left + 1) / 2;
      if (local_id + kept < left)
        localSums[local_id] += localSums[local_id + kept];
      left = kept;
   }

   // Write result into partialSums[nWorkGroups]
   if (local_id == 0)
     partialSums[group_id] = localSums[0];
}
)foo";
