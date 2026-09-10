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

module energy_opencl_surface_area;

import std;

import opencl;
import int3;
import uint3;
import float4;
import double4;
import skspacegroupdatabase;
import pair_interactions;
import crystal;
import units;
import surface_curvature;

import energy_shared_isosurface;
import energy_opencl_lewiner_isosurface;

EnergyOpenCLSurfaceArea::EnergyOpenCLSurfaceArea()
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char *energyGridShaderSourceCode = EnergyOpenCLSurfaceArea::energyGridKernelSource;
    energyGridProgram =
        clCreateProgramWithSource(OpenCL::clContext.value(), 1, &energyGridShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
    }

    err = clBuildProgram(energyGridProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      std::size_t len;
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
                                   sizeof(std::size_t), &energyGridWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));
    }

    const char *energyEnergyOpenCLSurfaceAreaShaderSourceCode =
        EnergyOpenCLSurfaceArea::marchingCubesKernelSource.c_str();
    energyEnergyOpenCLSurfaceAreaProgram = clCreateProgramWithSource(
        OpenCL::clContext.value(), 1, &energyEnergyOpenCLSurfaceAreaShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format(
          "EnergyOpenCLSurfaceArea: OpenCL clCreateProgramWithSource failed at {} line {}\n", __FILE__, __LINE__));
    }

    // Build the program executable
    err = clBuildProgram(energyEnergyOpenCLSurfaceAreaProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      std::size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(energyEnergyOpenCLSurfaceAreaProgram, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG,
                            sizeof(buffer), buffer, &len);
      throw std::runtime_error(std::format("Isorface: OpenCL Failed to build program at {}1 (line {} error: {})\n",
                                           __FILE__, __LINE__, buffer));
    }

    constructHPLevelKernel = clCreateKernel(energyEnergyOpenCLSurfaceAreaProgram, "constructHPLevel", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("OpenCL clCreateProgramWithSource failed at {} line {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(constructHPLevelKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(std::size_t), &constructHPLevelKernelWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));
    }

    classifyCubesKernel = clCreateKernel(energyEnergyOpenCLSurfaceAreaProgram, "classifyCubes", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("OpenCL clCreateProgramWithSource failed at {} line {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(constructHPLevelKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(std::size_t), &classifyCubesKernelWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));
    }

    traverseHPKernel[4] = clCreateKernel(energyEnergyOpenCLSurfaceAreaProgram, "traverseHP16", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(traverseHPKernel[4], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(std::size_t), &traverseHPKernelWorkGroupSize[4], nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("cOpenCL lGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));
    }
    traverseHPKernel[5] = clCreateKernel(energyEnergyOpenCLSurfaceAreaProgram, "traverseHP32", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(traverseHPKernel[5], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(std::size_t), &traverseHPKernelWorkGroupSize[5], nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));
    }
    traverseHPKernel[6] = clCreateKernel(energyEnergyOpenCLSurfaceAreaProgram, "traverseHP64", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(traverseHPKernel[6], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(std::size_t), &traverseHPKernelWorkGroupSize[6], nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));
    }
    traverseHPKernel[7] = clCreateKernel(energyEnergyOpenCLSurfaceAreaProgram, "traverseHP128", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(traverseHPKernel[7], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(std::size_t), &traverseHPKernelWorkGroupSize[7], nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));
    }
    traverseHPKernel[8] = clCreateKernel(energyEnergyOpenCLSurfaceAreaProgram, "traverseHP256", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(traverseHPKernel[8], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(std::size_t), &traverseHPKernelWorkGroupSize[8], nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));
    }
    traverseHPKernel[9] = clCreateKernel(energyEnergyOpenCLSurfaceAreaProgram, "traverseHP512", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateKernel failed at {} line {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(traverseHPKernel[9], OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(std::size_t), &traverseHPKernelWorkGroupSize[9], nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("OpenCL clGetKernelWorkGroupInfo failed at {} line {}\n", __FILE__, __LINE__));
    }
  }
};

EnergyOpenCLSurfaceArea::~EnergyOpenCLSurfaceArea()
{
  if (OpenCL::clContext.has_value())
  {
    clReleaseKernel(energyGridKernel);
    clReleaseProgram(energyGridProgram);

    clReleaseKernel(traverseHPKernel[9]);
    clReleaseKernel(traverseHPKernel[8]);
    clReleaseKernel(traverseHPKernel[7]);
    clReleaseKernel(traverseHPKernel[6]);
    clReleaseKernel(traverseHPKernel[5]);
    clReleaseKernel(traverseHPKernel[4]);
    clReleaseKernel(classifyCubesKernel);
    clReleaseKernel(constructHPLevelKernel);
    clReleaseProgram(energyEnergyOpenCLSurfaceAreaProgram);
  }
}

void EnergyOpenCLSurfaceArea::run(const PairInteractions &interactions, const Crystal &framework, double isoValue,
                         std::string probePseudoAtom, uint3 grid_size)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(std::format("MC_SurfaceArea: Unknown probe-atom type\n"));
  }

  double cutoff = interactions.cutOffVDW;
  double3x3 unitCell = framework.unitCell.cell;
  int3 numberOfReplicas = framework.unitCell.smallestNumberOfUnitCellsForMinimumImagesConvention(cutoff);
  std::vector<double3> positions = framework.fractionalPositions;
  std::vector<double2> potentialParameters = framework.lennardJonesParameters(interactions);
  std::chrono::steady_clock::time_point time_begin, time_end;

  time_begin = std::chrono::steady_clock::now();

  // Energy-grid computation step
  // ==================================================================================================================================================================

  std::size_t numberOfAtoms = positions.size();
  std::size_t temp = static_cast<std::size_t>(grid_size.x * grid_size.y * grid_size.z);
  cl_int err = 0;

  // make sure the the global work size is an multiple of the work group size
  // (detected on NVIDIA)
  std::size_t numberOfGridPoints = (temp + energyGridWorkGroupSize - 1) & ~(energyGridWorkGroupSize - 1);
  std::size_t energy_global_work_size = numberOfGridPoints;

  std::vector<cl_float4> pos(numberOfAtoms);
  std::vector<cl_float> epsilon(numberOfAtoms);
  std::vector<cl_float> sigma(numberOfAtoms);
  std::vector<cl_float> shift(numberOfAtoms);

  std::vector<cl_float4> gridPositions(numberOfGridPoints);
  std::vector<cl_float> output(numberOfGridPoints);

  double3 correction =
      double3(1.0 / double(numberOfReplicas.x), 1.0 / double(numberOfReplicas.y), 1.0 / double(numberOfReplicas.z));

  if (numberOfAtoms == 0)
  {
    return;
  }

  for (std::size_t i = 0; i < numberOfAtoms; i++)
  {
    double3 position = correction * positions[i];

    std::size_t atomType = framework.atoms[i].type;
    double size_parameter = interactions(probeType.value(), atomType).sizeParameter;
    double strength_parameter = interactions(probeType.value(), atomType).strengthParameter;

    // fill in the Cartesian position
    pos[i] = {{cl_float(position.x), cl_float(position.y), cl_float(position.z), 0.0f}};

    // use 4 x epsilon for a probe epsilon of unity
    epsilon[i] = cl_float(4.0 * strength_parameter);

    // mixing rule for the atom and the probe
    sigma[i] = cl_float(size_parameter);

    // A truncated potential is usually shifted so that it reaches the cutoff at zero rather than stepping
    // there. Leaving the step in matters here more than elsewhere: the default iso-surface is the one where
    // the energy is zero, and an unshifted potential puts that surface in the wrong place.
    shift[i] = cl_float(interactions(probeType.value(), atomType).shift);
  }

  std::size_t index = 0;
  for (std::size_t k = 0; k < grid_size.z; ++k)
  {
    for (std::size_t j = 0; j < grid_size.y; ++j)
    {
      // X various the fastest (contiguous in x)
      for (std::size_t i = 0; i < grid_size.x; ++i)
      {
        // Endpoint-exclusive, as everywhere else: the far face of a periodic cell is the near face of the
        // next one along, so including both counts one plane twice.
        double3 position = correction * double3(double(i) / double(grid_size.x), double(j) / double(grid_size.y),
                                                double(k) / double(grid_size.z));
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

  std::size_t totalNumberOfReplicas =
      static_cast<std::size_t>(numberOfReplicas.x * numberOfReplicas.y * numberOfReplicas.z);
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
  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), outputMemory, CL_TRUE, 0, sizeof(cl_float) * output.size(),
                             output.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
  }

  double3x3 replicaCell = double3x3(double(numberOfReplicas.x) * unitCell[0], double(numberOfReplicas.y) * unitCell[1],
                                    double(numberOfReplicas.z) * unitCell[2]);

  cl_float4 clCella = {
      {cl_float(replicaCell[0][0]), cl_float(replicaCell[1][0]), cl_float(replicaCell[2][0]), cl_float(0.0)}};
  cl_float4 clCellb = {
      {cl_float(replicaCell[0][1]), cl_float(replicaCell[1][1]), cl_float(replicaCell[2][1]), cl_float(0.0)}};
  cl_float4 clCellc = {
      {cl_float(replicaCell[0][2]), cl_float(replicaCell[1][2]), cl_float(replicaCell[2][2]), cl_float(0.0)}};

  cl_float clCutOffSquared = cl_float(cutoff * cutoff);

  std::size_t unitsOfWorkDone = 0;
  std::size_t sizeOfWorkBatch = 4096;
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
    err |= clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), energyGridKernel, 1, nullptr,
                                  &energy_global_work_size, &energyGridWorkGroupSize, 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel failed {} : {}\n", __FILE__, __LINE__));
    }

    clFinish(OpenCL::clCommandQueue.value());

    unitsOfWorkDone += sizeOfWorkBatch;
  }

  clReleaseMemObject(inputPos);
  clReleaseMemObject(inputGridPos);
  clReleaseMemObject(inputEpsilon);
  clReleaseMemObject(inputSigma);
  clReleaseMemObject(inputShift);
  clReleaseMemObject(replicaCellBuffer);

  // Hand the energy field to the iso-surface extraction, which no longer cares how the field was made.
  err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), outputMemory, CL_TRUE, 0, sizeof(cl_float) * output.size(),
                            output.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clEnqueueReadBuffer failed at {} line {}\n", __FILE__, __LINE__));
  }
  clFinish(OpenCL::clCommandQueue.value());
  clReleaseMemObject(outputMemory);

  IsosurfaceArea surface = areaOfIsosurface(framework, std::span(output).first(temp), grid_size, isoValue);
  double accumulated_surface_area = surface.area;

  time_end = std::chrono::steady_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".energy.sa.gpu.txt");
  std::print(myfile, "# Surface area using energy-based method\n");
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
  std::print(myfile, "# Probe atom: {} iso-value: {} sigma: {}\n", probePseudoAtom, isoValue, interactions[probeType.value()].sizeParameter);
  std::print(myfile, "# Grid size: {}x{}x{}\n", grid_size.x, grid_size.y, grid_size.z);
  std::print(myfile, "# Triangles: {} (discarded as implausibly large: {})\n", surface.numberOfTriangles,
             surface.numberOfRejectedTriangles);
  std::print(myfile, "# GPU Timing: {} [s]\n", timing.count());
  myfile << accumulated_surface_area << " [A^2]" << std::endl;
  myfile << accumulated_surface_area * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant /
                framework.mass
         << " [m^2/g]" << std::endl;
  myfile << 1.0e4 * accumulated_surface_area / framework.unitCell.volume << " [m^2/cm^3]" << std::endl;
  writeIsosurfaceCurvature(myfile, surface, "Curvature split:");
  myfile.close();
}

std::vector<double3> EnergyOpenCLSurfaceArea::trianglesOfIsosurface(std::span<const float> field, uint3 grid_size,
                                                                    double isoValue,
                                                                    std::vector<double3> *gradients)
{
  // Lewiner's tables and the same face / interior tests as the processor extractor. The
  // histo-pyramid kernels stay in this translation unit but are no longer the `--gpu` path;
  // this walk is not capped at 512³.
  return trianglesOfLewinerIsosurface(field, grid_size, isoValue, gradients);
}


IsosurfaceArea EnergyOpenCLSurfaceArea::areaOfIsosurface(const Crystal &framework, std::span<const float> field,
                                                         uint3 grid_size, double isoValue)
{
  std::vector<double3> gradients;
  std::vector<double3> corners = this->trianglesOfIsosurface(field, grid_size, isoValue, &gradients);
  return accumulateTriangleAreas(framework.unitCell.cell, grid_size, corners, gradients, FieldSense::GrowsIntoSolid);
}

const char *EnergyOpenCLSurfaceArea::energyGridKernelSource = R"foo(
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

      float rr = dot(dr,dr);
      if (rr < cutOffSquared)
      {
        float temp = size * size / rr;
        float rri3 = temp * temp * temp;
        value += epsilon[iatom] * (rri3 * (rri3 - 1.0f)) - shift[iatom];
      }
    }
  }

  output[ igrid ] += min(value,10000000.0f);
}
)foo";

std::string EnergyOpenCLSurfaceArea::marchingCubesKernelSource = std::string(R"foo(

#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_NEAREST;

// A note on the vertex normals the traverseHP kernels write, because there are two extractors here and they
// have to agree.
//
// What is written beside each vertex is the field's gradient there, by central differences on the grid and
// interpolated along the cube edge the vertex sits on, in the *field's own sense*: f(i+1) - f(i-1), so it
// points towards larger values of the field. On an energy field that is into the wall. The processor extractor
// forms it the same way round, so the two need no reconciling and a consumer can be written once.
//
// This used to be stored negated, which happened to point out of the wall on an energy field and so was right
// for lighting a rendered surface, and the host had to negate it back to compare the two extractors. The sense
// now belongs to whoever asks: `FieldSense` on the host says which way the field grows relative to the void and
// turns the gradient into an outward normal accordingly. Anything rendering this buffer has to negate it.
//
// The magnitude differs between the two extractors and always did: this one leaves the difference unscaled,
// twice the true derivative, while the processor one halves it and then normalises. Only the direction is ever
// used, so neither is wrong, and nothing downstream may rely on the length.


// Cube description:
//         7 ________ 6           _____6__             ________
//         /|       /|         7/|       /|          /|       /|
//       /  |     /  |        /  |     /5 |        /  6     /  |
//   4 /_______ /    |      /__4____ /    10     /_______3/    |
//    |     |  |5    |     |    11  |     |     |     |  |   2 |
//    |    3|__|_____|2    |     |__|__2__|     | 4   |__|_____|
//    |    /   |    /      8   3/   9    /      |    /   |    /
//    |  /     |  /        |  /     |  /1       |  /     5  /
//    |/_______|/          |/___0___|/          |/_1_____|/
//   0          1        0          1
//        Nodes                Borders               Faces


__constant int4 cubeOffsets[8] =
{
  {0, 0, 0, 0},
  {1, 0, 0, 0},
  {0, 0, 1, 0},
  {1, 0, 1, 0},
  {0, 1, 0, 0},
  {1, 1, 0, 0},
  {0, 1, 1, 0},
  {1, 1, 1, 0}
};

__constant char  offsets3[72] =
{
  // 0
  0,0,0,
  1,0,0,
  // 1
  1,0,0,
  1,0,1,
  // 2
  1,0,1,
  0,0,1,
  // 3
  0,0,1,
  0,0,0,
  // 4
  0,1,0,
  1,1,0,
  // 5
  1,1,0,
  1,1,1,
  // 6
  1,1,1,
  0,1,1,
  // 7
  0,1,1,
  0,1,0,
  // 8
  0,0,0,
  0,1,0,
  // 9
  1,0,0,
  1,1,0,
  // 10
  1,0,1,
  1,1,1,
  // 11
  0,0,1,
  0,1,1
};


// Look up table for the number of triangles produced for each of the 256 cases. There at most 5 triangular facets necessary.
__constant uchar numberOfTriangles[256] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 2, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 3, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 3, 2, 3, 3, 2, 3, 4, 4, 3, 3, 4, 4, 3, 4, 5, 5, 2, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 3, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 4, 2, 3, 3, 4, 3, 4, 2, 3, 3, 4, 4, 5, 4, 5, 3, 2, 3, 4, 4, 3, 4, 5, 3, 2, 4, 5, 5, 4, 5, 2, 4, 1, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 3, 2, 3, 3, 4, 3, 4, 4, 5, 3, 2, 4, 3, 4, 3, 5, 2, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 4, 3, 4, 4, 3, 4, 5, 5, 4, 4, 3, 5, 2, 5, 4, 2, 1, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 2, 3, 3, 2, 3, 4, 4, 5, 4, 5, 5, 2, 4, 3, 5, 4, 3, 2, 4, 1, 3, 4, 4, 5, 4, 5, 3, 4, 4, 5, 5, 2, 3, 4, 2, 1, 2, 3, 3, 2, 3, 4, 2, 1, 3, 2, 4, 1, 2, 1, 1, 0};

)foo") +

                                                                 std::string(R"foo(

// The last part of the algorithm involves forming the correct facets from the positions that the isosurface intersects the edges of the grid cell.
// Again a table (by Cory Gene Bloyd) is used which this time uses the same cubeindex but allows the vertex sequence to be looked up for as many triangular
// facets are necessary to represent the isosurface within the grid cell.
__constant int triTable[4096] =
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1,
  3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1,
  3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1,
  3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1,
  9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1,
  9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1,
  2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1,
  8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1,
  9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1,
  4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1,
  3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1,
  1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1,
  4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1,
  4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1,
  9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1,
  5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1,
  2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1,
  9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1,
  0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1,
  2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1,
  10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1,
  4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1,
  5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1,
  5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1,
  9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1,
  0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1,
  1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1,
  10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1,
  8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1,
  2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1,
  7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1,
  9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1,
  2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1,
  11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1,
  9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1,
  5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1,
  11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1,
  11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1,
  1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1,
  9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1,
  5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1,
  2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1,
  0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1,
  5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1,
  6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1,
  3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1,
  6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1,
  5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1,
  1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1,
  10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1,
  6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1,
  8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1,
  7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1,
  3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1,
  5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1,
  0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1,
  9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1,
  8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1,
  5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1,
  0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1,
  6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1,
  10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1,
  10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1,
  8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1,
  1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1,
  3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1,
  0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1,
  10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1,
  3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1,
  6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1,
  9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1,
  8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1,
  3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1,
  6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1,
  0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1,
  10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1,
  10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1,
  2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1,
  7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1,
  7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1,
  2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1,
  1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1,
  11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1,
  8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1,
  0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1,
  7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1,
  10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1,
  2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1,
  6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1,
  7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1,
  2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1,
  1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1,
  10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1,
  10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1,
  0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1,
  7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1,
  6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1,
  8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1,
  9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1,
  6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1,
  4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1,
  10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1,
  8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1,
  0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1,
  1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1,
  8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1,
  10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1,
  4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1,
  10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1,
  5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1,
  11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1,
  9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1,
  6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1,
  7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1,
  3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1,
  7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1,
  9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1,
  3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1,
  6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1,
  9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1,
  1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1,
  4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1,
  7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1,
  6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1,
  3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1,
  0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1,
  6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1,
  0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1,
  11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1,
  6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1,
  5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1,
  9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1,
  1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1,
  1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1,
  10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1,
  0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1,
  5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1,
  10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1,
  11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1,
  9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1,
  7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1,
  2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1,
  8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1,
  9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1,
  9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1,
  1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1,
  9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1,
  9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1,
  5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1,
  0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1,
  10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1,
  2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1,
  0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1,
  0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1,
  9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1,
  5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1,
  3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1,
  5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1,
  8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1,
  0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1,
  9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1,
  0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1,
  1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1,
  3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1,
  4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1,
  9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1,
  11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1,
  11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1,
  2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1,
  9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1,
  3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1,
  1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1,
  4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1,
  4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1,
  0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1,
  3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1,
  3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1,
  0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1,
  9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1,
  1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

)foo") +

                                                                 std::string(R"foo(
__kernel void constructHPLevel(
                               __read_only image3d_t readHistoPyramid,
                               __write_only image3d_t writeHistoPyramid
                               ) {

  int4 writePos = {get_global_id(0), get_global_id(1), get_global_id(2), 0};
  int4 readPos = writePos*2;
  int writeValue = read_imagei(readHistoPyramid, sampler, readPos).x + // 0
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[1])).x + // 1
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[2])).x + // 2
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[3])).x + // 3
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[4])).x + // 4
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[5])).x + // 5
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[6])).x + // 6
  read_imagei(readHistoPyramid, sampler, (readPos+cubeOffsets[7])).x;  // 7

  write_imagei(writeHistoPyramid, writePos, writeValue);
}

int4 scanHPLevel(int target, __read_only image3d_t hp, int4 current) {

  int8 neighbors = {
    read_imagei(hp, sampler, current).x,
    read_imagei(hp, sampler, (current + cubeOffsets[1])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[2])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[3])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[4])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[5])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[6])).x,
    read_imagei(hp, sampler, (current + cubeOffsets[7])).x
  };

  int acc = current.s3 + neighbors.s0;
  int8 cmp;
  cmp.s0 = acc <= target;
  acc += neighbors.s1;
  cmp.s1 = acc <= target;
  acc += neighbors.s2;
  cmp.s2 = acc <= target;
  acc += neighbors.s3;
  cmp.s3 = acc <= target;
  acc += neighbors.s4;
  cmp.s4 = acc <= target;
  acc += neighbors.s5;
  cmp.s5 = acc <= target;
  acc += neighbors.s6;
  cmp.s6 = acc <= target;
  cmp.s7 = 0;

  current += cubeOffsets[(cmp.s0+cmp.s1+cmp.s2+cmp.s3+cmp.s4+cmp.s5+cmp.s6+cmp.s7)];
  current.s0 = current.s0*2;
  current.s1 = current.s1*2;
  current.s2 = current.s2*2;
  current.s3 = current.s3 +
  cmp.s0*neighbors.s0 +
  cmp.s1*neighbors.s1 +
  cmp.s2*neighbors.s2 +
  cmp.s3*neighbors.s3 +
  cmp.s4*neighbors.s4 +
  cmp.s5*neighbors.s5 +
  cmp.s6*neighbors.s6 +
  cmp.s7*neighbors.s7;
  return current;
}
)foo") +

                                                                 std::string(R"foo(
__kernel void traverseHP16(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // the field's gradient at either end of the edge, in the field's own sense (see the note at the head
    // of this source)
    const float4 centralDifference0 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 centralDifference1 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = centralDifference0 + (centralDifference1 - centralDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}
)foo") +

                                                                 std::string(R"foo(
__kernel void traverseHP32(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t hp4,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp4, cubePosition);
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // the field's gradient at either end of the edge, in the field's own sense (see the note at the head
    // of this source)
    const float4 centralDifference0 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 centralDifference1 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = centralDifference0 + (centralDifference1 - centralDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}
)foo") +

                                                                 std::string(R"foo(
__kernel void traverseHP64(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t hp4,
                         __read_only image3d_t hp5,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp5, cubePosition);
  cubePosition = scanHPLevel(target, hp4, cubePosition);
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // the field's gradient at either end of the edge, in the field's own sense (see the note at the head
    // of this source)
    const float4 centralDifference0 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 centralDifference1 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = centralDifference0 + (centralDifference1 - centralDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}

)foo") +

                                                                 std::string(R"foo(
__kernel void traverseHP128(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t hp4,
                         __read_only image3d_t hp5,
                         __read_only image3d_t hp6,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp6, cubePosition);
  cubePosition = scanHPLevel(target, hp5, cubePosition);
  cubePosition = scanHPLevel(target, hp4, cubePosition);
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // the field's gradient at either end of the edge, in the field's own sense (see the note at the head
    // of this source)
    const float4 centralDifference0 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 centralDifference1 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = centralDifference0 + (centralDifference1 - centralDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}
)foo") +

                                                                 std::string(R"foo(
__kernel void traverseHP256(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t hp4,
                         __read_only image3d_t hp5,
                         __read_only image3d_t hp6,
                         __read_only image3d_t hp7,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp7, cubePosition);
  cubePosition = scanHPLevel(target, hp6, cubePosition);
  cubePosition = scanHPLevel(target, hp5, cubePosition);
  cubePosition = scanHPLevel(target, hp4, cubePosition);
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // the field's gradient at either end of the edge, in the field's own sense (see the note at the head
    // of this source)
    const float4 centralDifference0 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 centralDifference1 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = centralDifference0 + (centralDifference1 - centralDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}
)foo") +

                                                                 std::string(R"foo(
__kernel void traverseHP512(
                         __read_only image3d_t hp0, // Largest HP
                         __read_only image3d_t hp1,
                         __read_only image3d_t hp2,
                         __read_only image3d_t hp3,
                         __read_only image3d_t hp4,
                         __read_only image3d_t hp5,
                         __read_only image3d_t hp6,
                         __read_only image3d_t hp7,
                         __read_only image3d_t hp8,
                         __read_only image3d_t rawData,
                         __global float * VBOBuffer,
                         __private int4 dimensions,
                         __private float isolevel,
                         __private int sum
                         ) {

  int target = get_global_id(0);
  if(target >= sum)
    target = 0;

  int4 cubePosition = {0,0,0,0}; // x,y,z,sum
  cubePosition = scanHPLevel(target, hp8, cubePosition);
  cubePosition = scanHPLevel(target, hp7, cubePosition);
  cubePosition = scanHPLevel(target, hp6, cubePosition);
  cubePosition = scanHPLevel(target, hp5, cubePosition);
  cubePosition = scanHPLevel(target, hp4, cubePosition);
  cubePosition = scanHPLevel(target, hp3, cubePosition);
  cubePosition = scanHPLevel(target, hp2, cubePosition);
  cubePosition = scanHPLevel(target, hp1, cubePosition);
  cubePosition = scanHPLevel(target, hp0, cubePosition);
  cubePosition.x = cubePosition.x / 2;
  cubePosition.y = cubePosition.y / 2;
  cubePosition.z = cubePosition.z / 2;

  char vertexNr = 0;
  const int4 cubeData = read_imagei(hp0, sampler, cubePosition);

  // max 5 triangles
  for(int i = (target-cubePosition.s3)*3; i < (target-cubePosition.s3+1)*3; i++)
  {
    // for each vertex in triangle
    const uchar edge = triTable[cubeData.y*16 + i];
    const int3 point0 = (int3)(cubePosition.x + offsets3[edge*6], cubePosition.y + offsets3[edge*6+1], cubePosition.z + offsets3[edge*6+2]);
    const int3 point1 = (int3)(cubePosition.x + offsets3[edge*6+3], cubePosition.y + offsets3[edge*6+4], cubePosition.z + offsets3[edge*6+5]);

    // the field's gradient at either end of the edge, in the field's own sense (see the note at the head
    // of this source)
    const float4 centralDifference0 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x+1, point0.y,   point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x-1, point0.y,   point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y+1, point0.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y-1, point0.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point0.x,   point0.y,   point0.z-1, 0) % dimensions).x),
                                               0.0f
                                               );
    const float4 centralDifference1 = (float4)(
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x+1, point1.y,   point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x-1, point1.y,   point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y+1, point1.z,   0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y-1, point1.z,   0) % dimensions).x),
                                               (float)(read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z+1, 0) % dimensions).x-
                                                       read_imagef(rawData, sampler, (int4)(point1.x,   point1.y,   point1.z-1, 0) % dimensions).x),
                                               0.0f
                                               );


    const float value0 = read_imagef(rawData, sampler, (int4)(point0.x, point0.y, point0.z, 0) % dimensions).x;
    const float diff = native_divide(
                                     (float)(isolevel-value0),
                                     (float)(read_imagef(rawData, sampler, (int4)(point1.x, point1.y, point1.z, 0) % dimensions).x - value0));

    const float4 vertex = (float4)(point0.x, point0.y, point0.z, 1.0f) + ((float4)(point1.x, point1.y, point1.z,0.0f) - (float4)(point0.x, point0.y, point0.z,0.0f)) * diff;
    const float4 scaledVertex = (float4)(vertex.x/(float)(dimensions.x),vertex.y/(float)(dimensions.y),vertex.z/(float)(dimensions.z),1.0f);

    const float4 normal = centralDifference0 + (centralDifference1 - centralDifference0) * diff;

    vstore4(scaledVertex, target*9 + vertexNr*3, VBOBuffer);
    vstore4(normal, target*9 + vertexNr*3 + 1, VBOBuffer);

    vertexNr++;
  }
}
)foo") +

                                                                 std::string(R"foo(
// The first part of the algorithm uses a table (edgeTable) which maps the vertices under the isosurface to the intersecting edges.
// An 8 bit index is formed where each bit corresponds to a vertex.
__kernel void classifyCubes(__write_only image3d_t histoPyramid,
                            __read_only image3d_t rawData,
                            __private int4 dimensions,
                            __private float isolevel)
{
  int4 pos = {get_global_id(0), get_global_id(1), get_global_id(2), 0};

  if(any(pos>=dimensions))
  {
    write_imageui(histoPyramid, pos, (uint4)(0, 0, 0, 0));
    return;
  }

  // Find cube class nr
  const float first = read_imagef(rawData, sampler, pos).x;
  const uchar cubeindex =
       ((first > isolevel)) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[1]) % dimensions).x > isolevel) << 1) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[3]) % dimensions).x > isolevel) << 2) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[2]) % dimensions).x > isolevel) << 3) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[4]) % dimensions).x > isolevel) << 4) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[5]) % dimensions).x > isolevel) << 5) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[7]) % dimensions).x > isolevel) << 6) |
       ((read_imagef(rawData, sampler, (pos + cubeOffsets[6]) % dimensions).x > isolevel) << 7);

  // Store number of triangles
  write_imageui(histoPyramid, pos, (uint4)(numberOfTriangles[cubeindex], cubeindex, first, 0));
}
)foo");
