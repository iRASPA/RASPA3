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
  }
};

EnergyOpenCLSurfaceArea::~EnergyOpenCLSurfaceArea()
{
  if (OpenCL::clContext.has_value())
  {
    if (energyGridKernel != nullptr) clReleaseKernel(energyGridKernel);
    if (energyGridProgram != nullptr) clReleaseProgram(energyGridProgram);
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
