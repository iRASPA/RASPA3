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

module opencl_pore_size_distribution;

import std;

import opencl;
import int3;
import uint3;
import double3;
import double3x3;
import simulationbox;
import framework;
import forcefield;
import units;
import opencl_clearance_grid;
import opencl_connected_components;


GridPoreSizeDistribution::GridPoreSizeDistribution() {}


GridPoreSizeDistribution::~GridPoreSizeDistribution() {}


void GridPoreSizeDistribution::run(const ForceField &forceField, const Framework &framework,
                                   std::string probePseudoAtom, uint3 gridSize,
                                   std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins)
{
  ClearanceGrid grid = ClearanceGrid::compute(forceField, framework, gridSize);
  this->run(forceField, framework, probePseudoAtom, grid, maximumDiameter, numberOfBins);
}


void GridPoreSizeDistribution::run(const ForceField &forceField, const Framework &framework,
                                   std::string probePseudoAtom, const ClearanceGrid &grid,
                                   std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins)
{
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value() || !OpenCL::clCommandQueue.has_value())
  {
    throw std::runtime_error("GridPoreSizeDistribution: no OpenCL device found\n");
  }

  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(
        std::format("GridPoreSizeDistribution: probe atom '{}' not found in the force field\n", probePseudoAtom));
  }
  this->probeRadius = 0.5 * forceField[probeType.value()].sizeParameter();

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  const std::int32_t nx = static_cast<std::int32_t>(grid.gridSize.x);
  const std::int32_t ny = static_cast<std::int32_t>(grid.gridSize.y);
  const std::int32_t nz = static_cast<std::int32_t>(grid.gridSize.z);
  const std::size_t numberOfVoxels = grid.numberOfVoxels();
  const double3x3 cell = framework.simulationBox.cell;

  // A sphere is drawn about every point of the void. None of them can be dropped for sitting inside another:
  // the field is a distance, so it gains exactly one along its own gradient and less in every other direction,
  // and one ball holds another only when the two centres and the nearest atom are collinear, which on a grid
  // never happens exactly.
  std::vector<cl_int> centres;
  centres.reserve(numberOfVoxels / 2);
  for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
  {
    if (grid.clearance[voxel] >= 0.0f) centres.push_back(static_cast<cl_int>(voxel));
  }
  const std::size_t voidVoxels = centres.size();
  this->numberOfCentres = centres.size();
  this->numberOfVoidVoxels = voidVoxels;

  // How far a sphere is allowed to reach past its own radius and still be said to cover a grid point. A grid
  // point stands for its voxel, and a sphere that reaches into the voxel has covered what the point stands for.
  // The size of it is set by the shell of void within a voxel of an atom's surface: the largest sphere covering
  // such a point touches it from the inside, so its centre has to sit on the line from the atom through the
  // point, and the nearest grid point to where the centre belongs is off that line by half a voxel diagonal,
  // which costs half a diagonal in the sphere's radius and another half in the distance to be covered.
  double3 halfStep = 0.5 * (cell * double3(1.0 / static_cast<double>(nx), 1.0 / static_cast<double>(ny),
                                           1.0 / static_cast<double>(nz)));
  this->slack = 2.0 * halfStep.length();

  std::vector<float> poreRadius(numberOfVoxels, 0.0f);

  if (!centres.empty())
  {
    std::chrono::steady_clock::time_point gpu_begin = std::chrono::steady_clock::now();

    cl_int err;
    const char *source = GridPoreSizeDistribution::poreSizeKernelSource;
    cl_program program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &source, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("GridPoreSizeDistribution: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
    }

    err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      std::size_t length;
      char buffer[2048];
      clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
      throw std::runtime_error(std::format("GridPoreSizeDistribution: OpenCL failed to build program (error: {})\n",
                                           std::string(buffer)));
    }

    cl_kernel kernel = clCreateKernel(program, "PoreSize", &err);
    if (err != CL_SUCCESS)
    {
      clReleaseProgram(program);
      throw std::runtime_error(std::format("GridPoreSizeDistribution: OpenCL clCreateKernel failed at {}\n", __LINE__));
    }

    std::vector<cl_float> clearance(numberOfVoxels);
    for (std::size_t v = 0; v < numberOfVoxels; ++v) clearance[v] = grid.clearance[v];

    cl_int centreError, clearanceError, poreError;
    cl_mem centreBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_int) * centres.size(),
                                         nullptr, &centreError);
    cl_mem clearanceBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                            sizeof(cl_float) * numberOfVoxels, nullptr, &clearanceError);
    cl_mem poreBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, sizeof(cl_float) * numberOfVoxels,
                                       nullptr, &poreError);
    if (centreError != CL_SUCCESS || clearanceError != CL_SUCCESS || poreError != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("GridPoreSizeDistribution: OpenCL clCreateBuffer failed at {}\n", __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), centreBuffer, CL_TRUE, 0, sizeof(cl_int) * centres.size(),
                               centres.data(), 0, nullptr, nullptr);
    err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), clearanceBuffer, CL_TRUE, 0,
                                sizeof(cl_float) * numberOfVoxels, clearance.data(), 0, nullptr, nullptr);
    err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), poreBuffer, CL_TRUE, 0,
                                sizeof(cl_float) * numberOfVoxels, poreRadius.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("GridPoreSizeDistribution: OpenCL clEnqueueWriteBuffer failed at {}\n", __LINE__));
    }

    cl_float4 cellRow0 = {{cl_float(cell[0][0]), cl_float(cell[1][0]), cl_float(cell[2][0]), 0.0f}};
    cl_float4 cellRow1 = {{cl_float(cell[0][1]), cl_float(cell[1][1]), cl_float(cell[2][1]), 0.0f}};
    cl_float4 cellRow2 = {{cl_float(cell[0][2]), cl_float(cell[1][2]), cl_float(cell[2][2]), 0.0f}};

    // How many grid steps along each axis a ball of a given radius can possibly reach. The perpendicular width
    // is the right measure for an oblique cell: a step along one axis moves that far towards the opposite face
    // at least, so the box built from it holds the ball.
    double3 widths = framework.simulationBox.perpendicularWidths();
    cl_float4 stepsPerLength = {{cl_float(static_cast<double>(nx) / widths.x), cl_float(static_cast<double>(ny) /
                                                                                       widths.y),
                                 cl_float(static_cast<double>(nz) / widths.z), 0.0f}};

    cl_int numberOfCentresArgument = static_cast<cl_int>(centres.size());
    cl_int3 clGridSize = {{nx, ny, nz, 0}};
    cl_float clSlack = static_cast<cl_float>(this->slack);

    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &centreBuffer);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &clearanceBuffer);
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
      throw std::runtime_error(std::format("GridPoreSizeDistribution: OpenCL clSetKernelArg failed at {}\n", __LINE__));
    }

    std::size_t globalWorkSize = centres.size();
    err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 1, nullptr, &globalWorkSize, nullptr, 0,
                                 nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("GridPoreSizeDistribution: OpenCL clEnqueueNDRangeKernel failed at {}\n", __LINE__));
    }

    err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), poreBuffer, CL_TRUE, 0, sizeof(cl_float) * numberOfVoxels,
                              poreRadius.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("GridPoreSizeDistribution: OpenCL clEnqueueReadBuffer failed at {}\n", __LINE__));
    }

    clFinish(OpenCL::clCommandQueue.value());

    clReleaseMemObject(centreBuffer);
    clReleaseMemObject(clearanceBuffer);
    clReleaseMemObject(poreBuffer);
    clReleaseKernel(kernel);
    clReleaseProgram(program);

    std::chrono::duration<double> gpuElapsed = std::chrono::steady_clock::now() - gpu_begin;
    this->gpuSeconds = gpuElapsed.count();
  }

  // Which of the void a probe can reach: the pores it fits in are the channels, and the void around them is
  // reached by anything that gets into the channel, so the void divides by which void component holds a
  // channel and which does not.
  GridComponents probeComponents = GridComponents::compute(grid, this->probeRadius);
  GridComponents voidComponents = GridComponents::compute(grid, 0.0);

  std::vector<std::uint8_t> voidComponentReached(voidComponents.pores.size(), 0);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    std::int32_t probePore = probeComponents.voxelPore[v];
    std::int32_t voidPore = voidComponents.voxelPore[v];
    if (probePore < 0 || voidPore < 0) continue;
    if (probeComponents.pores[static_cast<std::size_t>(probePore)].isChannel)
    {
      voidComponentReached[static_cast<std::size_t>(voidPore)] = 1;
    }
  }

  this->largestDiameter = 0.0;
  std::size_t reachedVoxels = 0;
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    if (grid.clearance[v] < 0.0f) continue;
    this->largestDiameter = std::max(this->largestDiameter, 2.0 * static_cast<double>(poreRadius[v]));
    std::int32_t voidPore = voidComponents.voxelPore[v];
    if (voidPore >= 0 && voidComponentReached[static_cast<std::size_t>(voidPore)] != 0) ++reachedVoxels;
  }

  const double totalVoxels = static_cast<double>(numberOfVoxels);
  this->voidFraction = static_cast<double>(voidVoxels) / totalVoxels;
  this->accessibleVoidFraction = static_cast<double>(reachedVoxels) / totalVoxels;

  // The curve is a count of grid points against the size found at each of them, so it needs no fitting: the
  // cumulative curve is exact for the field it was measured on, and only its derivative depends on how finely
  // the sizes are divided up.
  const std::size_t bins = numberOfBins.value_or(200);
  const double largest = maximumDiameter.value_or(std::max(this->largestDiameter, 1.0));
  const double step = largest / static_cast<double>(bins);

  std::vector<double> voidAtLeast(bins + 1, 0.0);
  std::vector<double> reachedAtLeast(bins + 1, 0.0);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    if (grid.clearance[v] < 0.0f) continue;
    double diameter = 2.0 * static_cast<double>(poreRadius[v]);
    std::size_t bin = static_cast<std::size_t>(std::floor(diameter / step));
    if (bin > bins) bin = bins;

    voidAtLeast[bin] += 1.0;
    std::int32_t voidPore = voidComponents.voxelPore[v];
    if (voidPore >= 0 && voidComponentReached[static_cast<std::size_t>(voidPore)] != 0) reachedAtLeast[bin] += 1.0;
  }

  // Turn the counts per size into counts at or above a size, from the top down.
  for (std::size_t bin = bins; bin-- > 0;)
  {
    voidAtLeast[bin] += voidAtLeast[bin + 1];
    reachedAtLeast[bin] += reachedAtLeast[bin + 1];
  }

  double voidNormalisation = (voidVoxels > 0) ? static_cast<double>(voidVoxels) : 1.0;
  double reachedNormalisation = (reachedVoxels > 0) ? static_cast<double>(reachedVoxels) : 1.0;

  this->points.resize(bins);
  for (std::size_t bin = 0; bin < bins; ++bin)
  {
    GridPoreSizeDistributionPoint point;
    point.diameter = (static_cast<double>(bin) + 0.5) * step;
    point.cumulativeVoidFraction = voidAtLeast[bin] / voidNormalisation;
    point.cumulativeAccessibleFraction = reachedAtLeast[bin] / reachedNormalisation;
    point.distribution = (voidAtLeast[bin] - voidAtLeast[bin + 1]) / (voidNormalisation * step);
    point.accessibleDistribution = (reachedAtLeast[bin] - reachedAtLeast[bin + 1]) / (reachedNormalisation * step);
    this->points[bin] = point;
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  this->seconds = elapsed.count();

  double3 spacing = grid.spacing();
  const double volume = framework.simulationBox.volume;

  std::ofstream report;
  report.open(framework.name + ".grid.psd.gpu.txt");
  std::print(report, "# Pore-size distribution (clearance grid, Gelb-Gubbins)\n");
  std::print(report, "# Framework: {}\n", framework.name);
  std::print(report, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(report, "# Number of framework atoms: {}\n", framework.unitCellAtoms.size());
  std::print(report, "# Framework volume: {} [Å³]\n", volume);
  std::print(report, "# Framework mass: {} [g/mol]\n", framework.unitCellMass);
  std::print(report, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", grid.gridSize.x,
             grid.gridSize.y, grid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(report, "# Grid points in the void: {} of {}, void fraction {}\n", voidVoxels, numberOfVoxels,
             this->voidFraction);
  std::print(report, "# Probe: {}, radius {} [Å], diameter {} [Å]\n", probePseudoAtom, this->probeRadius,
             2.0 * this->probeRadius);
  std::print(report, "# Of the void, the probe reaches {} [Å³], {} of it\n", this->accessibleVoidFraction * volume,
             (this->voidFraction > 0.0) ? this->accessibleVoidFraction / this->voidFraction : 0.0);
  std::print(report, "# Spheres drawn: {}, one about every point of the void\n", this->numberOfCentres);
  std::print(report, "# A sphere covers a point when it comes within {:.5f} [Å] of it, half a voxel diagonal twice\n",
             this->slack);
  std::print(report, "# The largest sphere that fits in the void is {:.5f} [Å] across\n", this->largestDiameter);
  std::print(report, "# Diameters evaluated: {}, up to {} [Å]\n", bins, largest);
  std::print(report, "# GPU Timing: {} [s] for the spheres, {} [s] for the clearance field\n", this->gpuSeconds,
             grid.seconds);
  std::print(report, "# CPU Timing: {} [s] in total\n", this->seconds);
  std::print(report, "#\n");
  std::print(report, "# The size at a point is the largest sphere that fits in the void and covers the point, so a\n");
  std::print(report, "# narrow neck between two wide cavities is given the width of a cavity rather than its own:\n");
  std::print(report, "# that is what the definition says, and it is why the curve has no weight at all below the\n");
  std::print(report, "# narrowest sphere the void can hold anywhere.\n");
  std::print(report, "#\n");
  std::print(report, "# The void nearest the atoms is the hard part of it. The largest sphere covering a point a\n");
  std::print(report, "# hair from an atom's surface touches that point from the inside, so its centre lies on the\n");
  std::print(report, "# line from the atom through the point and no grid point sits on that line. Letting a sphere\n");
  std::print(report, "# reach a voxel past itself is what covers that shell; without it the shell reports its own\n");
  std::print(report, "# clearance, which is nearly nothing, and the curve grows a spurious tail at small sizes.\n");
  std::print(report, "#\n");
  std::print(report, "# column 1: diameter d [Å]\n");
  std::print(report, "# column 2: share of the void in pores at least d across\n");
  std::print(report, "# column 3: the distribution, the derivative of column 2 [1/Å]\n");
  std::print(report, "# column 4: share of the void the probe reaches, in pores at least d across\n");
  std::print(report, "# column 5: the derivative of column 4 [1/Å]\n");
  std::print(report, "#      d [Å]     cumulative   distribution     accessible   distribution\n");
  for (const GridPoreSizeDistributionPoint &point : this->points)
  {
    std::print(report, "{:11.6f} {:14.8f} {:14.8f} {:14.8f} {:14.8f}\n", point.diameter, point.cumulativeVoidFraction,
               point.distribution, point.cumulativeAccessibleFraction, point.accessibleDistribution);
  }
  report.close();
}


const char *GridPoreSizeDistribution::poreSizeKernelSource = R"foo(
__kernel void PoreSize(__global const int *centre,
                       __global const float *clearance,
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
  float radius = clearance[voxel];
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
