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
import unit_cell;
import crystal;
import pair_interactions;
import units;
import opencl_clearance_grid;
import grid_connected_components;
import grid_pore_size;
import grid_pore_size_opencl;


GridPoreSizeDistribution::GridPoreSizeDistribution() {}


GridPoreSizeDistribution::~GridPoreSizeDistribution() {}


void GridPoreSizeDistribution::run(const PairInteractions &interactions, const Crystal &framework,
                                   std::string probePseudoAtom, uint3 gridSize,
                                   std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins)
{
  ClearanceGrid grid = ClearanceGrid::compute(interactions, framework, gridSize);
  this->run(interactions, framework, probePseudoAtom, grid, maximumDiameter, numberOfBins);
}


void GridPoreSizeDistribution::run(const PairInteractions &interactions, const Crystal &framework,
                                   std::string probePseudoAtom, const ClearanceGrid &grid,
                                   std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins)
{
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value() || !OpenCL::clCommandQueue.has_value())
  {
    throw std::runtime_error("GridPoreSizeDistribution: no OpenCL device found\n");
  }

  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(
        std::format("GridPoreSizeDistribution: probe atom '{}' not found in the force field\n", probePseudoAtom));
  }
  this->probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  const std::size_t numberOfVoxels = grid.numberOfVoxels();

  // A sphere is drawn about every point of the void.
  std::size_t voidVoxels = 0;
  for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
  {
    if (grid.clearance[voxel] >= 0.0f) ++voidVoxels;
  }
  this->numberOfCentres = voidVoxels;
  this->numberOfVoidVoxels = voidVoxels;

  this->slack = coveringSlack(framework.unitCell, grid.gridSize);

  std::vector<float> poreRadius(numberOfVoxels, 0.0f);

  if (voidVoxels > 0)
  {
    std::chrono::steady_clock::time_point gpu_begin = std::chrono::steady_clock::now();
    poreRadius = poreRadiusFieldOpenCL(grid.gridSize, framework.unitCell, grid.clearance, this->slack);
    std::chrono::duration<double> gpuElapsed = std::chrono::steady_clock::now() - gpu_begin;
    this->gpuSeconds = gpuElapsed.count();
  }

  // Which of the void a probe can reach: the pores it fits in are the channels, and the void around them is
  // reached by anything that gets into the channel, so the void divides by which void component holds a
  // channel and which does not.
  GridComponents probeComponents = GridComponents::compute(grid.gridSize, grid.clearance, this->probeRadius);
  GridComponents voidComponents = GridComponents::compute(grid.gridSize, grid.clearance, 0.0);

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
  const double volume = framework.unitCell.volume;

  std::ofstream report;
  report.open(framework.name + ".grid.psd.gpu.txt");
  std::print(report, "# Pore-size distribution (clearance grid, Gelb-Gubbins)\n");
  std::print(report, "# Crystal: {}\n", framework.name);
  std::print(report, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(report, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(report, "# Crystal volume: {} [Å³]\n", volume);
  std::print(report, "# Crystal mass: {} [g/mol]\n", framework.mass);
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

