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
import structure_parallel;
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

  // The points of the void, and the ones among them a sphere is drawn about, which is all that have room of
  // some size: a point on a surface has none to draw with, and a point with no boundary in reach has no size
  // to give one. The sweep leaves those two out, so the count of spheres has to as well.
  std::size_t voidVoxels = 0;
  std::size_t centres = 0;
  for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
  {
    if (grid.clearance[voxel] >= 0.0f) ++voidVoxels;
    if (grid.clearance[voxel] > 0.0f && grid.clearance[voxel] != unlimitedRoom) ++centres;
  }
  this->numberOfCentres = centres;
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

  // Which of the void a probe can reach: the points its centre can be at in a channel, and the void its own
  // sphere covers from those. The void was divided before by which piece of the bare void held a channel, the
  // bare void being followed wherever the clearance is not below nought: a window narrower than the probe is
  // bare void all the same, and joined a cage sealed behind it to the channel in front, so that the whole of
  // the sodalite cages of LTA and FAU were counted as reached by a probe that cannot enter them.
  GridComponents probeComponents = GridComponents::compute(grid.gridSize, grid.clearance, this->probeRadius);

  std::vector<std::uint8_t> channelCentre(numberOfVoxels, 0);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    std::int32_t probePore = probeComponents.voxelPore[v];
    if (probePore >= 0 && probeComponents.pores[static_cast<std::size_t>(probePore)].isChannel) channelCentre[v] = 1;
  }

  std::vector<std::uint8_t> accessible(numberOfVoxels, 0);
  {
    const std::int64_t nx = grid.gridSize.x, ny = grid.gridSize.y, nz = grid.gridSize.z;
    // the grid's own slack for a centre that lies between its points
    const double reach = this->probeRadius + this->slack;
    const double3 widths = framework.unitCell.perpendicularWidths();
    const std::int64_t stepsX = static_cast<std::int64_t>(std::ceil(reach * static_cast<double>(nx) / widths.x));
    const std::int64_t stepsY = static_cast<std::int64_t>(std::ceil(reach * static_cast<double>(ny) / widths.y));
    const std::int64_t stepsZ = static_cast<std::int64_t>(std::ceil(reach * static_cast<double>(nz) / widths.z));
    const double3x3 cell = framework.unitCell.cell;

    forEachIndex(static_cast<std::size_t>(nz), workersAvailable(),
                 [&](std::size_t, std::size_t plane)
                 {
                   const std::int64_t k = static_cast<std::int64_t>(plane);
                   for (std::int64_t j = 0; j < ny; ++j)
                   {
                     for (std::int64_t i = 0; i < nx; ++i)
                     {
                       const std::size_t v = static_cast<std::size_t>((k * ny + j) * nx + i);
                       if (grid.clearance[v] < 0.0f) continue;
                       if (channelCentre[v] != 0)
                       {
                         accessible[v] = 1;
                         continue;
                       }
                       // a centre the probe fits at in a pocket is the pocket's, however near a channel
                       if (probeComponents.voxelPore[v] >= 0) continue;

                       bool covered = false;
                       for (std::int64_t dk = -stepsZ; dk <= stepsZ && !covered; ++dk)
                       {
                         for (std::int64_t dj = -stepsY; dj <= stepsY && !covered; ++dj)
                         {
                           for (std::int64_t di = -stepsX; di <= stepsX; ++di)
                           {
                             const std::size_t centre = static_cast<std::size_t>(
                                 ((((k + dk) % nz + nz) % nz) * ny + ((j + dj) % ny + ny) % ny) * nx +
                                 ((i + di) % nx + nx) % nx);
                             if (channelCentre[centre] == 0) continue;
                             const double3 dr =
                                 cell * double3(static_cast<double>(di) / static_cast<double>(nx),
                                                static_cast<double>(dj) / static_cast<double>(ny),
                                                static_cast<double>(dk) / static_cast<double>(nz));
                             if (double3::dot(dr, dr) <= reach * reach)
                             {
                               covered = true;
                               break;
                             }
                           }
                         }
                       }
                       if (covered) accessible[v] = 1;
                     }
                   }
                 });
  }

  // The curve is a count of grid points against the size found at each of them, binned by the one routine
  // every grid route bins with, so that where the top of the range falls is decided once.
  PoreSizeCurve curve = binPoreSizes(poreRadius, grid.clearance, accessible, {}, maximumDiameter, numberOfBins);
  this->largestDiameter = curve.largestDiameter;

  const double totalVoxels = static_cast<double>(numberOfVoxels);
  this->voidFraction = static_cast<double>(voidVoxels) / totalVoxels;
  this->accessibleVoidFraction = static_cast<double>(curve.numberOfAccessibleVoxels) / totalVoxels;

  this->points.resize(curve.points.size());
  for (std::size_t bin = 0; bin < curve.points.size(); ++bin)
  {
    const PoreSizeCurvePoint &from = curve.points[bin];
    this->points[bin] = GridPoreSizeDistributionPoint{from.diameter, from.cumulativeVoidFraction,
                                                      from.distribution, from.cumulativeAccessibleFraction,
                                                      from.accessibleDistribution};
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
  std::print(report, "# Spheres drawn: {}, one about every point of the void with room of some size\n",
             this->numberOfCentres);
  std::print(report, "# A sphere covers a point when it comes within {:.5f} [Å] of it, half a voxel diagonal twice\n",
             this->slack);
  std::print(report, "# The largest sphere that fits in the void is {} [Å] across\n",
             poreDiameterText(this->largestDiameter));
  std::print(report, "# Diameters evaluated: {}, up to {} [Å]\n", curve.points.size(),
             static_cast<double>(curve.points.size()) * curve.binWidth);
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

