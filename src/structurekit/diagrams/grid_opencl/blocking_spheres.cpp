module;

module opencl_blocking_spheres;

import std;

import int3;
import uint3;
import double3;
import double3x3;
import crystal;
import pair_interactions;
import unit_cell;
import opencl_clearance_grid;
import grid_connected_components;
import grid_blocking_cover;


GridBlockingSpheres::GridBlockingSpheres() {}


GridBlockingSpheres::~GridBlockingSpheres() {}


void GridBlockingSpheres::run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
                              uint3 gridSize)
{
  ClearanceGrid grid = ClearanceGrid::compute(interactions, framework, gridSize);
  this->run(interactions, framework, probePseudoAtom, grid);
}


void GridBlockingSpheres::run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
                              const ClearanceGrid &grid)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(
        std::format("GridBlockingSpheres: probe atom '{}' not found in the force field\n", probePseudoAtom));
  }
  this->probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;
  this->gridSeconds = grid.seconds;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  GridComponents components = GridComponents::compute(grid.gridSize, grid.clearance, this->probeRadius);

  const std::size_t numberOfVoxels = grid.numberOfVoxels();

  this->numberOfPockets = components.numberOfPockets;
  this->numberOfSinglePointPockets = components.numberOfSinglePointPockets;

  // Where the probe both fits and can get out, which is what a sphere must not reach into, and which of the
  // pores have to be covered: on this route the two are the same division, a pore being sealed exactly when
  // it does not run anywhere.
  std::vector<std::uint8_t> isChannel(numberOfVoxels, 0);
  for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
  {
    std::int32_t pore = components.voxelPore[voxel];
    if (pore < 0) continue;
    if (components.pores[static_cast<std::size_t>(pore)].isChannel) isChannel[voxel] = 1;
  }

  std::vector<std::uint8_t> needsCover(components.pores.size(), 0);
  for (std::size_t poreIndex = 0; poreIndex < components.pores.size(); ++poreIndex)
  {
    needsCover[poreIndex] = components.pores[poreIndex].isChannel ? 0 : 1;
  }

  GridBlockingCover cover = coverPockets(grid.gridSize, framework.unitCell, components.voxelPore, components.pores,
                                     needsCover, grid.clearance, isChannel);

  this->spheres = cover.spheres;
  this->numberOfClippedSpheres = cover.numberOfClippedSpheres;
  this->numberOfRefusedPoints = cover.numberOfRefusedPoints;
  this->numberOfPocketVoxels = cover.numberOfPocketVoxels;
  this->pocketFraction = static_cast<double>(cover.numberOfPocketVoxels) / static_cast<double>(numberOfVoxels);

  const double voxelVolume = grid.voxelVolume();
  const std::size_t pocketVoxels = cover.numberOfPocketVoxels;

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  this->seconds = elapsed.count();

  // What a simulation reads: nothing but the numbers, a comment being more than the format allows.
  std::ofstream blockFile;
  blockFile.open(framework.name + ".block");
  std::print(blockFile, "{}\n", this->spheres.size());
  for (const GridBlockingSphere &sphere : this->spheres)
  {
    std::print(blockFile, "{} {} {} {}\n", sphere.centreFractional.x, sphere.centreFractional.y,
               sphere.centreFractional.z, sphere.radius);
  }
  blockFile.close();

  double3 spacing = grid.spacing();

  std::ofstream report;
  report.open(framework.name + ".grid.block.gpu.txt");
  std::print(report, "# Blocking spheres (clearance grid)\n");
  std::print(report, "# Crystal: {}\n", framework.name);
  std::print(report, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, this->probeRadius);
  std::print(report, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", grid.gridSize.x,
             grid.gridSize.y, grid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(report, "# Channels: {}, pockets: {}, of those {} holding a single grid point\n",
             components.numberOfChannels, this->numberOfPockets, this->numberOfSinglePointPockets);
  std::print(report, "# Grid points in the pockets: {} of {}, {} [Å³], a fraction {} of the cell\n", pocketVoxels,
             numberOfVoxels, static_cast<double>(pocketVoxels) * voxelVolume, this->pocketFraction);
  std::print(report, "# Spheres written to {}.block: {}\n", framework.name, this->spheres.size());
  std::print(report, "# Of those, {} were stopped by a channel rather than by the pocket's own extent\n",
             this->numberOfClippedSpheres);
  std::print(report, "# {} points of a pocket had no sphere placed on them, lying nearer a channel than the grid\n",
             this->numberOfRefusedPoints);
  std::print(report, "# can resolve: a neck the grid has closed rather than a pocket, holding no volume\n");
  std::print(report, "# GPU Timing: {} [s] for the clearance field\n", this->gridSeconds);
  std::print(report, "# CPU Timing: {} [s] for the pockets and the spheres\n", this->seconds);
  std::print(report, "#\n");
  std::print(report, "# A pocket is a piece of the region where the probe fits that does not run anywhere, so what\n");
  std::print(report, "# is blocked here is decided by the shape of that region and not by sampling it. Each pocket\n");
  std::print(report, "# is covered from its widest point outwards, and a sphere stops at the nearest point the\n");
  std::print(report, "# probe can actually reach: a sphere past that would block pore volume the probe is owed,\n");
  std::print(report, "# and a simulation given such a file would report an uptake too low for a reason nothing in\n");
  std::print(report, "# its output would explain.\n");
  std::print(report, "#\n");
  std::print(report, "# The grid decides where the necks are. A passage narrower than the spacing reads as closed,\n");
  std::print(report, "# which turns a piece of channel into a pocket and blocks it; a finer grid opens it again.\n");
  std::print(report, "# So the count of pockets is the part of this to check against a finer grid, the spheres\n");
  std::print(report, "# themselves being no more than a covering of whatever the pockets turned out to be.\n");
  std::print(report, "#\n");
  std::print(report, "# column 1-3: centre, fractional\n");
  std::print(report, "# column 4: radius [Å]\n");
  std::print(report, "# column 5: which pocket\n");
  std::print(report, "# column 6: grid points of that pocket this sphere covered\n");
  std::print(report, "# column 7: what stopped it growing\n");
  std::print(report, "#       s_a         s_b         s_c   radius [Å]  pocket    points  stopped by\n");
  for (const GridBlockingSphere &sphere : this->spheres)
  {
    std::print(report, "Sphere: {:11.7f} {:11.7f} {:11.7f} {:11.5f} {:7} {:9}  {}\n", sphere.centreFractional.x,
               sphere.centreFractional.y, sphere.centreFractional.z, sphere.radius, sphere.pocket,
               sphere.voxelsCovered,
               sphere.clipped ? "a channel" : "the pocket");
  }
  report.close();
}
