module;

module opencl_void_fraction;

import std;

import uint3;
import double3;
import crystal;
import pair_interactions;
import units;
import opencl_clearance_grid;
import grid_connected_components;


GridVoidFraction::GridVoidFraction() {}


GridVoidFraction::~GridVoidFraction() {}


void GridVoidFraction::run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
                           uint3 gridSize)
{
  ClearanceGrid grid = ClearanceGrid::compute(interactions, framework, gridSize);
  this->run(interactions, framework, probePseudoAtom, grid);
}


void GridVoidFraction::run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
                           const ClearanceGrid &grid)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(
        std::format("GridVoidFraction: probe atom '{}' not found in the force field\n", probePseudoAtom));
  }
  this->probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  GridComponents components = GridComponents::compute(grid.gridSize, grid.clearance, this->probeRadius);

  const double volume = framework.unitCell.volume;
  const double numberOfVoxels = static_cast<double>(grid.numberOfVoxels());

  std::size_t channelVoxels = 0;
  std::size_t pocketVoxels = 0;
  for (const GridPore &pore : components.pores)
  {
    if (pore.isChannel)
      channelVoxels += pore.numberOfVoxels;
    else
      pocketVoxels += pore.numberOfVoxels;
  }

  this->numberOfChannels = components.numberOfChannels;
  this->numberOfPockets = components.numberOfPockets;

  this->voidFraction = static_cast<double>(components.numberOfOpenVoxels) / numberOfVoxels;
  this->accessibleVolumeFraction = static_cast<double>(channelVoxels) / numberOfVoxels;
  this->inaccessibleVolumeFraction = static_cast<double>(pocketVoxels) / numberOfVoxels;

  this->voidVolume = this->voidFraction * volume;
  this->accessibleVolume = this->accessibleVolumeFraction * volume;
  this->inaccessibleVolume = this->inaccessibleVolumeFraction * volume;

  double densityCrystal = framework.mass / (volume * Units::Angstrom * Units::Angstrom * Units::Angstrom *
                                                    1.0e6 * Units::AvogadroConstant);
  double toGravimetric = (Units::Angstrom * Units::Angstrom * Units::Angstrom * 1.0e6) * Units::AvogadroConstant /
                         framework.mass;

  double3 spacing = grid.spacing();

  std::ofstream myfile;
  myfile.open(framework.name + ".grid.av.gpu.txt");
  std::print(myfile, "# Void volume (clearance grid)\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, this->probeRadius);
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", grid.gridSize.x,
             grid.gridSize.y, grid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Grid points the probe centre fits at: {} of {}\n", components.numberOfOpenVoxels,
             grid.numberOfVoxels());
  std::print(myfile, "# Channels: {}, pockets: {}, of those {} holding a single grid point\n",
             components.numberOfChannels, components.numberOfPockets, components.numberOfSinglePointPockets);
  std::print(myfile, "# Pore system dimensionality: {}\n", components.dimensionality);
  std::print(myfile, "# Crystal volume: {} [Å³]\n", volume);
  std::print(myfile, "# Crystal density: {} [g/cm³]\n", densityCrystal);
  std::print(myfile, "# GPU Timing: {} [s] for the clearance field\n", grid.seconds);
  std::print(myfile, "# CPU Timing: {} [s] for the division into channels and pockets\n", components.seconds);
  std::print(myfile, "# The three volumes are counts of grid points, so the void is exact to the point and the\n");
  std::print(myfile, "# division of it is exact to the neck: a passage narrower than the grid spacing is closed\n");
  std::print(myfile, "# on a coarse grid and open on a fine one, which moves volume from a pocket to a channel\n");
  std::print(myfile, "# rather than changing the total.\n");
  std::print(myfile, "Void volume:         fraction {}  {} [Å³]  {} [cm³/g]\n", this->voidFraction, this->voidVolume,
             this->voidVolume * toGravimetric);
  std::print(myfile, "Accessible volume:   fraction {}  {} [Å³]  {} [cm³/g]\n", this->accessibleVolumeFraction,
             this->accessibleVolume, this->accessibleVolume * toGravimetric);
  std::print(myfile, "Inaccessible volume: fraction {}  {} [Å³]  {} [cm³/g]\n", this->inaccessibleVolumeFraction,
             this->inaccessibleVolume, this->inaccessibleVolume * toGravimetric);

  if (components.numberOfPockets > 0)
  {
    std::print(myfile, "\n");
    std::print(myfile, "# Each pocket, from the grid points in it. The centre is their mean, taken over the pocket\n");
    std::print(myfile, "# followed out of the cell, so a pocket lying across a cell boundary gets its own centre.\n");
    std::print(myfile, "# `widest` is the largest ball that fits anywhere inside the pocket and `equivalent` is\n");
    std::print(myfile, "# the radius of the ball of the pocket's own volume: they agree for a round cavity and\n");
    std::print(myfile, "# spread apart for one drawn out along a passage.\n");
    std::print(myfile, "#     s_a         s_b         s_c      volume [Å³]  points   widest [Å]  equivalent [Å]\n");

    double voxelVolume = grid.voxelVolume();
    for (const GridPore &pore : components.pores)
    {
      if (pore.isChannel) continue;

      double pocketVolume = static_cast<double>(pore.numberOfVoxels) * voxelVolume;
      double equivalentRadius = std::cbrt(3.0 * pocketVolume / (4.0 * std::numbers::pi));
      std::print(myfile, "Pocket: {:11.7f} {:11.7f} {:11.7f} {:11.5f} {:8} {:11.5f} {:11.5f}\n",
                 pore.centroidFractional.x, pore.centroidFractional.y, pore.centroidFractional.z, pocketVolume,
                 pore.numberOfVoxels, pore.largestOpenness, equivalentRadius);
    }
  }

  myfile.close();
}
