module;

module opencl_pore_analysis;

import std;

import int3;
import uint3;
import double3;
import crystal;
import pair_interactions;
import opencl_clearance_grid;
import grid_connected_components;
import grid_percolation;


GridPoreDiameters GridPoreDiameters::compute(const ClearanceGrid &grid)
{
  GridPoreDiameters diameters;
  if (grid.numberOfVoxels() == 0) return diameters;

  // The clearance is an openness as it stands: a probe fits where it is large. Only the void takes part,
  // since a probe has a radius of at least nothing and a point inside an atom is of no use to it.
  GridPercolation swept = sweepPercolation(grid.gridSize, grid.clearance, 0.0f);

  // A sphere is twice its radius across, and the sweep dealt in radii.
  diameters.includedSphereDiameter = 2.0 * std::max(0.0, grid.maximumClearance());
  diameters.numberOfVoidVoxels = swept.numberOfVoxels;
  diameters.percolates = swept.percolates;
  diameters.dimensionalityAtThreshold = swept.dimensionalityAtThreshold;
  diameters.seconds = swept.seconds;

  if (swept.percolates)
  {
    diameters.freeSphereDiameter = 2.0 * swept.percolationOpenness;
    diameters.includedAlongFreePathDiameter = 2.0 * swept.highestOpennessOnPath;
  }

  // A direction the pore system never runs in is reported as a diameter of zero rather than as the sweep's
  // own symbol for it, which is what the reader below is written against.
  for (std::size_t dimension = 0; dimension < 3; ++dimension)
  {
    double openness = swept.opennessByDimension[dimension];
    diameters.freeSphereDiameterByDimension[dimension] =
        (openness == std::numeric_limits<double>::lowest()) ? 0.0 : 2.0 * openness;
  }

  return diameters;
}


GridPoreAnalysis::GridPoreAnalysis() {}


GridPoreAnalysis::~GridPoreAnalysis() {}


void GridPoreAnalysis::run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
                           uint3 gridSize)
{
  ClearanceGrid grid = ClearanceGrid::compute(interactions, framework, gridSize);
  this->run(interactions, framework, probePseudoAtom, grid);
}


void GridPoreAnalysis::run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
                           const ClearanceGrid &grid)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(std::format("GridPoreAnalysis: probe atom '{}' not found in the force field\n",
                                         probePseudoAtom));
  }
  this->probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  this->diameters = GridPoreDiameters::compute(grid);
  this->channels = GridComponents::compute(grid.gridSize, grid.clearance, this->probeRadius);

  double3 spacing = grid.spacing();

  std::ofstream diameterFile;
  diameterFile.open(framework.name + ".grid.res.gpu.txt");
  std::print(diameterFile, "# Pore diameters from the clearance grid (Di, Df, Dif)\n");
  std::print(diameterFile, "# Crystal: {}\n", framework.name);
  std::print(diameterFile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(diameterFile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(diameterFile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", grid.gridSize.x,
             grid.gridSize.y, grid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(diameterFile, "# Grid points in the void: {} of {}\n", this->diameters.numberOfVoidVoxels,
             grid.numberOfVoxels());
  std::print(diameterFile, "# GPU Timing: {} [s] for the clearance field\n", grid.seconds);
  std::print(diameterFile, "# CPU Timing: {} [s] for the sweep\n", this->diameters.seconds);
  std::print(diameterFile, "# A diameter is read off the grid, so it is right to within the grid spacing and is\n");
  std::print(diameterFile, "# reached from below: the widest point of a pore and the narrowest point of a\n");
  std::print(diameterFile, "# passage both fall between grid points, and a coarse grid understates them.\n");
  std::print(diameterFile, "Di (largest included sphere):            {} [Å]\n",
             this->diameters.includedSphereDiameter);
  std::print(diameterFile, "Df (largest free sphere):                {} [Å]\n", this->diameters.freeSphereDiameter);
  std::print(diameterFile, "Dif (included sphere along free path):   {} [Å]\n",
             this->diameters.includedAlongFreePathDiameter);
  if (!this->diameters.percolates)
  {
    std::print(diameterFile, "The structure does not percolate: every pore closes inside one cell.\n");
  }
  std::print(diameterFile, "\n# The widest free sphere that still runs in a given number of directions. The sweep\n");
  std::print(diameterFile, "# passes each of these on its way down, so they come from the same pass as Df.\n");
  for (int dimension = 1; dimension <= 3; ++dimension)
  {
    double diameter = this->diameters.freeSphereDiameterByDimension[static_cast<std::size_t>(dimension - 1)];
    if (diameter > 0.0)
      std::print(diameterFile, "Df in {} direction(s) or more:            {} [Å]\n", dimension, diameter);
    else
      std::print(diameterFile, "Df in {} direction(s) or more:            none\n", dimension);
  }
  diameterFile.close();

  std::ofstream channelFile;
  channelFile.open(framework.name + ".grid.chan.gpu.txt");
  std::print(channelFile, "# Channel and pocket analysis from the clearance grid\n");
  std::print(channelFile, "# Crystal: {}\n", framework.name);
  std::print(channelFile, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, this->probeRadius);
  std::print(channelFile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", grid.gridSize.x,
             grid.gridSize.y, grid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(channelFile, "# Grid points the probe centre fits at: {} of {}\n", this->channels.numberOfOpenVoxels,
             grid.numberOfVoxels());
  std::print(channelFile, "# CPU Timing: {} [s]\n", this->channels.seconds);
  std::print(channelFile, "# A pore's dimensionality is the number of independent lattice translations that bring\n");
  std::print(channelFile, "# it back to itself: none for a pocket, and one, two or three for a channel, a layer,\n");
  std::print(channelFile, "# and a network. The translations are integers, so the count is decided rather than\n");
  std::print(channelFile, "# estimated, and the grid enters only through which points the probe fits at.\n");
  std::print(channelFile, "Number of channels: {}\n", this->channels.numberOfChannels);
  std::print(channelFile, "Number of pockets:  {}\n", this->channels.numberOfPockets);
  std::print(channelFile, "  of those, pockets holding a single grid point: {}\n",
             this->channels.numberOfSinglePointPockets);
  std::print(channelFile, "Pore system dimensionality: {}\n", this->channels.dimensionality);
  for (std::size_t i = 0; i < this->channels.pores.size(); ++i)
  {
    const GridPore &pore = this->channels.pores[i];
    std::print(channelFile, "  pore {}: {} dimensionality={} points={} widest={:.5f} [Å]\n", i,
               pore.isChannel ? "channel" : "pocket", pore.dimensionality, pore.numberOfVoxels,
               2.0 * pore.largestOpenness);
  }
  channelFile.close();
}
