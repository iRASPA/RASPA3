module;

module opencl_pore_analysis;

import std;

import int3;
import uint3;
import double3;
import framework;
import forcefield;
import opencl_clearance_grid;
import opencl_connected_components;


// Union-find over grid points that carries, beside the pore each point belongs to, how many cells away from
// its representative the point lies. Joining two points already in the same pore then says something: if the
// two routes to them disagree about the cell, their difference is a translation that brings the pore back to
// itself, and the pore runs through the crystal rather than closing inside one cell.
namespace
{
struct PeriodicUnionFind
{
  std::vector<std::int32_t> parent;
  std::vector<int3> offset;
  std::vector<std::int32_t> componentSize;
  std::vector<float> largestClearance;

  // Only the pores that have closed on themselves are in here, which is few of them.
  std::unordered_map<std::int32_t, std::vector<int3>> windings;

  explicit PeriodicUnionFind(std::size_t n)
      : parent(n), offset(n, int3(0, 0, 0)), componentSize(n, 1), largestClearance(n, 0.0f)
  {
    std::iota(this->parent.begin(), this->parent.end(), std::int32_t{0});
  }

  // The representative of `x`, and how far `x` lies from it in cells.
  std::int32_t find(std::int32_t x, int3 &displacement)
  {
    std::int32_t root = x;
    int3 accumulated(0, 0, 0);
    while (this->parent[static_cast<std::size_t>(root)] != root)
    {
      accumulated += this->offset[static_cast<std::size_t>(root)];
      root = this->parent[static_cast<std::size_t>(root)];
    }

    // Point everything along the way straight at the representative, with the displacement it needs.
    std::int32_t current = x;
    int3 remaining = accumulated;
    while (this->parent[static_cast<std::size_t>(current)] != current)
    {
      std::int32_t next = this->parent[static_cast<std::size_t>(current)];
      int3 rest = remaining - this->offset[static_cast<std::size_t>(current)];
      this->parent[static_cast<std::size_t>(current)] = root;
      this->offset[static_cast<std::size_t>(current)] = remaining;
      current = next;
      remaining = rest;
    }

    displacement = accumulated;
    return root;
  }
};


// Adds a translation to a pore's collection, but only when it widens the pore system, so that a channel of a
// million points does not carry a million copies of the same direction. Returns the rank afterwards.
int addWinding(std::vector<int3> &basis, const int3 &winding)
{
  int before = latticeVectorRank(basis);
  if (before == 3) return 3;

  std::vector<int3> trial = basis;
  trial.push_back(winding);
  int after = latticeVectorRank(trial);
  if (after > before) basis.push_back(winding);

  return after;
}
}  // namespace


GridPoreDiameters GridPoreDiameters::compute(const ClearanceGrid &grid)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  GridPoreDiameters diameters;

  const std::int32_t nx = static_cast<std::int32_t>(grid.gridSize.x);
  const std::int32_t ny = static_cast<std::int32_t>(grid.gridSize.y);
  const std::int32_t nz = static_cast<std::int32_t>(grid.gridSize.z);
  const std::size_t numberOfVoxels = grid.numberOfVoxels();
  if (numberOfVoxels == 0) return diameters;

  diameters.includedSphereDiameter = 2.0 * std::max(0.0, grid.maximumClearance());

  // Only the void takes part: a probe has a radius of at least nothing, so a point inside an atom is of no
  // use to it however the sweep is ordered.
  std::vector<std::int32_t> order;
  order.reserve(numberOfVoxels / 4);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    if (grid.clearance[v] >= 0.0f) order.push_back(static_cast<std::int32_t>(v));
  }
  diameters.numberOfVoidVoxels = order.size();
  if (order.empty())
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
    diameters.seconds = elapsed.count();
    return diameters;
  }

  std::ranges::sort(order, [&grid](std::int32_t a, std::int32_t b)
                    { return grid.clearance[static_cast<std::size_t>(a)] >
                             grid.clearance[static_cast<std::size_t>(b)]; });

  // A point's position in the sweep is its label, so a neighbour is already switched on exactly when its
  // label is the smaller one, and no separate record of what is on is needed.
  std::vector<std::int32_t> label(numberOfVoxels, -1);
  for (std::size_t i = 0; i < order.size(); ++i)
  {
    label[static_cast<std::size_t>(order[i])] = static_cast<std::int32_t>(i);
  }

  const std::array<int3, 6> directions{int3(1, 0, 0),  int3(-1, 0, 0), int3(0, 1, 0),
                                       int3(0, -1, 0), int3(0, 0, 1),  int3(0, 0, -1)};

  PeriodicUnionFind unionFind(order.size());

  for (std::size_t i = 0; i < order.size(); ++i)
  {
    const std::int32_t voxel = order[i];
    const float room = grid.clearance[static_cast<std::size_t>(voxel)];
    const std::int32_t current = static_cast<std::int32_t>(i);
    unionFind.largestClearance[i] = room;

    const std::int32_t k = voxel / (nx * ny);
    const std::int32_t j = (voxel / nx) % ny;
    const std::int32_t iv = voxel % nx;

    for (const int3 &direction : directions)
    {
      std::int32_t ni = iv + direction.x;
      std::int32_t nj = j + direction.y;
      std::int32_t nk = k + direction.z;
      int3 crossed(0, 0, 0);

      if (ni < 0)
      {
        ni += nx;
        crossed.x = -1;
      }
      else if (ni >= nx)
      {
        ni -= nx;
        crossed.x = 1;
      }
      if (nj < 0)
      {
        nj += ny;
        crossed.y = -1;
      }
      else if (nj >= ny)
      {
        nj -= ny;
        crossed.y = 1;
      }
      if (nk < 0)
      {
        nk += nz;
        crossed.z = -1;
      }
      else if (nk >= nz)
      {
        nk -= nz;
        crossed.z = 1;
      }

      const std::size_t neighbourVoxel = static_cast<std::size_t>((nk * ny + nj) * nx + ni);
      const std::int32_t neighbour = label[neighbourVoxel];
      if (neighbour < 0 || neighbour >= current) continue;

      int3 displacementHere(0, 0, 0);
      int3 displacementThere(0, 0, 0);
      std::int32_t rootHere = unionFind.find(current, displacementHere);
      std::int32_t rootThere = unionFind.find(neighbour, displacementThere);

      // Where the neighbour sits relative to this point, taken round by way of the representatives.
      int3 shift = crossed + displacementHere - displacementThere;

      std::int32_t root = rootHere;
      int rank = 0;
      if (rootHere == rootThere)
      {
        if (shift.x == 0 && shift.y == 0 && shift.z == 0) continue;
        rank = addWinding(unionFind.windings[rootHere], shift);
      }
      else
      {
        if (unionFind.componentSize[static_cast<std::size_t>(rootHere)] >=
            unionFind.componentSize[static_cast<std::size_t>(rootThere)])
        {
          unionFind.parent[static_cast<std::size_t>(rootThere)] = rootHere;
          unionFind.offset[static_cast<std::size_t>(rootThere)] = shift;
          root = rootHere;
        }
        else
        {
          unionFind.parent[static_cast<std::size_t>(rootHere)] = rootThere;
          unionFind.offset[static_cast<std::size_t>(rootHere)] = int3(0, 0, 0) - shift;
          root = rootThere;
        }

        std::int32_t other = (root == rootHere) ? rootThere : rootHere;
        unionFind.componentSize[static_cast<std::size_t>(root)] +=
            unionFind.componentSize[static_cast<std::size_t>(other)];
        unionFind.largestClearance[static_cast<std::size_t>(root)] =
            std::max(unionFind.largestClearance[static_cast<std::size_t>(root)],
                     unionFind.largestClearance[static_cast<std::size_t>(other)]);

        // The translations of the pore that is absorbed belong to the pore that absorbs it.
        std::unordered_map<std::int32_t, std::vector<int3>>::iterator absorbed = unionFind.windings.find(other);
        if (absorbed != unionFind.windings.end())
        {
          std::vector<int3> carried = absorbed->second;
          unionFind.windings.erase(absorbed);
          for (const int3 &winding : carried) rank = addWinding(unionFind.windings[root], winding);
        }
        std::unordered_map<std::int32_t, std::vector<int3>>::iterator kept = unionFind.windings.find(root);
        if (kept != unionFind.windings.end()) rank = latticeVectorRank(kept->second);
      }

      // The first probe size at which the pore system runs in one, two, and three directions. The sweep only
      // narrows, so whatever is seen first is the widest.
      for (int dimension = 1; dimension <= rank; ++dimension)
      {
        if (diameters.freeSphereDiameterByDimension[static_cast<std::size_t>(dimension - 1)] == 0.0)
        {
          diameters.freeSphereDiameterByDimension[static_cast<std::size_t>(dimension - 1)] =
              2.0 * static_cast<double>(room);
        }
      }

      if (!diameters.percolates && rank >= 1)
      {
        diameters.percolates = true;
        diameters.freeSphereDiameter = 2.0 * static_cast<double>(room);
        diameters.dimensionalityAtThreshold = rank;
        diameters.includedAlongFreePathDiameter =
            2.0 * static_cast<double>(unionFind.largestClearance[static_cast<std::size_t>(root)]);
      }
    }
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  diameters.seconds = elapsed.count();

  return diameters;
}


GridPoreAnalysis::GridPoreAnalysis() {}


GridPoreAnalysis::~GridPoreAnalysis() {}


void GridPoreAnalysis::run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom,
                           uint3 gridSize)
{
  ClearanceGrid grid = ClearanceGrid::compute(forceField, framework, gridSize);
  this->run(forceField, framework, probePseudoAtom, grid);
}


void GridPoreAnalysis::run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom,
                           const ClearanceGrid &grid)
{
  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(std::format("GridPoreAnalysis: probe atom '{}' not found in the force field\n",
                                         probePseudoAtom));
  }
  this->probeRadius = 0.5 * forceField[probeType.value()].sizeParameter();

  this->diameters = GridPoreDiameters::compute(grid);
  this->channels = GridComponents::compute(grid, this->probeRadius);

  double3 spacing = grid.spacing();

  std::ofstream diameterFile;
  diameterFile.open(framework.name + ".grid.res.gpu.txt");
  std::print(diameterFile, "# Pore diameters from the clearance grid (Di, Df, Dif)\n");
  std::print(diameterFile, "# Framework: {}\n", framework.name);
  std::print(diameterFile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(diameterFile, "# Number of framework atoms: {}\n", framework.unitCellAtoms.size());
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
  std::print(channelFile, "# Framework: {}\n", framework.name);
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
               2.0 * pore.largestClearance);
  }
  channelFile.close();
}
