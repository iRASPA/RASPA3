module;

module opencl_connected_components;

import std;

import int3;
import double3;
import opencl_clearance_grid;


int latticeVectorRank(const std::vector<int3> &vectors)
{
  using vector = std::array<std::int64_t, 3>;
  auto cross = [](const vector &u, const vector &v)
  { return vector{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]}; };
  auto nonZero = [](const vector &u) { return u[0] != 0 || u[1] != 0 || u[2] != 0; };

  vector along{}, normal{};
  int rank = 0;
  for (const int3 &v : vectors)
  {
    vector row{v.x, v.y, v.z};
    if (!nonZero(row)) continue;

    if (rank == 0)
    {
      along = row;
      rank = 1;
    }
    else if (rank == 1)
    {
      vector off = cross(along, row);
      if (nonZero(off))
      {
        normal = off;
        rank = 2;
      }
    }
    else if (normal[0] * row[0] + normal[1] * row[1] + normal[2] * row[2] != 0)
    {
      return 3;
    }
  }
  return rank;
}


GridComponents::GridComponents() {}


GridComponents::~GridComponents() {}


double GridComponents::occupiedFraction() const
{
  if (this->voxelPore.empty()) return 0.0;
  return static_cast<double>(this->numberOfOpenVoxels) / static_cast<double>(this->voxelPore.size());
}


double GridComponents::channelFraction() const
{
  if (this->voxelPore.empty()) return 0.0;

  std::size_t channelVoxels = 0;
  for (const GridPore &pore : this->pores)
  {
    if (pore.isChannel) channelVoxels += pore.numberOfVoxels;
  }
  return static_cast<double>(channelVoxels) / static_cast<double>(this->voxelPore.size());
}


GridComponents GridComponents::compute(const ClearanceGrid &grid, double probeRadius)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  GridComponents components;
  components.probeRadius = probeRadius;

  const std::int32_t nx = static_cast<std::int32_t>(grid.gridSize.x);
  const std::int32_t ny = static_cast<std::int32_t>(grid.gridSize.y);
  const std::int32_t nz = static_cast<std::int32_t>(grid.gridSize.z);
  const std::size_t numberOfVoxels = grid.numberOfVoxels();

  components.voxelPore.assign(numberOfVoxels, -1);
  if (numberOfVoxels == 0) return components;

  const float threshold = static_cast<float>(probeRadius);

  // Two open points touching only at an edge or a corner share a passage of no cross-section, so the void is
  // followed through faces alone. That is also what keeps a structure from being called percolating and
  // sealed at once: the void takes the six face neighbours and the solid gets the diagonals.
  const std::array<int3, 6> directions{int3(1, 0, 0),  int3(-1, 0, 0), int3(0, 1, 0),
                                       int3(0, -1, 0), int3(0, 0, 1),  int3(0, 0, -1)};

  // Only the points the probe fits at are labelled, and the walk below works on those labels, so the
  // bookkeeping is proportional to the open space rather than to the cell.
  std::vector<std::int32_t> label(numberOfVoxels, -1);
  std::vector<std::int32_t> labelVoxel;
  labelVoxel.reserve(numberOfVoxels / 4);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    if (grid.clearance[v] >= threshold)
    {
      label[v] = static_cast<std::int32_t>(labelVoxel.size());
      labelVoxel.push_back(static_cast<std::int32_t>(v));
    }
  }
  components.numberOfOpenVoxels = labelVoxel.size();
  if (labelVoxel.empty())
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
    components.seconds = elapsed.count();
    return components;
  }

  std::vector<int3> displacement(labelVoxel.size(), int3(0, 0, 0));
  std::vector<std::uint8_t> visited(labelVoxel.size(), 0);
  std::vector<std::int32_t> stack;

  for (std::size_t start = 0; start < labelVoxel.size(); ++start)
  {
    if (visited[start]) continue;

    GridPore pore;
    std::vector<int3> windings;
    std::size_t poreVoxels = 0;
    double3 centreSum(0.0, 0.0, 0.0);

    stack.clear();
    stack.push_back(static_cast<std::int32_t>(start));
    visited[start] = 1;
    displacement[start] = int3(0, 0, 0);

    std::int32_t poreId = static_cast<std::int32_t>(components.pores.size());

    {
      std::int32_t v = labelVoxel[start];
      std::int32_t k = v / (nx * ny);
      std::int32_t j = (v / nx) % ny;
      std::int32_t i = v % nx;
      pore.exampleVoxel = int3(i, j, k);
      pore.widestVoxel = pore.exampleVoxel;
      pore.largestClearance = static_cast<double>(grid.clearance[static_cast<std::size_t>(v)]);
    }

    while (!stack.empty())
    {
      const std::size_t current = static_cast<std::size_t>(stack.back());
      stack.pop_back();

      const std::int32_t voxel = labelVoxel[current];
      components.voxelPore[static_cast<std::size_t>(voxel)] = poreId;
      ++poreVoxels;

      const double room = static_cast<double>(grid.clearance[static_cast<std::size_t>(voxel)]);
      const std::int32_t k = voxel / (nx * ny);
      const std::int32_t j = (voxel / nx) % ny;
      const std::int32_t i = voxel % nx;
      if (room > pore.largestClearance)
      {
        pore.largestClearance = room;
        pore.widestVoxel = int3(i, j, k);
      }

      const int3 here = displacement[current];

      centreSum += double3((static_cast<double>(i) / static_cast<double>(nx)) + static_cast<double>(here.x),
                           (static_cast<double>(j) / static_cast<double>(ny)) + static_cast<double>(here.y),
                           (static_cast<double>(k) / static_cast<double>(nz)) + static_cast<double>(here.z));

      for (const int3 &direction : directions)
      {
        // Where the step leaves the cell it comes back in on the far side, one cell over, and that cell is
        // the translation the step crosses. Adding those up along a walk is what makes a loop back to the
        // starting point report the translation it returned by rather than nothing at all.
        std::int32_t ni = i + direction.x;
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
        if (label[neighbourVoxel] < 0) continue;
        const std::size_t neighbour = static_cast<std::size_t>(label[neighbourVoxel]);

        const int3 there = here + crossed;
        if (!visited[neighbour])
        {
          visited[neighbour] = 1;
          displacement[neighbour] = there;
          stack.push_back(static_cast<std::int32_t>(neighbour));
        }
        else
        {
          int3 loop = there - displacement[neighbour];
          if (loop.x != 0 || loop.y != 0 || loop.z != 0)
          {
            // Only translations that widen the pore system are kept, so a channel of a million points does
            // not carry a million copies of the same direction.
            if (windings.size() < 3)
            {
              std::vector<int3> trial = windings;
              trial.push_back(loop);
              if (latticeVectorRank(trial) > latticeVectorRank(windings)) windings.push_back(loop);
            }
          }
        }
      }
    }

    pore.numberOfVoxels = poreVoxels;
    if (poreVoxels > 0)
    {
      double3 centre = centreSum / static_cast<double>(poreVoxels);
      pore.centroidFractional = double3(centre.x - std::floor(centre.x), centre.y - std::floor(centre.y),
                                        centre.z - std::floor(centre.z));
    }
    pore.windingVectors = windings;
    pore.dimensionality = latticeVectorRank(windings);
    pore.isChannel = pore.dimensionality > 0;

    if (pore.isChannel)
      ++components.numberOfChannels;
    else
      ++components.numberOfPockets;

    if (!pore.isChannel && pore.numberOfVoxels == 1) ++components.numberOfSinglePointPockets;

    components.dimensionality = std::max(components.dimensionality, pore.dimensionality);
    components.pores.push_back(pore);
  }

  // Channels first, then the roomiest pores, so that what a report shows first is what holds the most.
  std::vector<std::int32_t> ordering(components.pores.size());
  std::iota(ordering.begin(), ordering.end(), std::int32_t{0});
  std::ranges::sort(ordering, [&components](std::int32_t a, std::int32_t b)
                    {
                      const GridPore &left = components.pores[static_cast<std::size_t>(a)];
                      const GridPore &right = components.pores[static_cast<std::size_t>(b)];
                      if (left.isChannel != right.isChannel) return left.isChannel;
                      if (left.dimensionality != right.dimensionality) return left.dimensionality > right.dimensionality;
                      return left.numberOfVoxels > right.numberOfVoxels;
                    });

  std::vector<std::int32_t> renumbered(components.pores.size());
  std::vector<GridPore> sorted(components.pores.size());
  for (std::size_t position = 0; position < ordering.size(); ++position)
  {
    renumbered[static_cast<std::size_t>(ordering[position])] = static_cast<std::int32_t>(position);
    sorted[position] = components.pores[static_cast<std::size_t>(ordering[position])];
  }
  components.pores = sorted;
  for (std::int32_t &pore : components.voxelPore)
  {
    if (pore >= 0) pore = renumbered[static_cast<std::size_t>(pore)];
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  components.seconds = elapsed.count();

  return components;
}
