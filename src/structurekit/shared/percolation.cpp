module;

module grid_percolation;

import std;

import int3;
import uint3;
import grid_connected_components;


// Union-find over grid points that carries, beside the region a point belongs to, how many cells away from
// its representative the point lies. Joining two points already in the same region then says something: if
// the two routes to them disagree about the cell, their difference is a translation that brings the region
// back to itself, and the region runs through the crystal rather than closing inside one cell.
namespace
{
struct PeriodicUnionFind
{
  std::vector<std::int32_t> parent;
  std::vector<int3> offset;
  std::vector<std::int32_t> componentSize;
  std::vector<float> highestOpenness;

  // Only the regions that have closed on themselves are in here, which is few of them.
  std::unordered_map<std::int32_t, std::vector<int3>> windings;

  // The members of each region that has not yet broken through, as a list threaded through the points
  // themselves so that two regions are joined by pointing the end of one at the start of the other. A region
  // that breaks through hands the level to everything on its list and then drops it, so no point is ever
  // written to twice and the lists between them are never longer than the points still trapped.
  bool recordEscape{false};
  std::vector<std::int32_t> nextMember;
  std::vector<std::int32_t> listHead;
  std::vector<std::int32_t> listTail;
  std::vector<std::uint8_t> brokenThrough;

  explicit PeriodicUnionFind(std::size_t n, bool recordEscape)
      : parent(n),
        offset(n, int3(0, 0, 0)),
        componentSize(n, 1),
        highestOpenness(n, 0.0f),
        recordEscape(recordEscape)
  {
    std::iota(this->parent.begin(), this->parent.end(), std::int32_t{0});
    if (recordEscape)
    {
      this->nextMember.assign(n, -1);
      this->listHead.assign(n, -1);
      this->listTail.assign(n, -1);
      this->brokenThrough.assign(n, 0);
    }
  }

  // Hands `level` to everything still trapped in `root` and marks the region as having broken through, after
  // which its points are answered as they arrive.
  void release(std::int32_t root, float level, std::span<const std::int32_t> order, std::span<float> escape)
  {
    if (!this->recordEscape) return;
    this->brokenThrough[static_cast<std::size_t>(root)] = 1;

    std::int32_t member = this->listHead[static_cast<std::size_t>(root)];
    while (member >= 0)
    {
      escape[static_cast<std::size_t>(order[static_cast<std::size_t>(member)])] = level;
      member = this->nextMember[static_cast<std::size_t>(member)];
    }
    this->listHead[static_cast<std::size_t>(root)] = -1;
    this->listTail[static_cast<std::size_t>(root)] = -1;
  }

  // Moves what is left of `other` onto `root`, the two having just become one region.
  void carryMembers(std::int32_t root, std::int32_t other)
  {
    if (!this->recordEscape) return;

    const std::int32_t head = this->listHead[static_cast<std::size_t>(other)];
    const std::int32_t tail = this->listTail[static_cast<std::size_t>(other)];
    this->listHead[static_cast<std::size_t>(other)] = -1;
    this->listTail[static_cast<std::size_t>(other)] = -1;
    if (head < 0) return;

    if (this->listHead[static_cast<std::size_t>(root)] < 0)
    {
      this->listHead[static_cast<std::size_t>(root)] = head;
    }
    else
    {
      this->nextMember[static_cast<std::size_t>(this->listTail[static_cast<std::size_t>(root)])] = head;
    }
    this->listTail[static_cast<std::size_t>(root)] = tail;
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


// Adds a translation to a region's collection, but only when it widens the pore system, so that a channel of
// a million points does not carry a million copies of the same direction. Returns the rank afterwards.
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


GridPercolation sweepPercolation(uint3 gridSize, std::span<const float> openness, float lowestOpenness,
                                 bool recordEscape)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  GridPercolation result;

  // No infinities: the build turns them off. The lowest finite double stands for a direction the pore
  // system never runs in, and survives the negation a caller reading energies will apply.
  const double never = std::numeric_limits<double>::lowest();
  result.opennessByDimension = {never, never, never};

  const std::int32_t nx = static_cast<std::int32_t>(gridSize.x);
  const std::int32_t ny = static_cast<std::int32_t>(gridSize.y);
  const std::int32_t nz = static_cast<std::int32_t>(gridSize.z);
  const std::size_t numberOfVoxels = openness.size();
  if (numberOfVoxels == 0) return result;

  result.highestOpenness = static_cast<double>(*std::ranges::max_element(openness));

  std::vector<std::int32_t> order;
  order.reserve(numberOfVoxels / 4);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    if (openness[v] >= lowestOpenness) order.push_back(static_cast<std::int32_t>(v));
  }
  result.numberOfVoxels = order.size();

  // A point that never joins a region running through the crystal has no way out at any level, and keeps
  // this. Points left out of the sweep for falling below the floor keep it too, no route having been offered
  // to them either.
  const float noWayOut = std::numeric_limits<float>::lowest();
  if (recordEscape) result.escapeOpenness.assign(numberOfVoxels, noWayOut);

  if (order.empty())
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
    result.seconds = elapsed.count();
    return result;
  }

  std::ranges::sort(order, [openness](std::int32_t a, std::int32_t b)
                    { return openness[static_cast<std::size_t>(a)] > openness[static_cast<std::size_t>(b)]; });

  // A point's position in the sweep is its label, so a neighbour is already switched on exactly when its
  // label is the smaller one, and no separate record of what is on is needed.
  std::vector<std::int32_t> label(numberOfVoxels, -1);
  for (std::size_t i = 0; i < order.size(); ++i)
  {
    label[static_cast<std::size_t>(order[i])] = static_cast<std::int32_t>(i);
  }

  const std::array<int3, 6> directions{int3(1, 0, 0),  int3(-1, 0, 0), int3(0, 1, 0),
                                       int3(0, -1, 0), int3(0, 0, 1),  int3(0, 0, -1)};

  PeriodicUnionFind unionFind(order.size(), recordEscape);

  for (std::size_t i = 0; i < order.size(); ++i)
  {
    const std::int32_t voxel = order[i];
    const float here = openness[static_cast<std::size_t>(voxel)];
    const std::int32_t current = static_cast<std::int32_t>(i);
    unionFind.highestOpenness[i] = here;
    if (recordEscape)
    {
      unionFind.listHead[i] = current;
      unionFind.listTail[i] = current;
    }

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

        // Whichever of the two already had a way out gives one to the other, here and now, at the level the
        // two of them met at. If neither has, what is trapped in the one joins what is trapped in the other
        // and they wait together.
        if (recordEscape)
        {
          const bool rootIsOut = unionFind.brokenThrough[static_cast<std::size_t>(root)] != 0;
          const bool otherIsOut = unionFind.brokenThrough[static_cast<std::size_t>(other)] != 0;
          if (rootIsOut && !otherIsOut)
          {
            unionFind.release(other, here, order, result.escapeOpenness);
          }
          else if (!rootIsOut && otherIsOut)
          {
            unionFind.release(root, here, order, result.escapeOpenness);
          }
          else if (!rootIsOut && !otherIsOut)
          {
            unionFind.carryMembers(root, other);
          }
        }

        unionFind.componentSize[static_cast<std::size_t>(root)] +=
            unionFind.componentSize[static_cast<std::size_t>(other)];
        unionFind.highestOpenness[static_cast<std::size_t>(root)] =
            std::max(unionFind.highestOpenness[static_cast<std::size_t>(root)],
                     unionFind.highestOpenness[static_cast<std::size_t>(other)]);

        // The translations of the region that is absorbed belong to the region that absorbs it.
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

      // A region that has just closed on itself runs through the crystal, so everything in it can now leave
      // by the level the closing happened at.
      if (recordEscape && rank >= 1 && unionFind.brokenThrough[static_cast<std::size_t>(root)] == 0)
      {
        unionFind.release(root, here, order, result.escapeOpenness);
      }

      // The openness at which the pore system first runs in one, two, and three directions. The sweep only
      // closes down, so whatever is seen first is the most open.
      for (int dimension = 1; dimension <= rank; ++dimension)
      {
        if (result.opennessByDimension[static_cast<std::size_t>(dimension - 1)] == never)
        {
          result.opennessByDimension[static_cast<std::size_t>(dimension - 1)] = static_cast<double>(here);
        }
      }

      if (!result.percolates && rank >= 1)
      {
        result.percolates = true;
        result.percolationOpenness = static_cast<double>(here);
        result.dimensionalityAtThreshold = rank;
        result.highestOpennessOnPath =
            static_cast<double>(unionFind.highestOpenness[static_cast<std::size_t>(root)]);
      }
    }
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  result.seconds = elapsed.count();

  return result;
}
