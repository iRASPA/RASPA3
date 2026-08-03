module;

module brute_force_diameters;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import brute_force_structure;
import brute_force_voxels;

namespace
{
// Walk uphill from a point until the step is smaller than anything that could matter. At each turn the
// twenty-six directions around the current point are tried at the current step; the best of them is taken
// if it finds more room, and the step is halved if none does. Deterministic, so that two runs of the check
// agree with each other as well as with what they are checking.
std::pair<double3, double> walkUphill(const BruteForceStructure &structure, double3 from, double startingStep)
{
  double3 best = from;
  double bestClearance = structure.clearance(best);

  double step = startingStep;

  while (step > 1.0e-10)
  {
    bool moved = false;

    for (std::int32_t k = -1; k <= 1; ++k)
    {
      for (std::int32_t j = -1; j <= 1; ++j)
      {
        for (std::int32_t i = -1; i <= 1; ++i)
        {
          if (i == 0 && j == 0 && k == 0) continue;

          double3 direction(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
          double3 trial = best + (step / direction.length()) * direction;

          double clearance = structure.clearance(trial);
          if (clearance > bestClearance)
          {
            bestClearance = clearance;
            best = trial;
            moved = true;
          }
        }
      }
    }

    if (!moved) step *= 0.5;
  }

  return {best, bestClearance};
}

// Union-find over voxels carrying the lattice translation to the root, so that a second route back to a
// voxel already reached shows up as the vector the piece repeats under.
struct PeriodicUnionFind
{
  std::vector<std::size_t> parent;
  std::vector<int3> offset;

  explicit PeriodicUnionFind(std::size_t count) : parent(count), offset(count, int3(0, 0, 0))
  {
    std::iota(parent.begin(), parent.end(), std::size_t{0});
  }

  std::pair<std::size_t, int3> find(std::size_t node)
  {
    std::size_t root = node;
    int3 total(0, 0, 0);
    while (parent[root] != root)
    {
      total = int3(total.x + offset[root].x, total.y + offset[root].y, total.z + offset[root].z);
      root = parent[root];
    }

    // Flatten, so that the next walk is short.
    std::size_t walk = node;
    int3 remaining = total;
    while (parent[walk] != walk)
    {
      std::size_t next = parent[walk];
      int3 step = offset[walk];

      parent[walk] = root;
      offset[walk] = remaining;

      remaining = int3(remaining.x - step.x, remaining.y - step.y, remaining.z - step.z);
      walk = next;
    }

    return {root, total};
  }

  // Joins two voxels; returns the lattice vector the piece is found to repeat under, when there is one.
  std::optional<int3> join(std::size_t a, std::size_t b, const int3 &crossing)
  {
    auto [rootA, toRootA] = find(a);
    auto [rootB, toRootB] = find(b);

    int3 between(toRootA.x + crossing.x - toRootB.x, toRootA.y + crossing.y - toRootB.y,
                 toRootA.z + crossing.z - toRootB.z);

    if (rootA == rootB)
    {
      if (between.x != 0 || between.y != 0 || between.z != 0) return between;
      return std::nullopt;
    }

    parent[rootB] = rootA;
    offset[rootB] = between;
    return std::nullopt;
  }
};

std::size_t latticeRank(const std::vector<int3> &vectors)
{
  std::array<std::array<std::int64_t, 3>, 3> rows{};
  std::size_t rank = 0;

  for (const int3 &vector : vectors)
  {
    std::array<std::int64_t, 3> row{vector.x, vector.y, vector.z};

    for (std::size_t pivot = 0; pivot < rank; ++pivot)
    {
      std::size_t column = 0;
      while (column < 3 && rows[pivot][column] == 0) ++column;
      if (column == 3 || row[column] == 0) continue;

      std::int64_t a = rows[pivot][column];
      std::int64_t b = row[column];
      for (std::size_t c = 0; c < 3; ++c) row[c] = row[c] * a - rows[pivot][c] * b;
    }

    if (row[0] != 0 || row[1] != 0 || row[2] != 0)
    {
      rows[rank] = row;
      ++rank;
      if (rank == 3) break;
    }
  }

  return rank;
}

// One hop between neighbouring voxels, and how wide a sphere it takes.
struct Hop
{
  std::size_t from;
  std::size_t to;
  int3 crossing;
  float width;
};
}  // namespace

BruteForceDiameters BruteForceDiameters::compute(const BruteForceStructure &structure,
                                                 const BruteForceVoxels &voxels)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  BruteForceDiameters self;
  self.spacing = voxels.spacing;
  self.numberOfVoidVoxels = voxels.numberOfVoidVoxels;

  std::size_t numberOfVoxels = voxels.clearance.size();
  if (numberOfVoxels == 0) return self;

  // Di: the roomiest few voxels, each walked uphill. Starting from more than one guards against the best
  // voxel sitting in a basin that is not the deepest.
  {
    std::vector<std::size_t> ranked;
    ranked.reserve(voxels.numberOfVoidVoxels);
    for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
    {
      if (voxels.regionOf[voxel] >= 0) ranked.push_back(voxel);
    }

    const std::size_t starts = std::min<std::size_t>(64, ranked.size());
    std::partial_sort(ranked.begin(), ranked.begin() + static_cast<std::ptrdiff_t>(starts), ranked.end(),
                      [&](std::size_t a, std::size_t b) { return voxels.clearance[a] > voxels.clearance[b]; });

    double bestOnTheGrid = ranked.empty() ? 0.0 : static_cast<double>(voxels.clearance[ranked.front()]);

    double bestClearance = 0.0;
    double3 bestPosition{};

    double step = std::max({voxels.spacing.x, voxels.spacing.y, voxels.spacing.z});

#pragma omp parallel for schedule(dynamic, 1)
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(starts); ++index)
    {
      auto [position, clearance] =
          walkUphill(structure, voxels.centre(structure, ranked[static_cast<std::size_t>(index)]), step);

#pragma omp critical
      {
        if (clearance > bestClearance)
        {
          bestClearance = clearance;
          bestPosition = position;
        }
      }
    }

    self.includedSphereDiameter = 2.0 * bestClearance;
    self.includedSphereCentre = bestPosition;
    self.walkGainedForDi = bestClearance - bestOnTheGrid;
  }

  // Df: every pair of neighbouring void voxels is a hop, and the hops are taken widest first until the
  // piece being built arrives back at a voxel it already holds by a route that crossed the boundary. The
  // width at which that happens is how wide a sphere can get all the way through the crystal.
  //
  // The width of a hop is the room at its narrower end. The room along the hop can be less than that, but
  // by no more than half its length, the clearance changing by no more than the distance moved; so the same
  // pass run with that subtracted gives a width every sphere is shown to manage, and the two together
  // bracket the answer.
  std::vector<Hop> hops;
  hops.reserve(voxels.numberOfVoidVoxels * 13);

  {
    const std::array<int3, 13> steps{int3(1, 0, 0),   int3(0, 1, 0),   int3(0, 0, 1),  int3(1, 1, 0),
                                     int3(1, -1, 0),  int3(1, 0, 1),   int3(1, 0, -1), int3(0, 1, 1),
                                     int3(0, 1, -1),  int3(1, 1, 1),   int3(1, 1, -1), int3(1, -1, 1),
                                     int3(1, -1, -1)};

    for (std::int32_t k = 0; k < voxels.counts.z; ++k)
    {
      for (std::int32_t j = 0; j < voxels.counts.y; ++j)
      {
        for (std::int32_t i = 0; i < voxels.counts.x; ++i)
        {
          std::size_t here = voxels.indexOf(i, j, k);
          if (voxels.regionOf[here] < 0) continue;

          for (const int3 &step : steps)
          {
            std::int32_t ni = i + step.x;
            std::int32_t nj = j + step.y;
            std::int32_t nk = k + step.z;

            int3 crossing(0, 0, 0);
            auto wrap = [](std::int32_t &value, std::int32_t count, std::int32_t &carry)
            {
              if (value >= count)
              {
                value -= count;
                carry = 1;
              }
              else if (value < 0)
              {
                value += count;
                carry = -1;
              }
            };
            wrap(ni, voxels.counts.x, crossing.x);
            wrap(nj, voxels.counts.y, crossing.y);
            wrap(nk, voxels.counts.z, crossing.z);

            std::size_t there = voxels.indexOf(ni, nj, nk);
            if (voxels.regionOf[there] < 0) continue;

            hops.push_back(Hop{.from = here,
                               .to = there,
                               .crossing = crossing,
                               .width = std::min(voxels.clearance[here], voxels.clearance[there])});
          }
        }
      }
    }
  }
  self.numberOfEdges = hops.size();

  std::ranges::sort(hops, [](const Hop &a, const Hop &b) { return a.width > b.width; });

  // The pass, run at whatever width each hop is credited with. Returns the width at which the void first
  // runs away through the crystal, and the vectors it runs away along.
  auto widestWayThrough = [&](double lengthPenalty)
  {
    PeriodicUnionFind pieces(numberOfVoxels);

    std::vector<int3> found;
    double at = 0.0;
    std::size_t closedAt = hops.size();

    for (std::size_t index = 0; index < hops.size(); ++index)
    {
      const Hop &hop = hops[index];

      if (std::optional<int3> repeat = pieces.join(hop.from, hop.to, hop.crossing); repeat.has_value())
      {
        if (found.empty())
        {
          at = static_cast<double>(hop.width) - lengthPenalty;
          closedAt = index;
        }
        found.push_back(repeat.value());
      }

      // Once it runs away, keep going only long enough to see in how many directions it does.
      if (!found.empty() && index > closedAt + hops.size() / 20 && latticeRank(found) == 3) break;
    }

    return std::pair<double, std::vector<int3>>{at, found};
  };

  // The hop lengths differ, so the penalty that makes every width one a sphere is shown to manage is half
  // the longest of them. Applying the worst case to all of them keeps the bracket honest and the pass
  // single.
  double longestHop = std::sqrt(voxels.spacing.x * voxels.spacing.x + voxels.spacing.y * voxels.spacing.y +
                                voxels.spacing.z * voxels.spacing.z);

  auto [width, repeats] = widestWayThrough(0.0);

  self.percolates = !repeats.empty();
  self.dimensionality = latticeRank(repeats);
  self.freeSphereDiameter = 2.0 * width;
  self.guaranteedFreeSphereDiameter = 2.0 * std::max(width - 0.5 * longestHop, 0.0);

  // Dif: the roomiest place that can be reached without ever passing anywhere narrower than the bottleneck,
  // walked uphill from the best voxel of that component.
  if (self.percolates)
  {
    PeriodicUnionFind pieces(numberOfVoxels);
    for (const Hop &hop : hops)
    {
      if (static_cast<double>(hop.width) < width) break;
      pieces.join(hop.from, hop.to, hop.crossing);
    }

    // Which piece the hop that closed the loop belongs to is the percolating one.
    std::size_t seed = hops.front().from;
    for (const Hop &hop : hops)
    {
      if (static_cast<double>(hop.width) < width) break;
      seed = hop.from;
    }
    std::size_t percolatingRoot = pieces.find(seed).first;

    std::size_t best = seed;
    double bestClearance = -std::numeric_limits<double>::max();
    for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
    {
      if (voxels.regionOf[voxel] < 0) continue;
      if (static_cast<double>(voxels.clearance[voxel]) < width) continue;
      if (pieces.find(voxel).first != percolatingRoot) continue;

      if (static_cast<double>(voxels.clearance[voxel]) > bestClearance)
      {
        bestClearance = voxels.clearance[voxel];
        best = voxel;
      }
    }

    double step = std::max({voxels.spacing.x, voxels.spacing.y, voxels.spacing.z});
    auto [position, clearance] = walkUphill(structure, voxels.centre(structure, best), step);

    self.includedAlongFreePathDiameter = 2.0 * clearance;
    self.includedAlongFreePathCentre = position;
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - begin;
  self.seconds = timing.count();

  return self;
}
