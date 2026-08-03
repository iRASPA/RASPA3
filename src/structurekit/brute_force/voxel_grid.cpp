module;

module brute_force_voxels;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import brute_force_structure;

namespace
{
// Union-find over the void voxels, carrying with each voxel the lattice translation between it and the root
// of its piece. Joining two voxels already in the same piece by a route with a different translation finds
// a vector the piece repeats under, which is what percolation is.
struct PeriodicUnionFind
{
  std::vector<std::size_t> parent;
  std::vector<int3> offset;  // from a voxel to its parent, in cells
  std::vector<std::vector<int3>> cycles;  // per root, the translations found

  explicit PeriodicUnionFind(std::size_t count) : parent(count), offset(count, int3(0, 0, 0)), cycles(count)
  {
    std::iota(parent.begin(), parent.end(), std::size_t{0});
  }

  // The root of a voxel's piece and the translation to it, flattening the tree as it goes.
  std::pair<std::size_t, int3> find(std::size_t node)
  {
    if (parent[node] == node) return {node, int3(0, 0, 0)};

    auto [root, toRoot] = find(parent[node]);
    int3 total(offset[node].x + toRoot.x, offset[node].y + toRoot.y, offset[node].z + toRoot.z);

    parent[node] = root;
    offset[node] = total;
    return {root, total};
  }

  void join(std::size_t a, std::size_t b, const int3 &crossing)
  {
    auto [rootA, toRootA] = find(a);
    auto [rootB, toRootB] = find(b);

    // Going a -> b costs `crossing`; the same trip through the roots costs -toRootA + (translation between
    // the roots) + toRootB, which is what fixes the translation when the roots differ.
    int3 between(toRootA.x + crossing.x - toRootB.x, toRootA.y + crossing.y - toRootB.y,
                 toRootA.z + crossing.z - toRootB.z);

    if (rootA == rootB)
    {
      // Two routes to the same place. Where they disagree, the difference repeats the piece.
      if (between.x != 0 || between.y != 0 || between.z != 0) cycles[rootA].push_back(between);
      return;
    }

    parent[rootB] = rootA;
    offset[rootB] = between;

    std::vector<int3> &into = cycles[rootA];
    into.insert(into.end(), cycles[rootB].begin(), cycles[rootB].end());
    cycles[rootB].clear();
    cycles[rootB].shrink_to_fit();
  }
};

// How many independent directions a set of integer vectors spans, by elimination over the integers so that
// no tolerance has to be chosen.
std::size_t latticeRank(std::vector<int3> vectors)
{
  std::array<std::array<std::int64_t, 3>, 3> rows{};
  std::size_t rank = 0;

  for (const int3 &vector : vectors)
  {
    std::array<std::int64_t, 3> row{vector.x, vector.y, vector.z};

    for (std::size_t pivot = 0; pivot < rank; ++pivot)
    {
      // Find the column this pivot row leads in and clear it, working in integers by cross-multiplying.
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
}  // namespace

double3 BruteForceVoxels::fractionalCentre(std::size_t index) const
{
  std::size_t nx = static_cast<std::size_t>(this->counts.x);
  std::size_t ny = static_cast<std::size_t>(this->counts.y);

  std::size_t i = index % nx;
  std::size_t j = (index / nx) % ny;
  std::size_t k = index / (nx * ny);

  return double3((static_cast<double>(i) + 0.5) / static_cast<double>(this->counts.x),
                 (static_cast<double>(j) + 0.5) / static_cast<double>(this->counts.y),
                 (static_cast<double>(k) + 0.5) / static_cast<double>(this->counts.z));
}

double3 BruteForceVoxels::centre(const BruteForceStructure &structure, std::size_t index) const
{
  return structure.unitCell.cell * this->fractionalCentre(index);
}

std::size_t BruteForceVoxels::voxelOf(const BruteForceStructure &structure, const double3 &position) const
{
  double3 s = structure.wrappedFractional(position);

  auto axis = [](double coordinate, std::int32_t count)
  {
    std::int32_t index = static_cast<std::int32_t>(coordinate * static_cast<double>(count));
    return std::clamp(index, std::int32_t{0}, count - 1);
  };

  return this->indexOf(axis(s.x, this->counts.x), axis(s.y, this->counts.y), axis(s.z, this->counts.z));
}

std::int32_t BruteForceVoxels::regionNear(const BruteForceStructure &structure, const double3 &position,
                                          std::int32_t rings) const
{
  std::size_t here = this->voxelOf(structure, position);
  if (this->regionOf[here] >= 0) return this->regionOf[here];

  std::size_t nx = static_cast<std::size_t>(this->counts.x);
  std::size_t ny = static_cast<std::size_t>(this->counts.y);

  std::int32_t i0 = static_cast<std::int32_t>(here % nx);
  std::int32_t j0 = static_cast<std::int32_t>((here / nx) % ny);
  std::int32_t k0 = static_cast<std::int32_t>(here / (nx * ny));

  double3 wanted = structure.wrappedFractional(position);

  std::int32_t best = -1;
  double nearest = std::numeric_limits<double>::max();

  for (std::int32_t dk = -rings; dk <= rings; ++dk)
  {
    for (std::int32_t dj = -rings; dj <= rings; ++dj)
    {
      for (std::int32_t di = -rings; di <= rings; ++di)
      {
        std::int32_t i = ((i0 + di) % this->counts.x + this->counts.x) % this->counts.x;
        std::int32_t j = ((j0 + dj) % this->counts.y + this->counts.y) % this->counts.y;
        std::int32_t k = ((k0 + dk) % this->counts.z + this->counts.z) % this->counts.z;

        std::size_t voxel = this->indexOf(i, j, k);
        std::int32_t region = this->regionOf[voxel];
        if (region < 0) continue;

        double3 to = structure.nearestImage(structure.unitCell.cell * wanted,
                                            this->centre(structure, voxel));
        double distance = double3::dot(to, to);
        if (distance < nearest)
        {
          nearest = distance;
          best = region;
        }
      }
    }
  }

  return best;
}

bool BruteForceVoxels::isInVoid(const BruteForceStructure &structure, const double3 &position) const
{
  return this->regionOf[this->voxelOf(structure, position)] >= 0;
}

bool BruteForceVoxels::isAccessible(const BruteForceStructure &structure, const double3 &position) const
{
  std::int32_t region = this->regionOf[this->voxelOf(structure, position)];
  return region >= 0 && this->regions[static_cast<std::size_t>(region)].percolates;
}

BruteForceVoxels BruteForceVoxels::build(const BruteForceStructure &structure, double targetSpacing,
                                         double threshold)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  BruteForceVoxels self;

  const UnitCell &box = structure.unitCell;
  double3 widths = box.perpendicularWidths();

  auto along = [&](double width)
  {
    return static_cast<std::int32_t>(
        std::clamp(std::round(width / std::max(targetSpacing, 1.0e-6)), 4.0, 1024.0));
  };
  self.counts = int3(along(widths.x), along(widths.y), along(widths.z));

  self.spacing = double3(widths.x / static_cast<double>(self.counts.x), widths.y / static_cast<double>(self.counts.y),
                         widths.z / static_cast<double>(self.counts.z));

  std::size_t numberOfVoxels = static_cast<std::size_t>(self.counts.x) * static_cast<std::size_t>(self.counts.y) *
                               static_cast<std::size_t>(self.counts.z);
  self.voxelVolume = box.volume / static_cast<double>(numberOfVoxels);

  // The clearance at every voxel centre, the long way round.
  self.clearance.assign(numberOfVoxels, 0.0f);

#pragma omp parallel for schedule(static)
  for (std::int64_t index = 0; index < static_cast<std::int64_t>(numberOfVoxels); ++index)
  {
    std::size_t voxel = static_cast<std::size_t>(index);
    self.clearance[voxel] = static_cast<float>(structure.clearance(self.centre(structure, voxel)));
  }

  // Flood the voxels with room in them, joining each to the three neighbours ahead of it so that every pair
  // is looked at once, and noting which joins crossed a cell boundary.
  self.regionOf.assign(numberOfVoxels, -1);

  std::vector<std::size_t> voidVoxels;
  for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
  {
    if (self.clearance[voxel] >= static_cast<float>(threshold)) voidVoxels.push_back(voxel);
  }
  self.numberOfVoidVoxels = voidVoxels.size();

  PeriodicUnionFind pieces(numberOfVoxels);

  const std::array<int3, 3> steps{int3(1, 0, 0), int3(0, 1, 0), int3(0, 0, 1)};

  for (std::int32_t k = 0; k < self.counts.z; ++k)
  {
    for (std::int32_t j = 0; j < self.counts.y; ++j)
    {
      for (std::int32_t i = 0; i < self.counts.x; ++i)
      {
        std::size_t here = self.indexOf(i, j, k);
        if (self.clearance[here] < static_cast<float>(threshold)) continue;

        for (const int3 &step : steps)
        {
          std::int32_t ni = i + step.x;
          std::int32_t nj = j + step.y;
          std::int32_t nk = k + step.z;

          // Which way the step left the cell, which is the translation the join carries.
          int3 crossing(0, 0, 0);
          if (ni >= self.counts.x)
          {
            ni -= self.counts.x;
            crossing.x = 1;
          }
          if (nj >= self.counts.y)
          {
            nj -= self.counts.y;
            crossing.y = 1;
          }
          if (nk >= self.counts.z)
          {
            nk -= self.counts.z;
            crossing.z = 1;
          }

          std::size_t there = self.indexOf(ni, nj, nk);
          if (self.clearance[there] < static_cast<float>(threshold)) continue;

          pieces.join(here, there, crossing);
        }
      }
    }
  }

  // The necks the grid stepped over.
  //
  // Two void voxels near each other but in different pieces are either side of something the flood could
  // not cross. Sometimes that something is a wall, and sometimes it is a passage the probe fits through
  // whose centre-line no voxel centre happened to land on. Which of the two it is can be settled without
  // refining anything: the room a sphere has along a straight segment is the distance from the nearest atom
  // to that segment less its radius, so one call says whether the whole line is passable.
  //
  // Only voxels against a wall start a search, since a voxel with void on all six sides is not next to
  // anything the flood failed to cross. Candidates are tried nearest first, because the nearest is the most
  // likely to be a passage and one success spares the rest.
  {
    double coarsest = std::max({self.spacing.x, self.spacing.y, self.spacing.z});
    double reach = std::max(1.0, 6.0 * coarsest);

    auto ringsAlong = [&](double step, std::int32_t count)
    {
      std::int32_t rings = static_cast<std::int32_t>(std::ceil(reach / std::max(step, 1.0e-9)));
      return std::max(1, std::min(rings, count / 2));
    };
    int3 rings(ringsAlong(self.spacing.x, self.counts.x), ringsAlong(self.spacing.y, self.counts.y),
               ringsAlong(self.spacing.z, self.counts.z));

    // How many lines one pair of pieces is worth trying. Two pieces separated by a wall thinner than the
    // reach offer a great many pairs and none of them will do, so the search is not allowed to spend the
    // run on a wall it has already failed to get through.
    const std::size_t mostTriesPerPair = 512;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> tried;

    struct Candidate
    {
      double distanceSquared;
      std::size_t voxel;
      int3 crossing;
    };
    std::vector<Candidate> candidates;

    const double3x3 &cell = structure.unitCell.cell;

    for (std::size_t here : voidVoxels)
    {
      std::int32_t i0 = static_cast<std::int32_t>(here % static_cast<std::size_t>(self.counts.x));
      std::int32_t j0 = static_cast<std::int32_t>((here / static_cast<std::size_t>(self.counts.x)) %
                                                  static_cast<std::size_t>(self.counts.y));
      std::int32_t k0 = static_cast<std::int32_t>(here / (static_cast<std::size_t>(self.counts.x) *
                                                          static_cast<std::size_t>(self.counts.y)));

      auto wrapped = [&](std::int32_t i, std::int32_t j, std::int32_t k)
      {
        return self.indexOf((i % self.counts.x + self.counts.x) % self.counts.x,
                            (j % self.counts.y + self.counts.y) % self.counts.y,
                            (k % self.counts.z + self.counts.z) % self.counts.z);
      };

      bool againstAWall = false;
      for (const int3 &step : steps)
      {
        if (self.clearance[wrapped(i0 + step.x, j0 + step.y, k0 + step.z)] < static_cast<float>(threshold) ||
            self.clearance[wrapped(i0 - step.x, j0 - step.y, k0 - step.z)] < static_cast<float>(threshold))
        {
          againstAWall = true;
          break;
        }
      }
      if (!againstAWall) continue;

      double3 from = cell * self.fractionalCentre(here);
      std::size_t rootHere = pieces.find(here).first;

      candidates.clear();

      for (std::int32_t dk = -rings.z; dk <= rings.z; ++dk)
      {
        for (std::int32_t dj = -rings.y; dj <= rings.y; ++dj)
        {
          for (std::int32_t di = -rings.x; di <= rings.x; ++di)
          {
            if (di == 0 && dj == 0 && dk == 0) continue;

            std::size_t there = wrapped(i0 + di, j0 + dj, k0 + dk);
            if (self.clearance[there] < static_cast<float>(threshold)) continue;
            if (pieces.find(there).first == rootHere) continue;

            // Which way round the cell the walk went, which is the translation the join carries.
            auto crossed = [](std::int32_t start, std::int32_t step, std::int32_t count)
            {
              return static_cast<std::int32_t>(std::floor(static_cast<double>(start + step) /
                                                          static_cast<double>(count)));
            };
            int3 crossing(crossed(i0, di, self.counts.x), crossed(j0, dj, self.counts.y),
                          crossed(k0, dk, self.counts.z));

            double3 to = cell * (self.fractionalCentre(there) +
                                 double3(static_cast<double>(crossing.x), static_cast<double>(crossing.y),
                                         static_cast<double>(crossing.z)));

            double3 gap = to - from;
            candidates.push_back(
                Candidate{.distanceSquared = double3::dot(gap, gap), .voxel = there, .crossing = crossing});
          }
        }
      }

      std::ranges::sort(candidates, {}, &Candidate::distanceSquared);

      for (const Candidate &candidate : candidates)
      {
        std::size_t rootThere = pieces.find(candidate.voxel).first;
        rootHere = pieces.find(here).first;
        if (rootThere == rootHere) continue;

        std::pair<std::size_t, std::size_t> pair{std::min(rootHere, rootThere), std::max(rootHere, rootThere)};
        std::size_t &spent = tried[pair];
        if (spent >= mostTriesPerPair) continue;
        ++spent;
        ++self.necksTried;

        double3 to = cell * (self.fractionalCentre(candidate.voxel) +
                             double3(static_cast<double>(candidate.crossing.x),
                                     static_cast<double>(candidate.crossing.y),
                                     static_cast<double>(candidate.crossing.z)));

        if (structure.segmentClearance(from, to - from) < threshold) continue;

        pieces.join(here, candidate.voxel, candidate.crossing);
        ++self.necksProved;
      }
    }
  }

  // Number the pieces in the order their roots are met, so that the numbering does not depend on the flood.
  std::unordered_map<std::size_t, std::int32_t> numbering;
  for (std::size_t voxel : voidVoxels)
  {
    std::size_t root = pieces.find(voxel).first;

    auto [entry, inserted] = numbering.try_emplace(root, static_cast<std::int32_t>(self.regions.size()));
    if (inserted)
    {
      BruteForceRegion region;
      region.dimensionality = latticeRank(pieces.cycles[root]);
      region.percolates = region.dimensionality > 0;
      region.deepestVoxel = voxel;
      region.deepestClearance = self.clearance[voxel];
      self.regions.push_back(region);
    }

    std::int32_t id = entry->second;
    self.regionOf[voxel] = id;

    BruteForceRegion &region = self.regions[static_cast<std::size_t>(id)];
    ++region.numberOfVoxels;
    if (self.clearance[voxel] > region.deepestClearance)
    {
      region.deepestClearance = self.clearance[voxel];
      region.deepestVoxel = voxel;
    }
  }

  for (BruteForceRegion &region : self.regions)
  {
    region.volume = static_cast<double>(region.numberOfVoxels) * self.voxelVolume;

    if (region.percolates)
      ++self.numberOfChannels;
    else
      ++self.numberOfPockets;
  }

  double perVoxel = 1.0 / static_cast<double>(numberOfVoxels);
  self.voidFraction = static_cast<double>(self.numberOfVoidVoxels) * perVoxel;

  std::size_t accessibleVoxels = 0;
  for (const BruteForceRegion &region : self.regions)
  {
    if (region.percolates) accessibleVoxels += region.numberOfVoxels;
  }
  self.accessibleFraction = static_cast<double>(accessibleVoxels) * perVoxel;
  self.blockedFraction = self.voidFraction - self.accessibleFraction;

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - begin;
  self.seconds = timing.count();

  return self;
}
