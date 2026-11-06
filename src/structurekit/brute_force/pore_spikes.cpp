module;

module brute_force_pore_spikes;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import brute_force_structure;
import brute_force_voxels;
import structure_parallel;

namespace
{
constexpr double diameterMergeTol = 1.0e-4;  // Å; walks converge far below this
constexpr double centreMergeTol = 1.0e-2;    // Å; same basin after wrapping
constexpr double weightFloor = 1.0e-5;       // same order as the exact spike floor
constexpr std::size_t maxFamiliesBeforeFilter = 12;
// One walk is started from a point thrown at the cell for every so many points of the draw: eight thousand of
// them at the four million points the check throws by default, at which a family whose basins take up a
// thousandth of the cell gets eight and is missed about once in three thousand structures.
constexpr std::size_t pointsPerRandomWalk = 500;

// The walk uphill is the structure's own, shared with the diameter check, so that the two agree on where a
// maximal sphere sits and how large it is.
std::pair<double3, double> walkUphill(const BruteForceStructure &structure, double3 from, double startingStep)
{
  return structure.walkUphill(from, startingStep);
}

double micDistance(const BruteForceStructure &structure, const double3 &a, const double3 &b)
{
  return structure.nearestImage(a, b).length();
}

double3 wrapCentre(const BruteForceStructure &structure, const double3 &position)
{
  return structure.unitCell.cell * structure.wrappedFractional(position);
}

// The grid in blocks of this many voxels along each axis, and the most room at the middle of any voxel of each
// block: what lets a search for room pass over the greater part of the grid without reading it.
constexpr std::int32_t blockEdge = 4;
struct RoomiestOfBlocks
{
  int3 blocks;
  std::vector<float> most;

  explicit RoomiestOfBlocks(const BruteForceVoxels &voxels)
      : blocks((voxels.counts.x + blockEdge - 1) / blockEdge, (voxels.counts.y + blockEdge - 1) / blockEdge,
               (voxels.counts.z + blockEdge - 1) / blockEdge),
        most(static_cast<std::size_t>(blocks.x) * static_cast<std::size_t>(blocks.y) *
                 static_cast<std::size_t>(blocks.z),
             std::numeric_limits<float>::lowest())
  {
    for (std::int32_t k = 0; k < voxels.counts.z; ++k)
    {
      for (std::int32_t j = 0; j < voxels.counts.y; ++j)
      {
        for (std::int32_t i = 0; i < voxels.counts.x; ++i)
        {
          float &kept = this->most[this->indexOf(i / blockEdge, j / blockEdge, k / blockEdge)];
          kept = std::max(kept, voxels.clearance[voxels.indexOf(i, j, k)]);
        }
      }
    }
  }

  std::size_t indexOf(std::int32_t i, std::int32_t j, std::int32_t k) const
  {
    return static_cast<std::size_t>(i) +
           static_cast<std::size_t>(this->blocks.x) *
               (static_cast<std::size_t>(j) + static_cast<std::size_t>(this->blocks.y) * static_cast<std::size_t>(k));
  }
};

// What one lane of the draw keeps from point to point, so that nothing is allocated for each.
struct Workspace
{
  // A stretch of voxels along one axis that lies in one block of one copy of the cell: the first of them as
  // counted from the home cell, which may be before its start or past its end, and the same one wrapped.
  struct Run
  {
    std::int32_t first;
    std::int32_t wrapped;
    std::int32_t length;
    std::int32_t block;
  };

  // A voxel the grid leaves open: the room at its middle, and the atoms that have a say in the room anywhere
  // in it.
  struct Voxel
  {
    double room;
    std::size_t firstAtom;
    std::size_t numberOfAtoms;
    std::size_t lastToRule;
  };

  // A voxel or a part of one: where its middle is from the point, in fractions of the cell, how many times the
  // voxel has been halved to make it, and which voxel that was.
  struct Box
  {
    double3 middle;
    std::int32_t level;
    std::uint32_t voxel;
  };

  std::array<std::vector<Run>, 3> runs;
  std::vector<BruteForceStructure::Neighbour> atoms;
  std::vector<BruteForceStructure::Neighbour> atomsOfVoxels;
  std::vector<Voxel> voxels;
  std::vector<Box> boxes;
};

void runsAlong(std::int32_t first, std::int32_t last, std::int32_t count, std::vector<Workspace::Run> &runs)
{
  runs.clear();
  for (std::int32_t i = first; i <= last;)
  {
    const std::int32_t wrapped = ((i % count) + count) % count;
    const std::int32_t block = wrapped / blockEdge;
    const std::int32_t endOfBlock = std::min((block + 1) * blockEdge, count);
    const std::int32_t length = std::min(last - i + 1, endOfBlock - wrapped);
    runs.push_back(Workspace::Run{i, wrapped, length, block});
    i += length;
  }
}

// No more boxes than this are made by halving for one point. It is far more than a point needs unless the room
// is level over a stretch, as along the middle of a straight tube, where the question has no margin to find.
// A point that uses them up is settled by climbing instead; see below.
constexpr std::size_t mostBoxesHalved = 4096;

// Whether a sphere wider than `radius` that fits in the void holds `position`.
//
// A sphere of any size that fits is a union of spheres of a smaller size that fit, so one wider than the
// radius holds the point exactly when there is a centre within the radius of the point with more room than
// the radius round it. The room changes by no more than the distance moved, so a box with its middle too
// far from the point, or with too little room at its middle, for either to be made up within the box has no
// such centre in it, and a box with such a centre at its middle settles the matter. The grid gives the room at
// the middle of each voxel, which settles most points outright. The voxels it leaves open are halved along
// each axis, the room at the middle of each part taken from the atoms, and so on down, a level at a time, until
// a centre is found, or no box is left, or the boxes are smaller than `resolution` and whatever sphere there
// may still be is wider than the radius by less than that. `slack` is half the diagonal of a voxel, and
// `roomAtPosition` the room at the point, which the caller has already.
//
// Nothing is left to a search that could stop short: every box is either shown to hold no centre or opened.
//
// Except where the boxes run out. That happens where the room stays within a hair of the radius over a stretch,
// which is one point in ten thousand on a zeolite, and the halving has by then shown only that no sphere much
// wider than the radius holds the point. Answered no, as it used to be, the question went one way for both of
// the things asked of it: a family lost the points it held by so little, and a point held by a slightly wider
// sphere was taken not to be and counted for the family below. So the roomiest centre the halving came across
// is walked uphill inside the sphere of the radius about the point, and the answer is whether that gets above
// the radius, looked at once more by every atom of the structure before it is believed. Where the room is level
// the climb is at the level at once, and which side of the radius that is is the answer.
bool aWiderSphereHolds(const BruteForceStructure &structure, const BruteForceVoxels &voxels,
                       const RoomiestOfBlocks &roomiest, const double3 &position, double roomAtPosition,
                       double radius, double resolution, double slack, const double3 &widths, Workspace &workspace)
{
  // the sphere about the point itself
  if (roomAtPosition > radius) return true;

  const double3 s = structure.unitCell.inverseCell * position;
  const double3 counts(static_cast<double>(voxels.counts.x), static_cast<double>(voxels.counts.y),
                       static_cast<double>(voxels.counts.z));
  const double reach = radius + slack;

  // A step of one voxel along an axis is at least the width across that axis over the count in Ångström.
  const std::int32_t spanX = static_cast<std::int32_t>(std::ceil(reach * counts.x / widths.x));
  const std::int32_t spanY = static_cast<std::int32_t>(std::ceil(reach * counts.y / widths.y));
  const std::int32_t spanZ = static_cast<std::int32_t>(std::ceil(reach * counts.z / widths.z));
  const std::int32_t homeX = static_cast<std::int32_t>(std::floor(s.x * counts.x));
  const std::int32_t homeY = static_cast<std::int32_t>(std::floor(s.y * counts.y));
  const std::int32_t homeZ = static_cast<std::int32_t>(std::floor(s.z * counts.z));
  runsAlong(homeX - spanX, homeX + spanX, voxels.counts.x, workspace.runs[0]);
  runsAlong(homeY - spanY, homeY + spanY, voxels.counts.y, workspace.runs[1]);
  runsAlong(homeZ - spanZ, homeZ + spanZ, voxels.counts.z, workspace.runs[2]);

  std::vector<Workspace::Box> &boxes = workspace.boxes;
  boxes.clear();
  workspace.voxels.clear();
  double furthest = 0.0;  // of the atoms that have a say in a voxel left open, from the point
  for (const Workspace::Run &runZ : workspace.runs[2])
  {
    for (const Workspace::Run &runY : workspace.runs[1])
    {
      for (const Workspace::Run &runX : workspace.runs[0])
      {
        const float most = roomiest.most[roomiest.indexOf(runX.block, runY.block, runZ.block)];
        if (static_cast<double>(most) + slack <= radius) continue;

        // nor is a block read that lies beyond reach altogether: its middle, and no voxel of it further from
        // that than a whole block's half diagonal
        const double3 middleOfRuns(
            (static_cast<double>(runX.first) + 0.5 * static_cast<double>(runX.length)) / counts.x - s.x,
            (static_cast<double>(runY.first) + 0.5 * static_cast<double>(runY.length)) / counts.y - s.y,
            (static_cast<double>(runZ.first) + 0.5 * static_cast<double>(runZ.length)) / counts.z - s.z);
        const double tooFar = reach + static_cast<double>(blockEdge) * slack;
        const double3 toRuns = structure.unitCell.cell * middleOfRuns;
        if (double3::dot(toRuns, toRuns) > tooFar * tooFar) continue;

        for (std::int32_t dk = 0; dk < runZ.length; ++dk)
        {
          for (std::int32_t dj = 0; dj < runY.length; ++dj)
          {
            for (std::int32_t di = 0; di < runX.length; ++di)
            {
              const double room = static_cast<double>(
                  voxels.clearance[voxels.indexOf(runX.wrapped + di, runY.wrapped + dj, runZ.wrapped + dk)]);
              if (room + slack <= radius) continue;

              const double3 middle((static_cast<double>(runX.first + di) + 0.5) / counts.x - s.x,
                                   (static_cast<double>(runY.first + dj) + 0.5) / counts.y - s.y,
                                   (static_cast<double>(runZ.first + dk) + 0.5) / counts.z - s.z);
              const double3 offset = structure.unitCell.cell * middle;
              const double distanceSquared = double3::dot(offset, offset);
              if (distanceSquared > reach * reach) continue;
              if (room > radius && distanceSquared <= radius * radius) return true;
              boxes.push_back(Workspace::Box{middle, 0, static_cast<std::uint32_t>(boxes.size())});
              workspace.voxels.push_back(Workspace::Voxel{room, 0, 0, 0});
              furthest = std::max(furthest, std::sqrt(distanceSquared) + room);
            }
          }
        }
      }
    }
  }
  if (boxes.empty()) return false;

  // Anywhere in a voxel the nearest atom is no further than the room at the voxel's middle and half its
  // diagonal, and an atom further from the middle than that by another half is nearer nowhere in the voxel
  // than the one that is nearest at the middle. So each voxel left open gets the few atoms within its own
  // diagonal of being nearest at its middle, the room among those is the room everywhere in it, and the parts
  // it is halved into ask only those. The grid keeps the room in single precision, which the margin allows for.
  constexpr double roundingOfTheGrid = 1.0e-5;
  std::vector<BruteForceStructure::Neighbour> &atoms = workspace.atoms;
  const bool nearbyOnly = structure.atomsWithin(position, furthest + 2.0 * slack + roundingOfTheGrid, atoms);
  if (nearbyOnly && atoms.empty()) return true;

  workspace.atomsOfVoxels.clear();
  if (nearbyOnly)
  {
    for (std::size_t v = 0; v < workspace.voxels.size(); ++v)
    {
      Workspace::Voxel &voxel = workspace.voxels[v];
      const double3 offset = structure.unitCell.cell * boxes[v].middle;
      const double within = voxel.room + 2.0 * slack + roundingOfTheGrid;
      voxel.firstAtom = workspace.atomsOfVoxels.size();
      for (const BruteForceStructure::Neighbour &atom : atoms)
      {
        const double3 fromCentre = atom.fromCentre + offset;
        const double limit = within + atom.radius;
        if (limit >= 0.0 && double3::dot(fromCentre, fromCentre) <= limit * limit)
        {
          workspace.atomsOfVoxels.push_back(atom);
        }
      }
      voxel.numberOfAtoms = workspace.atomsOfVoxels.size() - voxel.firstAtom;
    }
  }

  // The room at an offset from the point, in a part of the given voxel, or, once it is known to be no more than
  // `enough`, whatever it has come down to by then. The atom that ended the last look in the voxel is the first
  // asked in the next: boxes side by side are mostly ruled out by the same one.
  auto roomAt = [&](const double3 &offset, double enough, std::uint32_t voxel)
  {
    if (!nearbyOnly) return structure.clearance(position + offset);

    Workspace::Voxel &home = workspace.voxels[voxel];
    if (home.numberOfAtoms == 0) return std::numeric_limits<double>::max();
    const BruteForceStructure::Neighbour *list = workspace.atomsOfVoxels.data() + home.firstAtom;

    double least = (list[home.lastToRule].fromCentre + offset).length() - list[home.lastToRule].radius;
    if (least <= enough) return least;
    for (std::size_t n = 0; n < home.numberOfAtoms; ++n)
    {
      const double here = (list[n].fromCentre + offset).length() - list[n].radius;
      if (here < least)
      {
        least = here;
        if (least <= enough)
        {
          home.lastToRule = n;
          break;
        }
      }
    }
    return least;
  };

  // The roomiest centre inside the sphere that the halving comes across, for the climb if the boxes run out.
  double roomiestSeen = roomAtPosition;
  double3 roomiestOffset(0.0, 0.0, 0.0);

  const double3 quarter(0.25 / counts.x, 0.25 / counts.y, 0.25 / counts.z);
  const std::size_t voxelsLeftOpen = boxes.size();
  bool boxesRanOut = false;
  for (std::size_t next = 0; next < boxes.size(); ++next)
  {
    const Workspace::Box box = boxes[next];
    const double scale = std::ldexp(1.0, -box.level);
    const double halfDiagonal = 0.5 * slack * scale;  // of each of the eight parts
    for (std::int32_t part = 0; part < 8; ++part)
    {
      const double3 middle(box.middle.x + ((part & 1) ? quarter.x : -quarter.x) * scale,
                           box.middle.y + ((part & 2) ? quarter.y : -quarter.y) * scale,
                           box.middle.z + ((part & 4) ? quarter.z : -quarter.z) * scale);
      const double3 offset = structure.unitCell.cell * middle;
      const double distance = offset.length();
      if (distance - halfDiagonal > radius) continue;
      const double room = roomAt(offset, radius - halfDiagonal, box.voxel);
      if (room > radius && distance <= radius) return true;
      if (room + halfDiagonal <= radius) continue;
      if (distance <= radius && room > roomiestSeen)
      {
        roomiestSeen = room;
        roomiestOffset = offset;
      }
      if (halfDiagonal > resolution) boxes.push_back(Workspace::Box{middle, box.level + 1, box.voxel});
    }
    if (boxes.size() - voxelsLeftOpen > mostBoxesHalved)
    {
      boxesRanOut = true;
      break;
    }
  }
  if (!boxesRanOut) return false;

  // The room at any offset inside the sphere, by every atom gathered for it.
  auto roomAnywhere = [&](const double3 &offset)
  {
    if (!nearbyOnly) return structure.clearance(position + offset);
    double least = std::numeric_limits<double>::max();
    for (const BruteForceStructure::Neighbour &atom : atoms)
    {
      least = std::min(least, (atom.fromCentre + offset).length() - atom.radius);
    }
    return least;
  };

  double3 best = roomiestOffset;
  double bestRoom = roomAnywhere(best);
  for (double step = slack; step > 0.25 * resolution && bestRoom <= radius;)
  {
    bool moved = false;
    for (std::int32_t k = -1; k <= 1 && bestRoom <= radius; ++k)
    {
      for (std::int32_t j = -1; j <= 1 && bestRoom <= radius; ++j)
      {
        for (std::int32_t i = -1; i <= 1 && bestRoom <= radius; ++i)
        {
          if (i == 0 && j == 0 && k == 0) continue;
          const double3 direction(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
          const double3 trial = best + (step / direction.length()) * direction;
          if (double3::dot(trial, trial) > radius * radius) continue;
          const double room = roomAnywhere(trial);
          if (room > bestRoom)
          {
            bestRoom = room;
            best = trial;
            moved = true;
          }
        }
      }
    }
    if (!moved) step *= 0.5;
  }

  // The atoms gathered are those that have a say where the halving looked, and the climb may have left there:
  // so what it found is asked of the whole structure before it is taken for a sphere that fits.
  return bestRoom > radius && structure.clearance(position + best) > radius;
}

bool isGridLocalMaximum(const BruteForceVoxels &voxels, std::size_t voxel)
{
  if (voxels.regionOf[voxel] < 0) return false;

  const std::int32_t nx = voxels.counts.x;
  const std::int32_t ny = voxels.counts.y;
  const std::int32_t nz = voxels.counts.z;
  const std::int32_t i0 = static_cast<std::int32_t>(voxel % static_cast<std::size_t>(nx));
  const std::int32_t j0 =
      static_cast<std::int32_t>((voxel / static_cast<std::size_t>(nx)) % static_cast<std::size_t>(ny));
  const std::int32_t k0 = static_cast<std::int32_t>(voxel / (static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny)));
  const float self = voxels.clearance[voxel];

  for (std::int32_t dk = -1; dk <= 1; ++dk)
  {
    for (std::int32_t dj = -1; dj <= 1; ++dj)
    {
      for (std::int32_t di = -1; di <= 1; ++di)
      {
        if (di == 0 && dj == 0 && dk == 0) continue;
        std::int32_t i = i0 + di;
        std::int32_t j = j0 + dj;
        std::int32_t k = k0 + dk;
        // Periodic wrap.
        if (i < 0) i += nx;
        if (j < 0) j += ny;
        if (k < 0) k += nz;
        if (i >= nx) i -= nx;
        if (j >= ny) j -= ny;
        if (k >= nz) k -= nz;
        std::size_t neighbour = voxels.indexOf(i, j, k);
        if (voxels.regionOf[neighbour] >= 0 && voxels.clearance[neighbour] > self) return false;
      }
    }
  }
  return true;
}
}  // namespace

BruteForcePoreSpikes BruteForcePoreSpikes::compute(const BruteForceStructure &structure,
                                                   const BruteForceVoxels &voxels, std::size_t volumePoints,
                                                   std::size_t maxFamilies)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  BruteForcePoreSpikes self;

  const std::size_t numberOfVoxels = voxels.clearance.size();
  if (numberOfVoxels == 0 || voxels.numberOfVoidVoxels == 0) return self;

  // Seeds: every grid local maximum, plus the roomiest voxels (as for Di).
  std::vector<std::size_t> seeds;
  seeds.reserve(voxels.numberOfVoidVoxels);
  for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
  {
    if (isGridLocalMaximum(voxels, voxel)) seeds.push_back(voxel);
  }

  std::vector<std::size_t> ranked;
  ranked.reserve(voxels.numberOfVoidVoxels);
  for (std::size_t voxel = 0; voxel < numberOfVoxels; ++voxel)
  {
    if (voxels.regionOf[voxel] >= 0) ranked.push_back(voxel);
  }
  const std::size_t topStarts = std::min<std::size_t>(64, ranked.size());
  if (topStarts > 0)
  {
    std::partial_sort(ranked.begin(), ranked.begin() + static_cast<std::ptrdiff_t>(topStarts), ranked.end(),
                      [&](std::size_t a, std::size_t b) { return voxels.clearance[a] > voxels.clearance[b]; });
    for (std::size_t i = 0; i < topStarts; ++i) seeds.push_back(ranked[i]);
  }

  std::sort(seeds.begin(), seeds.end());
  seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());

  const double step = std::max({voxels.spacing.x, voxels.spacing.y, voxels.spacing.z});
  const std::size_t workers = workersAvailable();
  std::vector<std::vector<std::pair<double3, double>>> maximaOfWorker(workers);

  // Each seed is walked from twice. A walk that starts with a step of a whole voxel strides over whatever is
  // narrower than that, which is what takes it past the ripples of a wall to the sphere worth finding, and
  // also what takes it out of a small pore into the larger one beside it: a crevice half an Ångström across
  // is a family of its own and no stride of that length ends in it. So a second walk starts with a step too
  // short to leave the basin the seed is in.
  const std::array<double, 2> startingSteps{step, 0.125 * step};
  forEachIndex(seeds.size(), workers,
               [&](std::size_t worker, std::size_t index)
               {
                 for (double startingStep : startingSteps)
                 {
                   auto [position, clearance] =
                       walkUphill(structure, voxels.centre(structure, seeds[index]), startingStep);
                   if (clearance <= 0.0) continue;
                   const double3 centre = wrapCentre(structure, position);
                   maximaOfWorker[worker].emplace_back(centre, clearance);
                 }
               });

  // And walks from points thrown at the void, for the maxima the grid has no seed in. A voxel is a seed when
  // none of its neighbours has more room, and a rise of a few hundredths beside the slope up to something
  // larger never has such a voxel: the neighbour on the slope always has more. A point thrown into its basin
  // climbs to it all the same, so long as its first step is too short to leave.
  {
    const std::size_t lanes = 64;
    const std::size_t perLane = (volumePoints / pointsPerRandomWalk + lanes - 1) / lanes;
    std::vector<std::vector<std::pair<double3, double>>> maximaOfLane(lanes);
    forEachIndex(lanes, workers,
                 [&](std::size_t, std::size_t lane)
                 {
                   RandomNumber random{lane + 17 + 2 * 64};  // apart from the streams of the draws
                   for (std::size_t walk = 0; walk < perLane; ++walk)
                   {
                     const double3 from =
                         structure.unitCell.cell * double3(random.uniform(), random.uniform(), random.uniform());
                     if (structure.clearance(from) <= 0.0) continue;
                     auto [position, clearance] = walkUphill(structure, from, 0.25 * step);
                     if (clearance <= 0.0) continue;
                     maximaOfLane[lane].emplace_back(wrapCentre(structure, position), clearance);
                   }
                 });
    for (auto &bucket : maximaOfLane)
    {
      maximaOfWorker.front().insert(maximaOfWorker.front().end(), bucket.begin(), bucket.end());
    }
  }

  std::vector<std::pair<double3, double>> maxima;
  for (auto &bucket : maximaOfWorker)
  {
    maxima.insert(maxima.end(), bucket.begin(), bucket.end());
    bucket.clear();
  }

  // Deduplicate basins: same centre and same clearance after the walk.
  std::vector<std::pair<double3, double>> unique;
  unique.reserve(maxima.size());
  for (const auto &[centre, clearance] : maxima)
  {
    bool seen = false;
    for (auto &[keptCentre, keptClearance] : unique)
    {
      if (std::abs(clearance - keptClearance) <= 0.5 * diameterMergeTol &&
          micDistance(structure, centre, keptCentre) <= centreMergeTol)
      {
        if (clearance > keptClearance)
        {
          keptCentre = centre;
          keptClearance = clearance;
        }
        seen = true;
        break;
      }
    }
    if (!seen) unique.emplace_back(centre, clearance);
  }
  self.numberOfMaxima = unique.size();

  std::sort(unique.begin(), unique.end(),
            [](const auto &a, const auto &b) { return a.second > b.second; });

  // Cluster into diameter families (largest first).
  std::vector<BruteForceSpikeFamily> families;
  for (const auto &[centre, clearance] : unique)
  {
    const double diameter = 2.0 * clearance;
    if (!families.empty() && std::abs(families.back().diameter - diameter) <= diameterMergeTol)
    {
      families.back().centres.push_back(centre);
      // Keep the diameter of the deepest member.
      if (diameter > families.back().diameter) families.back().diameter = diameter;
    }
    else
    {
      BruteForceSpikeFamily family;
      family.diameter = diameter;
      family.centres.push_back(centre);
      families.push_back(std::move(family));
    }
  }
  self.numberOfFamilies = families.size();

  // One Monte Carlo draw over a set of families, largest diameter first: how many points land in the void,
  // and how many of those have each family's diameter for their pore size. A point's pore size is that of the
  // widest sphere that fits in the void and holds it, so it has a family's diameter when a sphere as wide as
  // the family's holds it and no wider one does, and it can have one family's at most.
  //
  // The centres the walks found are only where the question is answered at once. A sphere of a family's size
  // may stand where no walk ended: two maxima a tenth of an Ångström apart, or one whose rise is a few
  // thousandths, are less than the grid the walks start from can tell apart, and the void they hold is still
  // void of that pore size. So a point no centre of the family reaches is asked about directly, and what the
  // walks are relied on for is the diameters and not every place they occur. For the widest family of all
  // there is no wider sphere to ask about.
  const double3 widths = structure.unitCell.perpendicularWidths();
  const double halfVoxelDiagonal =
      0.5 * structure.unitCell.longestDiagonal(double3(1.0 / static_cast<double>(voxels.counts.x),
                                                       1.0 / static_cast<double>(voxels.counts.y),
                                                       1.0 / static_cast<double>(voxels.counts.z)));
  const RoomiestOfBlocks roomiest(voxels);
  struct Draw
  {
    std::size_t thrown{0};
    std::size_t voidPoints{0};
    std::vector<std::size_t> hits;
  };
  auto draw = [&](std::size_t points, const std::vector<BruteForceSpikeFamily> &drawn, std::size_t seedOffset)
  {
    const std::size_t lanes = 64;
    const std::size_t perLane = std::max<std::size_t>(1, points / lanes);
    std::vector<std::size_t> inVoid(lanes, 0);
    std::vector<std::vector<std::size_t>> inUnion(lanes, std::vector<std::size_t>(drawn.size(), 0));

    forEachIndex(lanes, workersAvailable(),
                 [&](std::size_t, std::size_t lane)
                 {
                   RandomNumber random{lane + seedOffset};
                   Workspace workspace;

                   // How finely the search for a wider sphere looks, and by how much a radius is eased so that a
                   // family is found by it: the spread of a family's radii and the resolution together.
                   constexpr double searchResolution = 0.5 * diameterMergeTol;
                   constexpr double radiusAllowance = 0.5 * diameterMergeTol + searchResolution;
                   for (std::size_t point = 0; point < perLane; ++point)
                   {
                     double3 fractional(random.uniform(), random.uniform(), random.uniform());
                     double3 position = structure.unitCell.cell * fractional;
                     const double room = structure.clearance(position);
                     if (room < 0.0) continue;
                     ++inVoid[lane];

                     // Whether a sphere as wide as the family's holds the point: as wide to within what the
                     // walks agree to, so that a family does not miss itself.
                     //
                     // The allowance is a radius, and it is made of two halves of the tolerance on a diameter. The
                     // spheres of one family differ in diameter by that tolerance at most, so in radius by half
                     // of it; and the search answers for spheres wider than the radius asked by more than its
                     // own resolution, which is the other half. A sphere of the family is therefore always wider
                     // than what is asked here by enough to be found, and no allowance smaller than this does that.
                     auto heldBy = [&](std::size_t f)
                     {
                       const double radius = 0.5 * drawn[f].diameter;
                       for (const double3 &centre : drawn[f].centres)
                       {
                         if (micDistance(structure, position, centre) <= radius) return true;
                       }
                       return aWiderSphereHolds(structure, voxels, roomiest, position, room,
                                                radius - radiusAllowance, searchResolution, halfVoxelDiagonal,
                                                widths, workspace);
                     };

                     // A point held by a sphere of one size is held by one of every smaller size, so down the
                     // list the families say no and then yes, and the first to say yes is found by halving.
                     std::size_t first = 0;
                     std::size_t beyond = drawn.size();
                     while (first < beyond)
                     {
                       const std::size_t middle = first + (beyond - first) / 2;
                       if (heldBy(middle))
                         beyond = middle;
                       else
                         first = middle + 1;
                     }
                     if (first == drawn.size()) continue;

                     // and wider than the family's by more than the walks agree to, so that it is not held by
                     // itself
                     const bool heldWider =
                         first > 0 && aWiderSphereHolds(structure, voxels, roomiest, position, room,
                                                        0.5 * drawn[first].diameter + radiusAllowance,
                                                        searchResolution, halfVoxelDiagonal, widths, workspace);
                     if (!heldWider) ++inUnion[lane][first];
                   }
                 });

    Draw result;
    result.thrown = perLane * lanes;
    result.voidPoints = std::reduce(inVoid.begin(), inVoid.end());
    result.hits.assign(drawn.size(), 0);
    for (std::size_t lane = 0; lane < lanes; ++lane)
    {
      for (std::size_t f = 0; f < drawn.size(); ++f) result.hits[f] += inUnion[lane][f];
    }
    return result;
  };

  // The weight is what decides which families are reported: the ones below the spike floor are wall
  // corrugation and noise, and of the rest the widest few are wanted. The full draw over every family would
  // be paid for mostly by the points that land in no ball and are held against every centre there is, so
  // when there are more families than are worth weighing in full, a lighter draw goes first and says how far
  // down the list the reported ones can lie.
  //
  // What can be left out of the full draw is only what lies past that. A point goes to the widest family that
  // covers it, so a family has no say in the weight of one wider than itself, and everything narrower than the
  // last one reported can go. Nothing wider can, however little it holds: its void would fall through to the
  // narrower families reaching over it, and a family holding little is still reported if it is wide enough. Nor
  // is the first draw fine enough to place a family close to the floor, where it has only a few points to go by,
  // so it is only taken at its word well above the floor and the full draw settles the rest.
  const std::size_t consider = std::max(maxFamiliesBeforeFilter, maxFamilies);
  if (families.size() > consider)
  {
    const std::size_t rankingPoints = std::min(volumePoints, std::max<std::size_t>(200000, volumePoints / 10));
    Draw ranking = draw(rankingPoints, families, 17 + 64);  // its own stream, apart from the full draw's

    constexpr double clearlyAboveTheFloor = 3.0 * weightFloor;
    const double enough = clearlyAboveTheFloor * static_cast<double>(ranking.voidPoints);
    std::size_t clearlyAbove = 0;
    for (std::size_t f = 0; f < families.size(); ++f)
    {
      if (ranking.voidPoints > 0 && static_cast<double>(ranking.hits[f]) >= enough && ++clearlyAbove == maxFamilies)
      {
        families.resize(f + 1);
        break;
      }
    }
  }

  Draw full = draw(volumePoints, families, 17);  // offset from the pore-volume stream
  const double cellVolume = structure.unitCell.volume;
  self.voidFraction = static_cast<double>(full.voidPoints) / static_cast<double>(full.thrown);
  self.voidVolume = self.voidFraction * cellVolume;

  if (full.voidPoints > 0)
  {
    for (std::size_t f = 0; f < families.size(); ++f)
    {
      const double w = static_cast<double>(full.hits[f]) / static_cast<double>(full.voidPoints);
      families[f].weight = w;
      families[f].weightError = std::sqrt(w * (1.0 - w) / static_cast<double>(full.voidPoints));
    }
  }

  // Drop families whose union is below the spike floor (wall corrugation / noise), keep the widest few. They are in
  // order of diameter, widest first, whatever went before.
  std::vector<BruteForceSpikeFamily> kept;
  kept.reserve(maxFamilies);
  for (auto &family : families)
  {
    if (family.weight < weightFloor) continue;
    kept.push_back(std::move(family));
    if (kept.size() == maxFamilies) break;
  }
  self.families = std::move(kept);

  self.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
  return self;
}
