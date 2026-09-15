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
constexpr double walkStepFloor = 1.0e-10;
constexpr double diameterMergeTol = 1.0e-4;  // Å; walks converge far below this
constexpr double centreMergeTol = 1.0e-2;    // Å; same basin after wrapping
constexpr double weightFloor = 1.0e-5;       // same order as the exact spike floor
constexpr std::size_t maxFamiliesBeforeFilter = 12;

std::pair<double3, double> walkUphill(const BruteForceStructure &structure, double3 from, double startingStep)
{
  double3 best = from;
  double bestClearance = structure.clearance(best);
  double step = startingStep;

  while (step > walkStepFloor)
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

double micDistance(const BruteForceStructure &structure, const double3 &a, const double3 &b)
{
  return structure.nearestImage(a, b).length();
}

double3 wrapCentre(const BruteForceStructure &structure, const double3 &position)
{
  return structure.unitCell.cell * structure.wrappedFractional(position);
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
                                                   std::size_t maxFamilies,
                                                   const BruteForceStructure *reachabilityStructure,
                                                   const BruteForceVoxels *reachabilityVoxels)
{
  const bool blockPockets =
      reachabilityStructure != nullptr && reachabilityVoxels != nullptr;
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

  forEachIndex(seeds.size(), workers,
               [&](std::size_t worker, std::size_t index)
               {
                 auto [position, clearance] =
                     walkUphill(structure, voxels.centre(structure, seeds[index]), step);
                 if (clearance <= 0.0) return;
                 const double3 centre = wrapCentre(structure, position);
                 if (blockPockets && !reachabilityVoxels->isAccessible(*reachabilityStructure, centre)) return;
                 maximaOfWorker[worker].emplace_back(centre, clearance);
               });

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

  const std::size_t consider =
      std::min(families.size(), std::max(maxFamiliesBeforeFilter, maxFamilies));
  families.resize(consider);

  // One Monte Carlo draw: void membership and union membership per family.
  const std::size_t lanes = 64;
  const std::size_t perLane = std::max<std::size_t>(1, volumePoints / lanes);
  std::vector<std::size_t> inVoid(lanes, 0);
  std::vector<std::vector<std::size_t>> inUnion(lanes, std::vector<std::size_t>(families.size(), 0));

  forEachIndex(lanes, workersAvailable(),
               [&](std::size_t, std::size_t lane)
               {
                 RandomNumber random{lane + 17};  // offset from the pore-volume stream
                 for (std::size_t point = 0; point < perLane; ++point)
                 {
                   double3 fractional(random.uniform(), random.uniform(), random.uniform());
                   double3 position = structure.unitCell.cell * fractional;
                   if (structure.clearance(position) < 0.0) continue;
                   if (blockPockets && !reachabilityVoxels->isAccessible(*reachabilityStructure, position))
                     continue;
                   ++inVoid[lane];

                   // A point's pore size is the largest covering maximal sphere, so it belongs to at most
                   // one spike family — the largest-diameter family whose ball contains it. Counting it in
                   // every covering family would charge the volume of a deep cavity to every shallower
                   // spike whose ball happens to overlap that cavity.
                   for (std::size_t f = 0; f < families.size(); ++f)
                   {
                     const double radius = 0.5 * families[f].diameter;
                     bool covered = false;
                     for (const double3 &centre : families[f].centres)
                     {
                       if (micDistance(structure, position, centre) <= radius)
                       {
                         covered = true;
                         break;
                       }
                     }
                     if (covered)
                     {
                       ++inUnion[lane][f];
                       break;
                     }
                   }
                 }
               });

  const std::size_t voidPoints = std::reduce(inVoid.begin(), inVoid.end());
  const double cellVolume = structure.unitCell.volume;
  self.voidFraction = static_cast<double>(voidPoints) / static_cast<double>(perLane * lanes);
  self.voidVolume = self.voidFraction * cellVolume;

  if (voidPoints > 0)
  {
    for (std::size_t f = 0; f < families.size(); ++f)
    {
      std::size_t hits = 0;
      for (std::size_t lane = 0; lane < lanes; ++lane) hits += inUnion[lane][f];
      const double w = static_cast<double>(hits) / static_cast<double>(voidPoints);
      families[f].weight = w;
      families[f].weightError = std::sqrt(w * (1.0 - w) / static_cast<double>(voidPoints));
    }
  }

  // Drop families whose union is below the spike floor (wall corrugation / noise), keep largest three.
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
