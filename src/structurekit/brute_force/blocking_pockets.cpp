module;

module brute_force_blocking_pockets;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import brute_force_structure;
import brute_force_voxels;
import structure_parallel;

BruteForceBlockingPockets BruteForceBlockingPockets::compute(const BruteForceStructure &structure,
                                                             const BruteForceVoxels &voxels)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  BruteForceBlockingPockets self;

  // Which voxels belong to which piece, so that each pocket can be walked without going over the grid once
  // per pocket.
  std::vector<std::vector<std::size_t>> voxelsOfRegion(voxels.regions.size());
  for (std::size_t voxel = 0; voxel < voxels.regionOf.size(); ++voxel)
  {
    std::int32_t region = voxels.regionOf[voxel];
    if (region >= 0) voxelsOfRegion[static_cast<std::size_t>(region)].push_back(voxel);
  }

  std::vector<std::size_t> channelVoxels;
  for (std::size_t region = 0; region < voxels.regions.size(); ++region)
  {
    if (voxels.regions[region].percolates)
    {
      channelVoxels.insert(channelVoxels.end(), voxelsOfRegion[region].begin(), voxelsOfRegion[region].end());
    }
  }

  for (std::size_t region = 0; region < voxels.regions.size(); ++region)
  {
    if (voxels.regions[region].percolates) continue;

    const std::vector<std::size_t> &members = voxelsOfRegion[region];
    if (members.empty()) continue;

    BruteForcePocket pocket;
    pocket.numberOfVoxels = members.size();
    pocket.volume = static_cast<double>(members.size()) * voxels.voxelVolume;

    // A pocket does not run away, but it can still straddle a boundary, so the centre is taken about one of
    // its own voxels with everything brought to the nearest image of that. Wrapping first and averaging
    // after would put the centre of a pocket sitting on a face somewhere in the middle of the cell.
    double3 anchor = voxels.centre(structure, members.front());

    double3 sum(0.0, 0.0, 0.0);
    for (std::size_t voxel : members)
    {
      sum += structure.nearestImage(anchor, voxels.centre(structure, voxel));
    }
    double3 centre = anchor + (1.0 / static_cast<double>(members.size())) * sum;

    pocket.centreFractional = structure.wrappedFractional(centre);
    pocket.centre = structure.unitCell.cell * pocket.centreFractional;
    pocket.freeRadius = std::max(structure.clearance(pocket.centre), 0.0);

    // How far the pocket reaches from its centre, and how near the nearest place a molecule belongs is. A
    // voxel stands for a little box, so the reach is taken to the far corner of it and the approach to the
    // near corner, which keeps the sphere covering what it should on a coarse grid.
    double halfDiagonal = 0.5 * std::sqrt(voxels.spacing.x * voxels.spacing.x + voxels.spacing.y * voxels.spacing.y +
                                          voxels.spacing.z * voxels.spacing.z);

    for (std::size_t voxel : members)
    {
      double distance = structure.nearestImage(pocket.centre, voxels.centre(structure, voxel)).length();
      pocket.coveringRadius = std::max(pocket.coveringRadius, distance + halfDiagonal);
    }

    for (std::size_t voxel : channelVoxels)
    {
      double distance = structure.nearestImage(pocket.centre, voxels.centre(structure, voxel)).length();
      pocket.channelRadius = std::min(pocket.channelRadius, std::max(distance - halfDiagonal, 0.0));
    }

    self.inaccessibleVolume += pocket.volume;
    self.pockets.push_back(pocket);
  }

  std::ranges::sort(self.pockets, [](const BruteForcePocket &a, const BruteForcePocket &b)
                    { return a.volume > b.volume; });

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - begin;
  self.seconds = timing.count();

  return self;
}

BlockingAudit auditBlockingSpheres(const BruteForceStructure &structure, const BruteForceVoxels &voxels,
                                   const std::vector<std::pair<double3, double>> &spheres)
{
  BlockingAudit audit;
  audit.numberOfSpheres = spheres.size();

  std::vector<double3> centres;
  centres.reserve(spheres.size());
  for (const auto &[fractional, radius] : spheres)
  {
    centres.push_back(structure.unitCell.cell * fractional);
  }

  std::size_t numberOfVoxels = voxels.regionOf.size();

  // Each worker audits a stretch of the voxels into counters of its own. What is collected is two counts and a
  // maximum, none of which depends on the order they were collected in, so this comes out the same to the bit
  // on any number of threads.
  const std::size_t workers = workersAvailable();
  std::vector<std::size_t> leftOpenOfWorker(workers, 0);
  std::vector<std::size_t> coveredUpOfWorker(workers, 0);
  std::vector<double> worstOfWorker(workers, 0.0);

  forEachBlock(numberOfVoxels, workers,
               [&](std::size_t worker, std::size_t begin, std::size_t end)
               {
                 for (std::size_t voxel = begin; voxel < end; ++voxel)
                 {
                   std::int32_t region = voxels.regionOf[voxel];
                   if (region < 0) continue;

                   bool percolates = voxels.regions[static_cast<std::size_t>(region)].percolates;
                   double3 position = voxels.centre(structure, voxel);

                   bool covered = false;
                   double deepest = 0.0;

                   for (std::size_t sphere = 0; sphere < centres.size(); ++sphere)
                   {
                     double distance = structure.nearestImage(position, centres[sphere]).length();
                     if (distance < spheres[sphere].second)
                     {
                       covered = true;
                       deepest = std::max(deepest, spheres[sphere].second - distance);
                     }
                   }

                   if (percolates)
                   {
                     if (covered)
                     {
                       ++coveredUpOfWorker[worker];
                       worstOfWorker[worker] = std::max(worstOfWorker[worker], deepest);
                     }
                   }
                   else if (!covered)
                   {
                     ++leftOpenOfWorker[worker];
                   }
                 }
               });

  std::size_t leftOpen = std::reduce(leftOpenOfWorker.begin(), leftOpenOfWorker.end());
  std::size_t coveredUp = std::reduce(coveredUpOfWorker.begin(), coveredUpOfWorker.end());
  double worst = std::ranges::max(worstOfWorker);

  audit.pocketVoxelsLeftOpen = leftOpen;
  audit.channelVoxelsCoveredUp = coveredUp;
  audit.volumeLeftOpen = static_cast<double>(leftOpen) * voxels.voxelVolume;
  audit.volumeCoveredUp = static_cast<double>(coveredUp) * voxels.voxelVolume;
  audit.worstOverreach = worst;

  return audit;
}
