module;

module brute_force_pore_volume;

import std;

import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import brute_force_structure;
import brute_force_voxels;

BruteForcePoreVolume BruteForcePoreVolume::compute(const BruteForceStructure &structure,
                                                   const BruteForceVoxels &voxels, std::size_t numberOfPoints)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  BruteForcePoreVolume self;
  self.numberOfPoints = numberOfPoints;

  const std::size_t lanes = 64;
  std::vector<std::size_t> inVoid(lanes, 0);
  std::vector<std::size_t> inChannel(lanes, 0);
  std::vector<std::size_t> inPocket(lanes, 0);
  std::vector<std::size_t> unresolved(lanes, 0);

  std::size_t perLane = numberOfPoints / lanes;

#pragma omp parallel for schedule(dynamic, 1)
  for (std::int64_t index = 0; index < static_cast<std::int64_t>(lanes); ++index)
  {
    std::size_t lane = static_cast<std::size_t>(index);
    RandomNumber random{lane};

    for (std::size_t point = 0; point < perLane; ++point)
    {
      double3 fractional(random.uniform(), random.uniform(), random.uniform());
      double3 position = structure.unitCell.cell * fractional;

      if (structure.clearance(position) < 0.0) continue;
      ++inVoid[lane];

      // Which piece of the void the point is in comes from the flood, which is the only thing here that can
      // tell a channel from a pocket. A point can have room around it and still fall in a voxel whose own
      // centre does not, so the nearest labelled voxel is what answers it; where even that fails the point
      // is in void too thin for the grid to have found at all, and it is counted as such rather than given
      // to either side.
      std::int32_t region = voxels.regionNear(structure, position);
      if (region < 0)
      {
        ++unresolved[lane];
        continue;
      }

      if (voxels.regions[static_cast<std::size_t>(region)].percolates)
        ++inChannel[lane];
      else
        ++inPocket[lane];
    }
  }

  std::size_t drawn = perLane * lanes;
  std::size_t voidPoints = std::reduce(inVoid.begin(), inVoid.end());
  std::size_t channelPoints = std::reduce(inChannel.begin(), inChannel.end());
  std::size_t pocketPoints = std::reduce(inPocket.begin(), inPocket.end());
  std::size_t unresolvedPoints = std::reduce(unresolved.begin(), unresolved.end());

  self.numberOfPoints = drawn;

  double perPoint = 1.0 / static_cast<double>(drawn);
  self.voidFraction = static_cast<double>(voidPoints) * perPoint;
  self.accessibleFraction = static_cast<double>(channelPoints) * perPoint;
  self.inaccessibleFraction = static_cast<double>(pocketPoints) * perPoint;
  self.unresolvedFraction = static_cast<double>(unresolvedPoints) * perPoint;

  self.voidFractionError = std::sqrt(self.voidFraction * (1.0 - self.voidFraction) * perPoint);

  double volume = structure.unitCell.volume;
  self.voidVolume = self.voidFraction * volume;
  self.accessibleVolume = self.accessibleFraction * volume;
  self.inaccessibleVolume = self.inaccessibleFraction * volume;
  self.unresolvedVolume = self.unresolvedFraction * volume;

  self.voidFractionFromVoxels = voxels.voidFraction;
  self.accessibleFractionFromVoxels = voxels.accessibleFraction;
  self.inaccessibleFractionFromVoxels = voxels.blockedFraction;

  self.numberOfChannels = voxels.numberOfChannels;
  self.numberOfPockets = voxels.numberOfPockets;
  self.necksProved = voxels.necksProved;
  self.necksTried = voxels.necksTried;

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - begin;
  self.seconds = timing.count();

  return self;
}
