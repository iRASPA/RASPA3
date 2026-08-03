module;

export module brute_force_blocking_pockets;

import std;

import double3;
import brute_force_structure;
import brute_force_voxels;

// One sealed pocket, as the flooded voxels have it.
export struct BruteForcePocket
{
  double3 centre{};             // Cartesian, in the home cell
  double3 centreFractional{};

  double volume{0.0};  // Å³
  std::size_t numberOfVoxels{0};

  double freeRadius{0.0};      // Å, the room at the centre
  double coveringRadius{0.0};  // Å, out to the farthest of the pocket's own voxels
  double channelRadius{std::numeric_limits<double>::max()};  // Å, in to the nearest voxel of a channel

  bool hasChannel() const { return channelRadius < 0.5 * std::numeric_limits<double>::max(); }
  double blockingRadius() const { return std::min(coveringRadius, channelRadius); }
  bool coversPocket() const { return coveringRadius <= channelRadius; }
};

// What a set of blocking spheres, from wherever, does to the void the flood found. This is the check that
// matters about a blocking sphere, and it is a check the exact route cannot make on itself: it says whether
// the spheres in fact cover what a molecule must be kept out of and in fact leave alone what it is allowed.
export struct BlockingAudit
{
  std::size_t numberOfSpheres{0};

  std::size_t pocketVoxelsLeftOpen{0};   // sealed void no sphere covers, so a molecule could be put there
  double volumeLeftOpen{0.0};            // Å³

  std::size_t channelVoxelsCoveredUp{0};  // void a molecule belongs in that the spheres refuse
  double volumeCoveredUp{0.0};            // Å³

  double worstOverreach{0.0};  // Å, how far the deepest covered channel voxel is inside its sphere

  bool clean() const { return pocketVoxelsLeftOpen == 0 && channelVoxelsCoveredUp == 0; }
};

export struct BruteForceBlockingPockets
{
  std::vector<BruteForcePocket> pockets;

  double inaccessibleVolume{0.0};  // Å³, over all of them
  double seconds{0.0};

  static BruteForceBlockingPockets compute(const BruteForceStructure &structure, const BruteForceVoxels &voxels);
};

// Holds a set of spheres, given as fractional centres and radii in Å, against the flood.
export BlockingAudit auditBlockingSpheres(const BruteForceStructure &structure, const BruteForceVoxels &voxels,
                                          const std::vector<std::pair<double3, double>> &spheres);
