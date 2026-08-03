module;

export module brute_force_validation;

import std;

import double3;
import pair_interactions;
import crystal;
import brute_force_structure;
import brute_force_voxels;
import brute_force_diameters;
import brute_force_surface_area;
import brute_force_solvent_excluded;
import brute_force_pore_volume;
import brute_force_blocking_pockets;

// How hard the checks are run. Everything here trades time for how tightly the brute force can pin the
// exact answer down, and nothing here changes what is being checked.
export struct BruteForceSettings
{
  double spacing{0.15};                // Å, of the grid the void is flooded on
  std::size_t samplesPerAtom{20000};   // directions thrown at each atom for the surface
  std::size_t volumePoints{4000000};   // points thrown at the cell for the void
  std::size_t subdivisions{1};         // panels per smooth piece, for the exact route it is held against

  std::size_t creaseSteps{256};    // steps round each crease, for the saddle sweep
  std::size_t cornerSamples{4096}; // directions drawn at each corner, for the concave patches

  bool skipSolventExcluded{false};  // the slowest of them, and the only one that is a quadrature
};

// One number worked out twice and the difference between the two.
export struct BruteForceCheck
{
  std::string property;
  std::string units;

  double exact{0.0};
  double bruteForce{0.0};

  // How far apart the two are allowed to be, and where that came from. A brute force that samples has a
  // standard error and a brute force that counts voxels has a spacing, so the tolerance is not a wish about
  // how close the answer ought to be: it is what this particular check can resolve.
  double tolerance{0.0};
  std::string basis;

  bool applicable{true};

  double difference() const { return bruteForce - exact; }
  bool agrees() const { return !applicable || std::abs(difference()) <= tolerance; }
};

// The exact routes and the brute force, run on the same structure and set against each other.
//
// What this is for is the one thing a test of the exact routes against themselves cannot do. Those routes
// share a decomposition of the boundary into patches, arcs and vertices, and a classifier that says which
// pore each patch faces; a mistake in either is invisible to any two of them compared with each other,
// because both inherit it. The checks here are built on none of that -- a grid, a flood and a lot of random
// points -- so they are wrong in unrelated ways, and where they agree the agreement means something.
//
// They are correspondingly blunt. A sampled area has a standard error of parts in a thousand where the
// exact one has parts in a trillion. So this says whether the exact answer is right to three or four
// digits, and nothing about the digits after that.
export struct BruteForceValidation
{
  std::vector<BruteForceCheck> checks;

  BruteForceDiameters diameters;
  BruteForceSurfaceArea surfaceArea;
  BruteForceSolventExcluded solventExcluded;
  BruteForcePoreVolume poreVolume;
  BruteForceBlockingPockets blockingPockets;
  BlockingAudit blockingAudit;

  // What the exact route said about its own vertices, so that a disagreement over the concave patches can
  // be read against how many corners each side thought there were.
  std::size_t exactVertices{0};
  std::size_t exactClippedVertices{0};
  std::size_t exactDegenerateVertices{0};
  std::size_t exactVanishedVertices{0};
  std::size_t exactDiscardedCorners{0};

  // Necks the flood on the probe's own grid stepped over and a straight line then proved passable. The one
  // for the void grid is carried by `poreVolume`; this is the other grid, the one the accessible and sealed
  // surface are split on.
  std::size_t surfaceNecksProved{0};
  std::size_t surfaceNecksTried{0};

  std::size_t numberOfDisagreements() const;

  double seconds{0.0};

  // Runs both routes and writes `{framework}.brute-force.txt`. `surfaceProbe` is the probe the surface area
  // and the blocking spheres are measured with and `voidProbe` the one the void fraction is; the pore
  // diameters take no probe, being about the void itself.
  void run(const PairInteractions &interactions, const Crystal &framework, const std::string &surfaceProbe,
           const std::string &voidProbe, const BruteForceSettings &settings);
};
