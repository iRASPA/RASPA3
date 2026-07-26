module;

export module voronoi_blocking_spheres;

import std;

import double3;
import simulationbox;
import framework;
import forcefield;
import voronoi_accessibility;

// Generation of blocking spheres for inaccessible pockets. Points of the void that a probe cannot
// reach are sampled by Monte Carlo and grouped by pocket, and each pocket is covered by spheres taking
// the widest first. The result is written in the RASPA `.block` format: a count followed by one line
// per sphere holding the fractional centre and the radius.
//
// A sphere reaches as far as the pocket does and no further than the nearest position the probe can
// occupy from a channel. Both are measured between sample points, which are positions of the probe's
// *centre*, so neither needs a further allowance for the probe's size -- an earlier version added one
// to the first and subtracted one from the second, and since the two limits are typically only a few
// Å apart that left nothing in between: the radii collapsed and a pocket came back covered by hundreds
// of spheres too small to block anything, or was abandoned altogether.
//
// Below both limits sits a floor that holds regardless: a sphere of the clearance at its centre holds
// no atom, so it lies in the free space, and being connected it cannot cross the neck that made the
// pocket a pocket. Such a sphere is inside the pocket by construction.
export struct BlockingSphere
{
  double3 centerFractional;
  double radius;  // Å
};

// The sampling and covering itself, over whatever accessibility classifier is handed to it, so that
// the same spheres can be found from a network taken from the radical diagram or from the Apollonius
// diagram. The probe radius is not a parameter: it is already in the inflated radii the classifier
// carries, which is what the clearance is measured against.
export std::vector<BlockingSphere> computeBlockingSpheres(const VoronoiAccessibility& accessibility,
                                                          std::size_t numberOfSamples);

export struct VoronoiBlockingSpheres
{
  std::vector<BlockingSphere> spheres;

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfSamples = std::nullopt);
};
