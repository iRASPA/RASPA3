module;

export module voronoi_blocking_spheres;

import std;

import double3;
import framework;
import forcefield;

// Generation of blocking spheres for inaccessible pockets. Inaccessible void points are
// sampled by Monte Carlo, grouped by pocket, and each pocket is covered greedily by
// spheres. The result is written in the RASPA `.block` format: a count followed by one
// line per sphere holding the fractional centre and the radius.
export struct BlockingSphere
{
  double3 centerFractional;
  double radius;  // Å
};

export struct VoronoiBlockingSpheres
{
  std::vector<BlockingSphere> spheres;

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfSamples = std::nullopt);
};
