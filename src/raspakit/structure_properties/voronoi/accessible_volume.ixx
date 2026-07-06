module;

export module voronoi_accessible_volume;

import std;

import framework;
import forcefield;

// Monte-Carlo accessible volume split into accessible and inaccessible (pocket) void,
// using the Voronoi accessibility classifier. Points are sampled uniformly in the unit
// cell; points inside inflated atoms are solid, the rest are accessible or inaccessible
// void.
export struct VoronoiAccessibleVolume
{
  double accessibleVolumeFraction{0.0};
  double inaccessibleVolumeFraction{0.0};
  double accessibleVolume{0.0};    // Å³
  double inaccessibleVolume{0.0};  // Å³

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfSamples = std::nullopt);
};
