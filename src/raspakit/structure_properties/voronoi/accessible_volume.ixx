module;

export module voronoi_accessible_volume;

import std;

import framework;
import forcefield;
import voronoi_accessibility;

// Monte-Carlo accessible volume split into accessible and inaccessible (pocket) void,
// using the Voronoi accessibility classifier. Points are sampled uniformly in the unit
// cell; points inside inflated atoms are solid, the rest are accessible or inaccessible
// void.
export struct VolumeSample
{
  double accessibleFraction{0.0};
  double inaccessibleFraction{0.0};
};

// The sampling itself, over whatever accessibility classifier is handed to it, so that the same
// estimate can be made of a network taken from the radical diagram or from the Apollonius diagram.
export VolumeSample sampleAccessibleVolume(const VoronoiAccessibility& accessibility, std::size_t numberOfSamples);

export struct VoronoiAccessibleVolume
{
  double accessibleVolumeFraction{0.0};
  double inaccessibleVolumeFraction{0.0};
  double accessibleVolume{0.0};    // Å³
  double inaccessibleVolume{0.0};  // Å³

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfSamples = std::nullopt);
};
