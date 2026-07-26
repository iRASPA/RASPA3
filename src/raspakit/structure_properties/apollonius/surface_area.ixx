module;

export module apollonius_surface_area;

import std;

import framework;
import forcefield;

// Accessible and inaccessible surface area, sampled against the Apollonius diagram.
//
// The estimate is the one made on the radical network: points on each probe-inflated sphere, those
// buried in another atom thrown away, the rest asked whether they sit in a channel or in a pocket.
// Only the classifier behind that question changes, so a difference in the answer is a difference in
// which pockets count as closed.
export struct ApolloniusSurfaceArea
{
  double accessibleSurfaceArea{0.0};    // Å²
  double inaccessibleSurfaceArea{0.0};  // Å²

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<std::size_t> samplesPerAtom = std::nullopt);
};
