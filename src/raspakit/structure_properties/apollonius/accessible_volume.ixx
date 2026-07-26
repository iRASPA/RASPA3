module;

export module apollonius_accessible_volume;

import std;

import framework;
import forcefield;

// Accessible and inaccessible void volume, sampled against the Apollonius diagram.
//
// Points are drawn uniformly through the cell and asked the same three-way question as on the radical
// network: solid, void a probe can reach, or void closed off in a pocket. The accessible fraction is
// the probe-accessible void fraction, zeo++'s -volpo.
export struct ApolloniusAccessibleVolume
{
  double accessibleVolumeFraction{0.0};
  double inaccessibleVolumeFraction{0.0};
  double accessibleVolume{0.0};    // Å³
  double inaccessibleVolume{0.0};  // Å³

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfSamples = std::nullopt);
};
