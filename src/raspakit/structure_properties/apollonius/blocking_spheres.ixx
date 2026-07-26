module;

export module apollonius_blocking_spheres;

import std;

import framework;
import forcefield;
import voronoi_blocking_spheres;

// Blocking spheres for the pockets a probe cannot reach, found against the Apollonius diagram.
//
// The inaccessible void is sampled and covered greedily as on the radical network; what the diagram
// decides is which void is inaccessible in the first place. The result is written in the RASPA
// `.block` format that a simulation reads back.
export struct ApolloniusBlockingSpheres
{
  std::vector<BlockingSphere> spheres;

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfSamples = std::nullopt);
};
