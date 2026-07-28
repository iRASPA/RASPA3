module;

export module apollonius_blocking_spheres;

import std;

import framework;
import forcefield;
import exact_void_split;
import voronoi_blocking_spheres;

// Blocking spheres for the pockets a probe cannot reach, found against the Apollonius diagram.
//
// Where each pocket is, how far it reaches and how near the channels come are all read off the surface, exactly
// as on the radical network and by the same argument; see voronoi_blocking_spheres. What the diagram is for is
// the one question the surface leaves open, which is whether the void around a cluster of atoms standing in it
// is sealed, and on a framework there is no such cluster. So the two networks give the same spheres here, and
// the choice between them costs only what building the diagram costs.
export struct ApolloniusBlockingSpheres
{
  std::vector<BlockingSphere> spheres;

  // Whether the spheres came from the surface or from sampling, and why not from the surface when they did not.
  bool measured{false};
  std::string fallbackReason;

  // The pockets the spheres stand for, when they were measured.
  std::vector<PocketGeometry> pockets;

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfSamples = std::nullopt);
};
