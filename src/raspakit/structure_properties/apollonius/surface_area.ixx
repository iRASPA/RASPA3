module;

export module apollonius_surface_area;

import std;

import framework;
import forcefield;

// The measured route and the sampled one live in modules of their own, neither of which knows about any
// diagram; what this adds is the Apollonius classifier they are run against.
export import exact_surface_patches;
import exact_solvent_excluded;

// Accessible and inaccessible surface area against the Apollonius diagram, either measured or sampled.
//
// Both routes split the boundary of the union of the probe-inflated atoms the same way, into the part
// reachable from outside the crystal and the part sealed in pockets, and both take the same classifier
// from the Apollonius diagram to decide which is which. What differs is how the area itself is arrived
// at: the sampled route is zeo++'s, points thrown on each inflated sphere and those buried in another
// atom thrown away, whose error is statistical and falls off as one over the square root of their
// number; the measured route integrates each patch and does not sample the sphere at all.
export struct ApolloniusSurfaceArea
{
  enum class Method
  {
    Exact,   // the latitude integral, accurate to round-off
    Sampled  // points on the inflated spheres, as zeo++
  };

  double accessibleSurfaceArea{0.0};    // Å²
  double inaccessibleSurfaceArea{0.0};  // Å²
  double undecidedSurfaceArea{0.0};     // Å², only ever nonzero for the exact method

  // The wall the probe touches, rather than the sheet its centre traces: the excluded surface at the same probe,
  // by the kind of patch and by the side it faces. Filled by the exact method alone.
  SolventExcludedGeometry excludedSurface;

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           Method method = Method::Exact, std::optional<std::size_t> samplesPerAtom = std::nullopt,
           std::optional<std::size_t> subdivisions = std::nullopt);
};
