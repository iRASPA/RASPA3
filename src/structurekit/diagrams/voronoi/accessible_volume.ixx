module;

export module voronoi_accessible_volume;

import std;

import crystal;
import pair_interactions;
import pore_accessibility;
import exact_void_split;

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
export VolumeSample sampleAccessibleVolume(const PoreAccessibility& accessibility, std::size_t numberOfSamples);

// Void volume against the radical network: how much of it there is, and how much of that a molecule can
// reach.
//
// As with the surface area, how much void there is belongs to the union of the inflated atoms and not to
// any diagram, so `voidFraction` is measured and is the same number the Apollonius analysis reports. The
// division between channels and pockets is measured too, from the arcs of the surface that face each
// sealed pore; see exact_void_split. Which pores those are is the radical network's answer, and it reaches
// further into pockets than a probe can.
export struct VoronoiAccessibleVolume
{
  enum class Method
  {
    Exact,   // the void fraction and its division between pores both measured
    Sampled  // both sampled, as zeo++
  };

  double voidFraction{0.0};  // the whole of it, solid subtracted
  double accessibleVolumeFraction{0.0};
  double inaccessibleVolumeFraction{0.0};
  double voidVolume{0.0};          // Å³
  double accessibleVolume{0.0};    // Å³
  double inaccessibleVolume{0.0};  // Å³

  // Whether the division above was measured or sampled. The exact method falls back on the sampling where
  // the surface leaves a pore undecided or a pocket's boundary does not close, so that the two are not
  // reported as if they were the same thing.
  bool splitMeasured{false};
  std::size_t numberOfPockets{0};
  double closureDefect{0.0};
  std::string splitRejection;  // why the measured division was not used, when it was not

  // How the division was arrived at: the connected surfaces of the boundary, how many of them were asked
  // which side the void is on with a free ball to prove it by, and how many bounded ones the percolation
  // analysis disagrees with the enclosed volume's own sign about. See exact_void_split.
  std::size_t numberOfSurfaces{0};
  std::size_t provedSurfaces{0};
  std::size_t provedPockets{0};
  std::size_t numberOfEnclosedSolids{0};
  std::size_t signDisagreements{0};
  double signDisagreementVolume{0.0};  // Å³

  // Where each pocket is and how large, from the surface around it. The boundary that gives a pocket its
  // volume gives the first moment of the region as well, so the centre costs one more sum per arc.
  std::vector<PocketGeometry> pockets;

  void run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom,
           Method method = Method::Exact, std::optional<std::size_t> numberOfSamples = std::nullopt,
           std::optional<std::size_t> subdivisions = std::nullopt);
};
