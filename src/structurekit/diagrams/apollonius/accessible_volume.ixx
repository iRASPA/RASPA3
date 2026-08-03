module;

export module apollonius_accessible_volume;

import std;

import crystal;
import pair_interactions;
import exact_void_split;

// Void volume against the Apollonius diagram: how much of it there is, and how much of that a molecule
// can reach.
//
// Neither question is answered by sampling. The whole void volume is the cell less the volume of the
// union of the probe-inflated atoms, which is a finite sum of closed-form pieces; see
// apollonius_union_volume. The division between the channels and the pockets comes from the surface: a
// sealed pore is bounded and its boundary is exactly the arcs of the surface that face it, so the
// divergence theorem gives its volume from those arcs, and the channels take what is left of the void.
// See exact_void_split.
//
// The division used to be sampled, on the argument that the pieces of a volume are not pieces of one pore:
// a piece of surface is an arc, which is connected, so one question settles all of it, whereas the void
// part of one cell is a three-dimensional region that can be disconnected, an atom in the wall between two
// channels having a cell that reaches into both. That argument is sound and it is about cells. Going by
// the boundary instead needs no cell, and so no connected pieces of one either.
//
// What remains sampled is nothing, unless the surface leaves a pore undecided or a pocket's boundary fails
// to close, in which case the division falls back on the sampling and says so. The accessible fraction is
// the one zeo++ reports as -volpo.
export struct ApolloniusAccessibleVolume
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

  // Whether the division above was measured or sampled, and the two numbers that decide it: how many
  // sealed pores carried any surface, and how far the worst of their boundaries was from closing.
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
