module;

export module mc_pore_diameters;

import std;

import pore_diameters;
import sampled_structure;
import sampled_roadmap;

// Di, Df and Dif by sampling: the void mapped by throwing points at it, and the three diameters read off
// the map.
//
// The other routes to these numbers build a diagram of the void and read the diameters from its vertices
// and arcs. That is exact, and it is the reason those routes exist, but it is also a construction whose
// cost and whose failure modes have nothing to do with the question being asked. This one asks the question
// directly, off the roadmap of `SampledRoadmap`, by the same arithmetic that reads the diameters off a
// diagram -- `PoreDiameters::compute` does not know which of the two it was handed.
//
// Everything it reports is therefore a lower bound, and this is the property worth understanding about it.
// A sphere the roadmap says fits, fits, because a point of the void was found with that much room around
// it. A path the roadmap says is open, is open, because it is a chain of straight hops each of which was
// checked. What sampling can miss is a deeper pocket that no point landed in and a wider way round that no
// chain of hops happened to take, and both of those can only make the answer smaller. So the numbers rise
// towards the truth as the sample grows, from below and never past it, and two runs at different sample
// sizes bracket how far there is left to go. The exact routes have no such handle on their own error.
//
// Di comes from the walks uphill and Df from the ways between the pockets those walks end at, so what Df
// converges by is how completely the pockets have been found, not how finely the void has been covered.
export struct MC_PoreDiameters
{
  PoreDiameters result;  // Di, Df, Dif [Å]

  // The roadmap the diameters were read from, kept for the analyses that follow and so that a caller can
  // see what the answer rests on.
  SampledRoadmap roadmap;

  bool percolates{false};  // whether any sphere at all can cross the crystal
  double seconds{0.0};

  void run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
           std::optional<std::size_t> numberOfInnerSteps);

  // The same, off a roadmap already built. Building it costs everything and reading it costs nothing, so a
  // caller wanting more than one property of a structure should build once and ask each in turn.
  void run(const SampledStructure &structure, const SampledRoadmap &roadmap);
};
