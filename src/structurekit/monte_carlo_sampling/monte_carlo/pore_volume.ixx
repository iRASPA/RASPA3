module;

export module mc_pore_volume;

import std;

import sampled_structure;
import sampled_roadmap;

// How much room there is in a crystal, and how much of it a probe let in from outside can reach, by
// throwing points at it.
//
// The fraction is the oldest Monte Carlo measurement there is and needs no roadmap: throw a point, ask
// whether it is clear of every atom, count. What it cannot do on its own is say how much of that room is
// shut away in pockets, because a single point carries no information about what it is connected to. The
// other routes here answer that from a diagram of the void, exactly where the diagram closes and by falling
// back on sampling where it does not. This one answers it from the roadmap: a point is reachable when the
// component of the roadmap it belongs to runs off across the crystal, and shut away when its component
// closes on itself.
//
// The two readings err in opposite directions and that is worth keeping in mind. The void fraction is
// unbiased -- it is a binomial proportion, with the standard error of one, and neither route to it can be
// nearer the truth than that. The split is not: hops the sample failed to find can break a channel into
// pieces that look like pockets, but no hop is ever invented, so the accessible share is a lower bound and
// the inaccessible share an upper one. A run that reports many more pockets than the structure has is a run
// whose sample was too thin, and the count is printed for exactly that reason.
export struct MC_PoreVolume
{
  // The share of the cell a probe's centre can occupy at all, whether or not it can get there. This is a
  // binomial proportion over the points thrown, and `voidFractionError` is its standard error.
  double voidFraction{0.0};
  double voidFractionError{0.0};

  double accessibleVolumeFraction{0.0};    // of it, the part connected to a channel
  double inaccessibleVolumeFraction{0.0};  // the part shut in a pocket the sample resolved
  double unsettledVolumeFraction{0.0};     // and what the sample could not place either way

  double voidVolume{0.0};          // Å³, per unit cell
  double accessibleVolume{0.0};    // Å³
  double inaccessibleVolume{0.0};  // Å³

  // The same three as a simulation would quote them, which is per gram of framework rather than per cell.
  double gravimetricVoidVolume{0.0};          // cm³/g
  double gravimetricAccessibleVolume{0.0};    // cm³/g
  double gravimetricInaccessibleVolume{0.0};  // cm³/g

  double density{0.0};  // kg/m³

  std::size_t numberOfChannels{0};  // components of the roadmap that run off across the crystal
  std::size_t numberOfPockets{0};   // and those that close on themselves

  SampledRoadmap roadmap;
  double seconds{0.0};

  void run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
           std::optional<std::size_t> numberOfInnerSteps);

  void run(const SampledStructure &structure, const SampledRoadmap &roadmap);
};
