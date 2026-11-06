module;

export module apollonius_pore_size_distribution;

import std;

import crystal;
import pair_interactions;
import blocking_spheres;
import exact_pore_size_distribution;

// The same pore-size distribution, with the pores taken from the Apollonius diagram.
//
// The curve is identical to the one the radical network gives, and has to be: it is the volume enclosed by the
// solvent excluded surface, which is a property of the atoms alone. Which diagram is used changes only how the
// distribution is divided between the void the probe can reach and the void it cannot, and on a framework it
// does not change that either. The choice costs what building the diagram costs and nothing else.
export struct ApolloniusPoreSizeDistribution
{
  PoreSizeDistributionCurve curve;

  // The spheres the default curve was measured on: the framework's atoms at their own radii, and one ball over
  // each pocket sealed to the probe, which is how those pockets are left out. Empty where the curve was asked
  // for with a floor at the probe, which is measured on the atoms alone. It is kept so that the curve can be
  // checked by something that takes the same spheres and works the pore sizes out another way.
  std::vector<BlockingSphere> sweptSpheres;

  // Whether those are more than the atoms: whether any pocket got a ball. Their number does not say, a ball
  // that swallows an atom taking one off the list as it adds itself.
  bool pocketsFilled{false};

  // `probePseudoAtom` is the accessibility / blocking probe. `floorRadius` is the diameter floor of the
  // primary curve (default 0 = hybrid: block sealed pockets, keep bare pore sizes on the rest).
  // With `peaksOnly`, only room-size spikes are evaluated (see `exactPoreSizePeaks`).
  void run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom,
           std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins,
           std::size_t subdivisions = 1, double floorRadius = 0.0, bool peaksOnly = false);
};
