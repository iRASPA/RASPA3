module;

export module voronoi_pore_size_distribution;

import std;

import crystal;
import pair_interactions;
import exact_pore_size_distribution;

// The pore-size distribution of Gelb and Gubbins in closed form, with the pores taken from the radical
// (Voronoi) network.
//
// The curve itself asks nothing of any diagram: it is the volume of the union of the balls of a given radius
// that fit in the void, and that is the volume enclosed by the solvent excluded surface of the framework. What
// the network is for is the one thing the geometry leaves open, which is whether a pocket the probe cannot
// reach ought to count, and that only divides the curve rather than changing it.
export struct VoronoiPoreSizeDistribution
{
  PoreSizeDistributionCurve curve;

  // `probePseudoAtom` is the accessibility / blocking probe. `floorRadius` is the diameter floor of the
  // primary curve (default 0 = hybrid: block sealed pockets, keep bare pore sizes on the rest).
  // With `peaksOnly`, only room-size spikes are evaluated (see `exactPoreSizePeaks`).
  void run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom,
           std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins,
           std::size_t subdivisions = 1, double floorRadius = 0.0, bool peaksOnly = false);
};

// The report both diagrams write, the two differing only in the name in it.
export void writePoreSizeDistribution(const Crystal& framework, const std::string& diagramName,
                                      const std::string& probePseudoAtom, const PoreSizeDistributionCurve& curve);

// Peaks-only report: room diameters and the void fraction that sits at each, without continuous rows.
export void writePoreSizePeaks(const Crystal& framework, const std::string& diagramName,
                               const std::string& probePseudoAtom, const PoreSizeDistributionCurve& curve);

// Promote a bare PSD of the framework with He pockets filled by blocking spheres into the hybrid report
// shape. The open-framework whole-void curve is not computed: only the blocked network is physical for
// adsorption. `accessibilityRadius` is recorded as the blocking probe.
export PoreSizeDistributionCurve hybridFromBlockedCurve(PoreSizeDistributionCurve blocked,
                                                        double accessibilityRadius);
