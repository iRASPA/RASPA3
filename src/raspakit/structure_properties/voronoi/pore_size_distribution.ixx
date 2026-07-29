module;

export module voronoi_pore_size_distribution;

import std;

import framework;
import forcefield;
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

  // `probePseudoAtom` is the probe the accessible curve is reported for, beside the curve of the whole void.
  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins,
           std::size_t subdivisions = 1);
};

// The report both diagrams write, the two differing only in the name in it.
export void writePoreSizeDistribution(const Framework& framework, const std::string& diagramName,
                                      const std::string& probePseudoAtom, const PoreSizeDistributionCurve& curve);
