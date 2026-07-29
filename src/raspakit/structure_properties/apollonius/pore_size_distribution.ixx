module;

export module apollonius_pore_size_distribution;

import std;

import framework;
import forcefield;
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

  // `probePseudoAtom` is the probe the accessible curve is reported for, beside the curve of the whole void.
  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins,
           std::size_t subdivisions = 1);
};
