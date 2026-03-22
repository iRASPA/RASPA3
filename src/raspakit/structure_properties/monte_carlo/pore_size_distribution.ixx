module;

export module mc_pore_size_distribution;

import std;

import framework;
import forcefield;

export struct MC_PoreSizeDistribution
{
  std::size_t numberOfBins;
  std::vector<double> histogram;
  std::vector<double> histogram_cummulative;

  MC_PoreSizeDistribution(std::size_t numberOfBins) : numberOfBins(numberOfBins), histogram(numberOfBins), histogram_cummulative(numberOfBins) {};

  void run(const ForceField &forceField, const Framework &framework, double wellDepthFactor,
           std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps,
           std::optional<double> maximumRange);
};
