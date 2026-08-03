module;

export module mc_pore_size_distribution;

import std;

import sampled_structure;

// The pore-size distribution of Gelb and Gubbins, sampled.
//
// The size at a point of the void is the diameter of the largest sphere that fits in the void and still
// covers that point. Here it is estimated: a point is drawn in the void, `numberOfInnerSteps` centres are
// drawn after it, and the largest sphere among those that fits and reaches the point gives its size. The
// whole is repeated `numberOfIterations` times and the sizes gathered into bins.
//
// Both counts are approximations of their own. Too few centres and the sphere that covers the point is
// never drawn, so every size comes out too small; too few points and the histogram is noise. Neither error
// shows up as scatter in the curve, which is what makes the exact route worth having beside it.
export struct MC_PoreSizeDistribution
{
  std::size_t numberOfBins;
  std::vector<double> histogram;
  std::vector<double> histogram_cummulative;

  double seconds{0.0};

  MC_PoreSizeDistribution(std::size_t numberOfBins)
      : numberOfBins(numberOfBins), histogram(numberOfBins), histogram_cummulative(numberOfBins) {};

  void run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
           std::optional<std::size_t> numberOfInnerSteps, std::optional<double> maximumRange);
};
