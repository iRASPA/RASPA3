module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#endif

export module mc_opencl_pore_size_distribution;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <print>;
import <string>;
import <vector>;
#endif

import framework;
import forcefield;

export struct MC_OpenCL_PoreSizeDistribution
{
  std::vector<double> data;

  MC_OpenCL_PoreSizeDistribution(size_t numberOfBins):
    data(numberOfBins)
  {
  };

  void run(const ForceField &forceField, const Framework &framework, double well_depth_factor, size_t number_of_iterations);
};
