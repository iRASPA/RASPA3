module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#endif

export module mc_opencl_void_fraction;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <print>;
import <string>;
import <vector>;
#endif

import framework;
import forcefield;

export struct MC_OpenCL_VoidFraction
{
  std::vector<double> data;

  MC_OpenCL_VoidFraction()
  {
  };

  void run(const ForceField &forceField, const Framework &framework, size_t number_of_iterations);

};
