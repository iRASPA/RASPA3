module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#endif

export module mc_opencl_surface_area;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <print>;
import <string>;
import <vector>;
#endif

import framework;
import forcefield;

export struct MC_OpenCL_SurfaceArea
{
  std::vector<double> data;

  MC_OpenCL_SurfaceArea()
  {
  };

  void run(const ForceField &forceField, const Framework &framework, double probe_distance, 
           std::string probe_pseudo_atom, size_t number_of_iterations) const;

};

