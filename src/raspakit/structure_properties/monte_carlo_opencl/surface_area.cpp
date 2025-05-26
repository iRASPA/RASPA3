module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#include <optional>
#include <numbers>
#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <exception>
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
  #include <OpenCL/cl.h>
#elif _WIN32
  #include <CL/cl.h>
#else
  #include <CL/opencl.h>
#endif
#endif

module mc_opencl_surface_area;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <print>;
import <string>;
import <vector>;
import <optional>;
import <numbers>;
import <limits>;
import <algorithm>;
import <iostream>;
import <fstream>;
import <exception>;
#endif

import double3;
import atom;
import randomnumbers;
import framework;
import forcefield;
import units;

void MC_OpenCL_SurfaceArea::run([[maybe_unused]]const ForceField &forceField, [[maybe_unused]]const Framework &framework, [[maybe_unused]]double well_depth_factor, 
                             [[maybe_unused]]std::string probe_pseudo_atom, [[maybe_unused]]size_t number_of_iterations) const
{
  RandomNumber random{std::nullopt};

}
