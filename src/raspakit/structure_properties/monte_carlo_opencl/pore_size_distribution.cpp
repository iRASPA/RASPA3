module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <print>
#include <random>
#include <string>
#include <vector>
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif
#endif

module mc_opencl_pore_size_distribution;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import randomnumbers;
import framework;
import forcefield;
import atom;

void MC_OpenCL_PoreSizeDistribution::run([[maybe_unused]] const ForceField &forceField,
                                         [[maybe_unused]] const Framework &framework,
                                         [[maybe_unused]] double well_depth_factor,
                                         [[maybe_unused]] std::size_t number_of_iterations)
{
  RandomNumber random{std::nullopt};
}
