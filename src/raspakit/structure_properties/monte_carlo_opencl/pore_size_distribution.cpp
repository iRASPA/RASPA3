module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <print>
#include <string>
#include <vector>
#pragma push_macro("__SSE3__")
#undef __SSE3__
#include <random>
#pragma pop_macro("__SSE3__")
#define CL_TARGET_OPENCL_VERSION 120
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif
#endif

module mc_opencl_pore_size_distribution;

#ifdef USE_STD_IMPORT
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
