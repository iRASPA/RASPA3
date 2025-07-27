module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <print>
#include <string>
#include <tuple>
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

module mc_opencl_void_fraction;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import randomnumbers;
import atom;
import framework;
import forcefield;
import component;
import system;
import mc_moves_widom;

void MC_OpenCL_VoidFraction::run([[maybe_unused]] const ForceField &forceField,
                                 [[maybe_unused]] const Framework &framework,
                                 [[maybe_unused]] std::size_t number_of_iterations)
{
  RandomNumber random{std::nullopt};
}
