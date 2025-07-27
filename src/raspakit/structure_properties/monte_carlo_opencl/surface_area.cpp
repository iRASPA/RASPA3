module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <print>
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

module mc_opencl_surface_area;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import atom;
import randomnumbers;
import framework;
import forcefield;
import units;

void MC_OpenCL_SurfaceArea::run([[maybe_unused]] const ForceField &forceField,
                                [[maybe_unused]] const Framework &framework, [[maybe_unused]] double well_depth_factor,
                                [[maybe_unused]] std::string probe_pseudo_atom,
                                [[maybe_unused]] std::size_t number_of_iterations) const
{
  RandomNumber random{std::nullopt};
}
