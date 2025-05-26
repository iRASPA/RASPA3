module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#include <optional>
#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <tuple>
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
import <cstddef>;
import <print>;
import <string>;
import <vector>;
import <optional>;
import <limits>;
import <algorithm>;
import <iostream>;
import <fstream>;
import <tuple>;
#endif

import double3;
import randomnumbers;
import atom;
import framework;
import forcefield;
import component;
import system;
import mc_moves_widom;



void MC_OpenCL_VoidFraction::run([[maybe_unused]]const ForceField &forceField, [[maybe_unused]]const Framework &framework, [[maybe_unused]]size_t number_of_iterations)
{
  RandomNumber random{std::nullopt};

}
