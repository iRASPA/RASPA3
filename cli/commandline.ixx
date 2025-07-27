module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <span>
#include <deque>
#include <optional>
#include <semaphore>
#include <mutex>
#include <complex>
#include <locale>
#include <ranges>
#include <string_view>
#include <filesystem>
#include <bitset>
#endif

export module commandline;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import threadpool;
import input_reader;
import forcefield;
import monte_carlo;
import monte_carlo_transition_matrix;
import molecular_dynamics;
import breakthrough;
import breakthrough_simulation;
import mixture_prediction_simulation;
import isotherm_fitting_simulation;
import multi_site_isotherm;
import opencl;
import getopt;
#ifdef BUILD_LIBTORCH
import libtorch_test;
#endif

export namespace CommandLine
{
  enum State : std::uint8_t
  {
    None = 0,
    Help = 1,
    OpenCL = 2,
    Input = 3,
    SurfaceArea = 4,
    VoidFraction = 5,
    TessellationComputation = 6,
    PSD = 7,
    EnergyGrid = 8,
    Last = 9
  };

  ForceField defaultForceFieldZeolite(double rc = 12.0, bool shifted = false, bool tailCorrections = false, bool useEwald = false);
  ForceField defaultForceFieldMOF(double rc = 12.0, bool shifted = false, bool tailCorrections = false, bool useEwald = false);

  void run(int argc, char* argv[]);
}
