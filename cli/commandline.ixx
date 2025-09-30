module;

#ifdef USE_LEGACY_HEADERS
#include <bitset>
#include <complex>
#include <cstddef>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <locale>
#include <mutex>
#include <optional>
#include <ranges>
#include <semaphore>
#include <span>
#include <string_view>
#include <vector>
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
  PSD_BV = 9, // Pore Size Distribution using Ban, Vlugt method
  Last = 10
};

ForceField defaultForceFieldZeolite(double rc = 12.0, bool shifted = false, bool tailCorrections = false,
                                    bool useEwald = false);
ForceField defaultForceFieldMOF(double rc = 12.0, bool shifted = false, bool tailCorrections = false,
                                bool useEwald = false);

void run(int argc, char* argv[]);
}  // namespace CommandLine
