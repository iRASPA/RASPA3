#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

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

#ifdef USE_STD_IMPORT
#include <locale.h>
#endif

#ifdef USE_STD_IMPORT
import std;
#endif

import archive;
import threadpool;
import input_reader;
import monte_carlo;
import monte_carlo_transition_matrix;

import archive;
import threadpool;
import input_reader;
import monte_carlo;
import monte_carlo_transition_matrix;
import molecular_dynamics;
import mixture_prediction_simulation;
import isotherm_fitting_simulation;
import multi_site_isotherm;
import opencl;
import getopt;
import commandline;
#ifdef BUILD_LIBTORCH
import libtorch_test;
#endif

int main(int argc, char* argv[])
{
  using namespace std::literals;

  setlocale(LC_ALL, "en-US");

  OpenCL::initialize();

  try
  {
    CommandLine::run(argc, argv);
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what();
    std::exit(-1);
  }
}
