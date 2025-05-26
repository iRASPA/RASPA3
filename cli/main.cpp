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

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <exception>;
import <iostream>;
import <fstream>;
import <vector>;
import <span>;
import <deque>;
import <optional>;
import <semaphore>;
import <mutex>;
import <complex>;
import <locale>;
import <string_view>;
#endif

import archive;
import threadpool;
import input_reader;
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
  catch (const std::exception &e)
  {
    std::cerr << e.what();
    exit(-1);
  }
}
