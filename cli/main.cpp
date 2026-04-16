#include <locale.h>

import std;

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
//import mixture_prediction_simulation;
//import isotherm_fitting_simulation;
//import multi_site_isotherm;
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
