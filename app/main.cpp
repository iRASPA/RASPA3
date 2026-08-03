#include <locale.h>
import std;

import graceful_shutdown;
import input_reader;
import run_simulation;
import opencl;
#ifdef BUILD_LIBTORCH
import libtorch_test;
#endif

int main(int argc, char* argv[])
{
  using namespace std::literals;

  setlocale(LC_ALL, "en-US");

  // SIGTERM/SIGINT/SIGUSR1 request a final binary restart file at the next cycle boundary
  GracefulShutdown::install();

  OpenCL::initialize();

#ifdef BUILD_LIBTORCH
  test_libtorch();
#endif

  std::vector<std::string_view> args(argv, argv + argc);

  for (auto it = args.begin(); it != args.end(); ++it)
  {
    if (*it == "--help"sv || *it == "-h"sv)
    {
      std::cout << "RASPA is a software package for simulating adsorption and\n"
                   "diffusion of molecules in flexible nanoporous materials.\n"
                   "The code implements the latest state-of-the-art algorithms\n"
                   "for Molecular Dynamics and Monte Carlo in various ensembles\n"
                   "including symplectic/measure-preserving integrators, Ewald\n"
                   "summation, Configurational-Bias Monte Carlo, Continuous\n"
                   "Fractional Component Monte Carlo, and Reactive Monte Carlo.\n";
      return 0;
    }
    else if (*it == "--opencl"sv)
    {
      std::cout << OpenCL::printBestOpenCLDevice();
      return 0;
    }
  }

  try
  {
    InputReader inputReader("simulation.json");

    runSimulation(inputReader);
  }
  catch (std::exception const& e)
  {
    std::cerr << e.what();
    std::exit(-1);
  }
  catch (...)
  {
    std::cerr << "Exception caught" << std::endl;
    std::exit(-1);
  }
}
