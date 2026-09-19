#include <locale.h>
import std;

import graceful_shutdown;
import input_reader;
import run_simulation;
import opencl;
import commandline;
#ifdef BUILD_LIBTORCH
import libtorch_test;
#endif

int main(int argc, char* argv[])
{
  using namespace std::literals;

  setlocale(LC_ALL, "en-US");

  OpenCL::initialize();

#ifdef BUILD_LIBTORCH
  test_libtorch();
#endif

  // Any non-empty argv besides a lone --opencl routes to the analysis CLI.
  // Bare `raspa3` (and `raspa3 --opencl`) keep the simulation.json driver.
  const bool analysis_cli =
      argc > 1 && !(argc == 2 && std::string_view(argv[1]) == "--opencl"sv);

  if (analysis_cli)
  {
    try
    {
      CommandLine::run(argc, argv);
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
    return 0;
  }

  // SIGTERM/SIGINT/SIGUSR1 request a final binary restart file at the next cycle boundary
  GracefulShutdown::install();

  std::vector<std::string_view> args(argv, argv + argc);

  for (auto it = args.begin(); it != args.end(); ++it)
  {
    if (*it == "--opencl"sv)
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
