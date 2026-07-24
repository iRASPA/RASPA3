#include <locale.h>
import std;

import archive;
import graceful_shutdown;
import threadpool;
import input_reader;
import monte_carlo;
import monte_carlo_transition_matrix;
import molecular_dynamics;
import minimization;
import thermodynamic_integration;
import parallel_thermodynamic_integration;
import parallel_tempering;
import hyper_parallel_tempering;
import reweighted_histogram;
import parallel_tmmc;
//import breakthrough;
//import breakthrough_simulation;
//import mixture_prediction_simulation;
//import isotherm_fitting_simulation;
//import multi_site_isotherm;
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

    auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type>::instance();
    pool.init(inputReader.numberOfThreads, inputReader.threadingType);

    switch (inputReader.simulationType)
    {
      case InputReader::SimulationType::MonteCarlo:
      {
        MonteCarlo mc(inputReader);
        if (inputReader.restartFromBinary)
        {
          readBinaryRestartFile(mc, inputReader.restartFromBinaryFileName);
          // the output files are opened in append mode (the header is already in the file), and
          // the interpolation grids are not stored in the restart file and must be rebuilt
          mc.createOutputFiles();
          mc.createInterpolationGrids();
        }

        mc.run();
        break;
      }
      case InputReader::SimulationType::MonteCarloTransitionMatrix:
      {
        MonteCarloTransitionMatrix mc(inputReader);
        if (inputReader.restartFromBinary)
        {
          readBinaryRestartFile(mc, inputReader.restartFromBinaryFileName);
          // the resumed stage jumps past the header block, so the output streams must exist
          mc.createOutputFiles();
        }
        mc.run();
        break;
      }
      case InputReader::SimulationType::MolecularDynamics:
      {
        MolecularDynamics md(inputReader);
        if (inputReader.restartFromBinary)
        {
          readBinaryRestartFile(md, inputReader.restartFromBinaryFileName);
          md.createOutputFiles();
          // the grids are not stored in the restart file and must be rebuilt
          md.createInterpolationGrids();
        }

        md.run();
        break;
      }
      case InputReader::SimulationType::Minimization:
      {
        Minimization minimization(inputReader);
        if (inputReader.restartFromBinary)
        {
          readBinaryRestartFile(minimization, inputReader.restartFromBinaryFileName);
        }
        minimization.run();
        break;
      }
      case InputReader::SimulationType::ThermodynamicIntegration:
      {
        ThermodynamicIntegration ti(inputReader);
        ti.run();
        break;
      }
      case InputReader::SimulationType::ParallelThermodynamicIntegration:
      {
        ParallelThermodynamicIntegration parallel_ti(inputReader);
        if (inputReader.restartFromBinary)
        {
          readBinaryRestartFile(parallel_ti, inputReader.restartFromBinaryFileName);
        }
        parallel_ti.run();
        break;
      }
      case InputReader::SimulationType::ParallelTempering:
      {
        ParallelTempering parallel_tempering(inputReader);
        if (inputReader.restartFromBinary)
        {
          readBinaryRestartFile(parallel_tempering, inputReader.restartFromBinaryFileName);
        }
        parallel_tempering.run();
        break;
      }
      case InputReader::SimulationType::HyperParallelTempering:
      {
        HyperParallelTempering hyper_parallel_tempering(inputReader);
        if (inputReader.restartFromBinary)
        {
          readBinaryRestartFile(hyper_parallel_tempering, inputReader.restartFromBinaryFileName);
        }
        hyper_parallel_tempering.run();
        break;
      }
      case InputReader::SimulationType::ReweightedHistogram:
      {
        ReweightedHistogram reweighted_histogram(inputReader);
        if (inputReader.restartFromBinary)
        {
          readBinaryRestartFile(reweighted_histogram, inputReader.restartFromBinaryFileName);
        }
        reweighted_histogram.run();
        break;
      }
      case InputReader::SimulationType::ParallelTMMC:
      {
        ParallelTMMC parallel_tmmc(inputReader);
        if (inputReader.restartFromBinary)
        {
          readBinaryRestartFile(parallel_tmmc, inputReader.restartFromBinaryFileName);
        }
        parallel_tmmc.run();
        break;
      }
      default:
        break;
    }
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
