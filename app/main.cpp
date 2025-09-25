#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <cstddef>
#include <deque>
#include <exception>
#include <fstream>
#include <iostream>
#include <locale>
#include <mutex>
#include <optional>
#include <semaphore>
#include <span>
#include <string_view>
#include <vector>
#endif

#ifndef USE_LEGACY_HEADERS
#include <locale.h>
import std;
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
          std::ifstream ifile(inputReader.restartFromBinaryFileName, std::ios::binary);
          if (!ifile.is_open())
          {
            throw std::runtime_error("Restart file doesn't exist..\n");
          }
          Archive<std::ifstream> archive(ifile);
          archive >> mc;
          mc.createOutputFiles();
          mc.writeOutputHeader();
          mc.createInterpolationGrids();
        }

        mc.run();
        break;
      }
      case InputReader::SimulationType::MonteCarloTransitionMatrix:
      {
        MonteCarloTransitionMatrix mc(inputReader);
        mc.run();
        break;
      }
      case InputReader::SimulationType::MolecularDynamics:
      {
        MolecularDynamics md(inputReader);
        if (inputReader.restartFromBinary)
        {
          std::ifstream ifile("restart_data.bin", std::ios::binary);
          if (!ifile.is_open())
          {
            throw std::runtime_error("Restart file doesn't exist..\n");
          }
          Archive<std::ifstream> archive(ifile);
          archive >> md;
          md.createOutputFiles();
        }

        md.run();
        break;
      }
      case InputReader::SimulationType::Breakthrough:
      {
        BreakthroughSimulation breakthrough(inputReader);
        breakthrough.run();
        break;
      }
      case InputReader::SimulationType::MixturePrediction:
      {
        MixturePredictionSimulation mixture(inputReader);
        mixture.run();
        break;
      }
      case InputReader::SimulationType::Fitting:
      {
        IsothermFittingSimulation fitting(inputReader);
        fitting.run();
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
}
