#ifdef USE_LEGACY_HEADERS
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
#endif

#ifndef USE_LEGACY_HEADERS
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

int main(int argc, char* argv[])
{  
  setlocale(LC_ALL, "en-US");

  std::vector<std::string> cmdLineArgs(argv, argv+argc);

  for(auto& arg : cmdLineArgs)
  {
    if(arg == "--help" || arg == "-help")
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
  }

  try
  {
    InputReader inputReader("simulation.json");

    auto &pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type>::instance();
    pool.init(inputReader.numberOfThreads, inputReader.threadingType);

    switch (inputReader.simulationType)
    {
      case InputReader::SimulationType::MonteCarlo:
      {
        MonteCarlo mc(inputReader);
        if(inputReader.restartFromBinary)
        {
          std::ifstream ifile("restart_data.bin", std::ios::binary);
          if(!ifile.is_open())
          {
            throw std::runtime_error("Restart file doesn't exist..\n");
          }
          Archive<std::ifstream> archive(ifile);
          archive >> mc;
          mc.createOutputFiles();
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
        if(inputReader.restartFromBinary)
        {
          std::ifstream ifile("restart_data.bin", std::ios::binary);
          if(!ifile.is_open())
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
    exit(-1);
  }
}
