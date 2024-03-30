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

import archive;
import threadpool;
import input_reader;
import monte_carlo;
import monte_carlo_transition_matrix;
import breakthrough;
import breakthrough_simulation;
import mixture_prediction_simulation;
import isotherm_fitting_simulation;
import multi_site_isotherm;

//unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

int main()
{  
  try
  {
    InputReader inputReader("simulation.input");

    ThreadPool::createPool(inputReader.numberOfThreads, inputReader.threadingType);

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
