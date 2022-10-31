
import <exception>;
import <iostream>;
import <fstream>;
import <vector>;
import <span>;
import <deque>;
import <optional>;
import <semaphore>;

import threadpool;
//import threading;
import input_reader;
import monte_carlo;
import breakthrough;
import breakthrough_simulation;
import mixture_prediction_simulation;
import isotherm_fitting_simulation;
import multi_site_isotherm;

int main()
{
  try
  {
    InputReader inputReader{};

    ThreadPool::createPool(inputReader.numberOfThreads, inputReader.threadingType);

    switch (inputReader.simulationType)
    {
      case InputReader::SimulationType::MonteCarlo:
      {
        MonteCarlo mc(inputReader);
        mc.run();
        break;
      }
      case InputReader::SimulationType::Breakthrough:
      {
        BreakthroughSimulation breakthrough(inputReader);
        breakthrough.run();
        break;
      }
      case InputReader::SimulationType::IAST:
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
