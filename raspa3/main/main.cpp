
import <exception>;
import <iostream>;
import <vector>;
import <deque>;
import <optional>;
import <semaphore>;

import monte_carlo;
import input_reader;
import threadpool;
import threading;

int main()
{
  try
  {
    InputReader inputReader("pseudo_atoms.def", "force_field_mixing_rules.def",
                          "force_field.def", "simulation.input");

    ThreadPool::createPool(inputReader.numberOfThreads, inputReader.threadingType);

    switch (inputReader.simulationType)
    {
      case InputReader::SimulationType::MonteCarlo:
      {
        MonteCarlo mc(inputReader);
        mc.run();
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
