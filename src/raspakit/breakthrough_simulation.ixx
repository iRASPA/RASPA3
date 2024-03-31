export module breakthrough_simulation;

import <vector>;
import <span>;
import <tuple>;
import <string>;
import <fstream>;

import input_reader;
import component;
import system;
import mixture_prediction;

export struct BreakthroughSimulation
{
  public:
    BreakthroughSimulation(InputReader &inputreader);

    void run();

  private:
    std::vector<System> systems;
};
