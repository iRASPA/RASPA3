module;

export module mixture_prediction_simulation;

import std;

import input_reader;
import component;
import system;
import multi_site_isotherm;
import mixture_prediction;

export struct MixturePredictionSimulation
{
 public:
  MixturePredictionSimulation(InputReader &inputreader);

  void run();

 private:
  std::vector<System> systems;
};
