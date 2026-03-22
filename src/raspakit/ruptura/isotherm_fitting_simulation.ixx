module;

export module isotherm_fitting_simulation;

import std;

import input_reader;
import component;
import system;
import multi_site_isotherm;
import isotherm_fitting;

export struct IsothermFittingSimulation
{
 public:
  IsothermFittingSimulation(InputReader &inputreader);

  void run();

 private:
  std::vector<System> systems;
};
