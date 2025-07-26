module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <fstream>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#endif

export module isotherm_fitting_simulation;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

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
