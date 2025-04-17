module;

#ifdef USE_LEGACY_HEADERS
#include <memory>
#include <vector>
#endif

export module isotherm_fitting_simulation;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <memory>;
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
  std::vector<std::shared_ptr<System>> systems;
};
