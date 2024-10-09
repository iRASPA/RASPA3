module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#endif

export module breakthrough_simulation;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <span>;
import <tuple>;
import <string>;
import <fstream>;
#endif

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
