module;

#ifdef USE_LEGACY_HEADERS
#include <memory>
#include <vector>
#endif

export module breakthrough_simulation;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <memory>;
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
  std::vector<std::shared_ptr<System>> systems;
};
