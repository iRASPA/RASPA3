module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <span>
#include <tuple>
#include <string>
#include <fstream>
#endif

export module mixture_prediction_simulation;

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
