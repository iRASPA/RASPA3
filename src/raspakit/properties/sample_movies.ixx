module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <numeric>
#include <span>
#include <string>
#include <utility>
#include <vector>
#endif

export module sample_movies;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <span>;
import <numeric>;
import <fstream>;
import <utility>;
import <string>;
#endif

import atom;
import simulationbox;
import forcefield;

export struct SampleMovie
{
  SampleMovie(size_t systemId, size_t sampleEvery);

  void update(const ForceField &forceField, size_t systemId, const SimulationBox simulationBox,
              const std::span<Atom> atomPositions, size_t currentCycle);

  size_t sampleEvery{10};

  int modelNumber{1};
};
