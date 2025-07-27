module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <fstream>
#include <numeric>
#include <span>
#include <string>
#include <utility>
#include <vector>
#endif

export module sample_movies;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import simulationbox;
import forcefield;

export struct SampleMovie
{
  SampleMovie(std::size_t systemId, std::size_t sampleEvery);

  void update(const ForceField &forceField, std::size_t systemId, const SimulationBox simulationBox,
              const std::span<Atom> atomPositions, std::size_t currentCycle);

  std::size_t sampleEvery{10};

  int modelNumber{1};
};
