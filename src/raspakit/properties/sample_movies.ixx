module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <fstream>
#include <numeric>
#include <ostream>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#endif

export module sample_movies;

#ifdef USE_STD_IMPORT
import std;
#endif

import archive;
import atom;
import simulationbox;
import forcefield;

export struct SampleMovie
{
  std::uint64_t versionNumber{1};

  SampleMovie() {};
  SampleMovie(std::size_t systemId, std::size_t sampleEvery);

  void update(const ForceField& forceField, std::size_t systemId, const SimulationBox simulationBox,
              const std::span<Atom> atomData, std::size_t currentCycle);

  std::size_t sampleEvery{10};

  int modelNumber{1};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const SampleMovie& box);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, SampleMovie& box);
};
