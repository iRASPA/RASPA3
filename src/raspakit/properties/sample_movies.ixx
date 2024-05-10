module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <span>
#include <numeric>
#include <fstream>
#include <utility>
#include <string>
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
    SampleMovie(size_t sampleEvery):
      sampleEvery(sampleEvery)
    {}

    void update(const ForceField &forceField, size_t systemId, const SimulationBox simulationBox, const std::span<Atom> atomPositions, size_t currentCycle);

    size_t sampleEvery{ 10 };

    int modelNumber{ 1 };
};
