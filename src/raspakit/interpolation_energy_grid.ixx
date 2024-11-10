module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <istream>
#include <ostream>
#include <print>
#include <sstream>
#include <type_traits>
#endif

export module interpolation_energy_grid;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <cmath>;
import <cstddef>;
import <istream>;
import <ostream>;
import <fstream>;
import <sstream>;
import <type_traits>;
import <print>;
#endif

import archive;
import int3;
import double3;
import stringutils;

import scaling;
import forcefield;
import framework;

export struct InterpolationEnergyGrid
{
  int3 numberOfGridPoints;
  std::vector<double> data;

  InterpolationEnergyGrid(int3 numberOfGridPoints): 
    numberOfGridPoints(numberOfGridPoints),
    data(static_cast<size_t>(8 * numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z))
  {
  }

  //std::vector<double> createVDWGrid(const ForceField &forcefield, const std::vector<Framework> & frameworkComponents, size_t typeA);
};
