module;

#ifdef USE_LEGACY_HEADERS
#include <array>
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
import <array>;
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
import double3x3;
import stringutils;

import scaling;
import forcefield;
import framework;
import interactions_framework_molecule;

export struct InterpolationEnergyGrid
{
  int3 numberOfCells;
  int3 numberOfGridPoints;
  std::vector<double> data;

  InterpolationEnergyGrid(int3 numberOfCells): 
    numberOfCells(numberOfCells),
    numberOfGridPoints(numberOfCells.x + 1, numberOfCells.y + 1, numberOfCells.z + 1),
    //data(static_cast<size_t>(8 * numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z))
    data(static_cast<size_t>(27 * numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z))
  {
  }

  void makeVDWGrid(const ForceField &forceField, const Framework &framework, size_t pseudo_atom_index);
  double interpolateVDWGrid(double3 s);
};
