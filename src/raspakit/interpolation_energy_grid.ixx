module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <istream>
#include <ostream>
#include <print>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <vector>
#include <utility>
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
import simulationbox;
import scaling;
import forcefield;
import framework;
import interactions_framework_molecule_grid;

export struct InterpolationEnergyGrid
{
  uint64_t versionNumber{1};

  SimulationBox unitCellBox;
  int3 numberOfCells;
  int3 numberOfGridPoints;
  ForceField::InterpolationScheme order;
  std::vector<double> data;

  InterpolationEnergyGrid() {}

  InterpolationEnergyGrid(const SimulationBox unitCellBox, int3 numberOfCells, ForceField::InterpolationScheme order)
      : unitCellBox(unitCellBox),
        numberOfCells(numberOfCells),
        numberOfGridPoints(numberOfCells.x + 1, numberOfCells.y + 1, numberOfCells.z + 1),
        order(order),
        data(std::to_underlying(order) *
             static_cast<size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z))
  {
  }

  void makeInterpolationGrid(std::ostream &stream, ForceField::InterpolationGridType interpolationGridType,
                             const ForceField &forceField, const Framework &framework, 
                             double cutOff, size_t pseudo_atom_index);

  double interpolate(double3 pos) const;
  std::pair<double, double3> interpolateGradient(double3 pos) const;
  std::tuple<double, double3, double3x3> interpolateHessian(double3 pos) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const InterpolationEnergyGrid &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, InterpolationEnergyGrid &s);

};
