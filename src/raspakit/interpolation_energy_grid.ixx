module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

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
#include <utility>
#include <vector>
#endif

export module interpolation_energy_grid;

#ifdef USE_STD_IMPORT
import std;
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
  std::uint64_t versionNumber{1};

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
             static_cast<std::size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z))
  {
  }

  constexpr static std::make_signed_t<std::size_t> num_points_interpolation{6};

  void makeInterpolationGrid(std::ostream &stream, ForceField::InterpolationGridType interpolationGridType,
                             const ForceField &forceField, const Framework &framework, double cutOff,
                             std::size_t pseudo_atom_index);

  double interpolate(double3 pos) const;
  std::pair<double, double3> interpolateGradient(double3 pos) const;
  std::tuple<double, double3, double3x3> interpolateHessian(double3 pos) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const InterpolationEnergyGrid &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, InterpolationEnergyGrid &s);
};
