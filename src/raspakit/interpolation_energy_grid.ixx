module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#include "mdspanwrapper.h"
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
#include <string>
#include "mdspanwrapper.h"
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
  double3 origin;
  int3 numberOfGridPoints;
  int3 numberOfCells;
  ForceField::InterpolationScheme order;
  std::vector<double> data;

  InterpolationEnergyGrid() {}

  InterpolationEnergyGrid(const SimulationBox unitCellBox, double3 origin, int3 numberOfGridPoints, ForceField::InterpolationScheme order);

  constexpr static std::make_signed_t<std::size_t> num_points_interpolation{6};

  void makeExternalFieldInterpolationGrid(std::ostream& stream, const ForceField& forceField, const SimulationBox &simulationBox);

  void makeFrameworkInterpolationGrid(std::ostream &stream, ForceField::InterpolationGridType interpolationGridType,
                             const ForceField &forceField, const Framework &framework, double cutOff,
                             std::size_t pseudo_atom_index);

  double interpolate(double3 pos) const;
  std::pair<double, double3> interpolateGradient(double3 pos) const;
  std::tuple<double, double3, double3x3> interpolateHessian(double3 pos) const;

  void writeOutput(std::size_t systemId, const SimulationBox &simulationBox, const ForceField &forceField);

  static const int3 parseExternalFieldGridDimensions(const std::string& filename);

  /**
   * \brief Parses external field grid based on the provided CUBE file.
   *
   * This function extracts parameters specific to external field grids,
   * allowing to use this data in construction of external field interpolation grid in simulations.
   *
   * \param filename The name of the CUBE file containing the external field data.
   * \return A pair containing the 3D grid of external field values and the grid dimensions.
   */
  static const std::vector<double> parseExternalFieldGridCube(const std::string& filename);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const InterpolationEnergyGrid &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, InterpolationEnergyGrid &s);
};
