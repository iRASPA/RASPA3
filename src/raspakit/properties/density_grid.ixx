module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <vector>
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

export module property_density_grid;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import int3;
import double3;

import atom;
import simulationbox;
import forcefield;
import framework;
import component;
#if !(defined(__has_include) && __has_include(<mdspan>))
import mdspan;
#endif

export struct PropertyDensityGrid
{
  enum class Normalization : std::size_t
  {
    Max = 0,
    NumberDensity = 1
  };

  PropertyDensityGrid() {}

  PropertyDensityGrid(std::size_t numberOfFrameworks, std::size_t numberOfComponents, int3 numberOfGridPoints,
                      std::size_t sampleEvery, std::size_t writeEvery,
                      std::vector<std::size_t> densityGridPseudoAtomsList, Normalization normType)
      : numberOfFrameworks(numberOfFrameworks),
        numberOfComponents(numberOfComponents),
        grid_cell(numberOfComponents *
                  static_cast<std::size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z)),
        grid_unitcell(numberOfComponents *
                      static_cast<std::size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z)),
        totalGridSize(static_cast<std::size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z)),
        numberOfGridPoints(numberOfGridPoints),
        gridSize(static_cast<double>(numberOfGridPoints.x), static_cast<double>(numberOfGridPoints.y),
                 static_cast<double>(numberOfGridPoints.z)),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        densityGridPseudoAtomsList(densityGridPseudoAtomsList),
        normType(normType)
  {
  }

  std::uint64_t versionNumber{2};

  std::size_t numberOfFrameworks;
  std::size_t numberOfComponents;
  std::vector<double> grid_cell;
  std::vector<double> grid_unitcell;
  std::size_t totalGridSize;
  int3 numberOfGridPoints;
  double3 gridSize;
  std::size_t sampleEvery;
  std::size_t writeEvery;
  std::vector<std::size_t> densityGridPseudoAtomsList;
  Normalization normType{Normalization::Max};
  std::size_t numberOfSamples{0};

  void sample(const std::optional<Framework> &frameworks, const SimulationBox &simulationBox,
              std::span<const Atom> moleculeAtoms, std::size_t currrentCycle);
  void writeOutput(std::size_t systemId, const SimulationBox &simulationBox, const ForceField &forceField,
                   const std::optional<Framework> &frameworkComponents, const std::vector<Component> &components,
                   std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyDensityGrid &temp);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyDensityGrid &temp);
};
