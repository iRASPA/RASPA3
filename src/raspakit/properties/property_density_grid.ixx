module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <mdspan>
#include <optional>
#include <span>
#include <string>
#include <vector>
#endif

export module property_density_grid;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <array>;
import <optional>;
import <cmath>;
import <string>;
import <span>;
import <algorithm>;
import <mdspan>;
#endif

import archive;
import int3;
import double3;

import atom;
import simulationbox;
import forcefield;
import framework;
import component;

export struct PropertyDensityGrid
{
  PropertyDensityGrid() {}

  PropertyDensityGrid(size_t numberOfFrameworks, size_t numberOfComponents, int3 numberOfGridPoints, size_t sampleEvery,
                      size_t writeEvery, std::vector<size_t> densityGridPseudoAtomsList)
      : numberOfFrameworks(numberOfFrameworks),
        numberOfComponents(numberOfComponents),
        grid_cell(numberOfComponents *
                  static_cast<size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z)),
        grid_unitcell(std::min(1uz, numberOfFrameworks) * numberOfComponents *
                      static_cast<size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z)),
        totalGridSize(static_cast<size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z)),
        numberOfGridPoints(numberOfGridPoints),
        gridSize(static_cast<double>(numberOfGridPoints.x), static_cast<double>(numberOfGridPoints.y),
                 static_cast<double>(numberOfGridPoints.z)),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        densityGridPseudoAtomsList(densityGridPseudoAtomsList)
  {
  }

  uint64_t versionNumber{1};

  size_t numberOfFrameworks;
  size_t numberOfComponents;
  std::vector<double> grid_cell;
  std::vector<double> grid_unitcell;
  size_t totalGridSize;
  int3 numberOfGridPoints;
  double3 gridSize;
  size_t sampleEvery;
  size_t writeEvery;
  std::vector<size_t> densityGridPseudoAtomsList;

  void sample(const std::vector<Framework> &frameworks, const SimulationBox &simulationBox,
              std::span<const Atom> moleculeAtoms, size_t currrentCycle);
  void writeOutput(size_t systemId, const SimulationBox &simulationBox, const ForceField &forceField,
                   const std::vector<Framework> &frameworkComponents, const std::vector<Component> &components,
                   size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyDensityGrid &temp);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyDensityGrid &temp);
};
