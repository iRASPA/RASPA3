export module property_density_grid;

import <vector>;
import <array>;
import <optional>;
import <cmath>;
import <string>;
import <span>;
import <algorithm>;
#if defined(__has_include) && __has_include(<mdspan>)
  import <mdspan>;
#else
  import mdspan;
#endif

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

  PropertyDensityGrid(size_t numberOfFrameworks, size_t numberOfComponents, int3 numberOfGridPoints, size_t sampleEvery, size_t writeEvery) :
    grid_cell(numberOfComponents * static_cast<size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z)),
    data_cell(grid_cell.data(), numberOfComponents, numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z),
    grid_unitcell(std::min(1uz, numberOfFrameworks) * numberOfComponents * static_cast<size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z)),
    data_unitcell(grid_unitcell.data(), numberOfComponents, std::min(1uz, numberOfFrameworks), numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z),
    totalGridSize(static_cast<size_t>(numberOfGridPoints.x * numberOfGridPoints.y * numberOfGridPoints.z)),
    numberOfGridPoints(numberOfGridPoints),
    gridSize(static_cast<double>(numberOfGridPoints.x), static_cast<double>(numberOfGridPoints.y), static_cast<double>(numberOfGridPoints.z)),
    sampleEvery(sampleEvery),
    writeEvery(writeEvery)
  {
  }

  std::vector<double> grid_cell;
  #if defined(__has_include) && __has_include(<mdspan>)
    std::mdspan<double, std::dextents<size_t, 4>> data_cell;
  #else
    std::experimental::mdspan<double, std::experimental::dextents<size_t, 4>> data_cell;
  #endif
  std::vector<double> grid_unitcell;
  #if defined(__has_include) && __has_include(<mdspan>)
    std::mdspan<double, std::dextents<size_t, 5>> data_unitcell;
  #else
    std::experimental::mdspan<double, std::experimental::dextents<size_t, 5>> data_unitcell;
  #endif
  size_t totalGridSize;
  int3 numberOfGridPoints;
  double3 gridSize;

  size_t sampleEvery;
  size_t writeEvery;

  void sample(const std::vector<Framework> &frameworks, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms, size_t currrentCycle);
  void writeOutput(size_t systemId, const SimulationBox &simulationBox,
                   const ForceField &forceField, 
                   const std::vector<Framework> &frameworkComponents,
                   const std::vector<Component> &components,
                   size_t currentCycle);
};


