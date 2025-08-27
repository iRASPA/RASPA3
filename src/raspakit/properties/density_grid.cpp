module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <complex>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <print>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

module property_density_grid;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import int3;
import double3;
import double3x3;
import atom;
import simulationbox;
import forcefield;
import averages;
import stringutils;
import framework;
import component;
#if !(defined(__has_include) && __has_include(<mdspan>))
import mdspan;
#endif

// Gaussian cube file are stored row-order (std::layout_right)
// The grid is arranged with the x axis as the outer loop and the z axis as the inner loop

void PropertyDensityGrid::sample(const std::optional<Framework> &framework, const SimulationBox &simulationBox,
                                 std::span<const Atom> moleculeAtoms, std::size_t currentCycle)
{
  if (currentCycle % sampleEvery != 0uz) return;

  std::mdspan<double, std::dextents<std::size_t, 4>> data_cell(
      grid_cell.data(), numberOfComponents, numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);
  std::mdspan<double, std::dextents<std::size_t, 4>> data_unitcell(
      grid_unitcell.data(), numberOfComponents, numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);

  for (std::span<const Atom>::iterator it = moleculeAtoms.begin(); it != moleculeAtoms.end(); ++it)
  {
    if (it->isFractional) continue;
    std::size_t comp = static_cast<std::size_t>(it->componentId);
    double3 pos = it->position;
    double3 s = (simulationBox.inverseCell * pos).fract();
    data_cell[comp, static_cast<std::size_t>(s.x * gridSize.x), static_cast<std::size_t>(s.y * gridSize.y),
              static_cast<std::size_t>(s.z * gridSize.z)]++;

    if (framework.has_value())
    {
      double3 t = (framework->simulationBox.inverseCell * pos).fract();
      data_unitcell[comp, static_cast<std::size_t>(t.x * gridSize.x), static_cast<std::size_t>(t.y * gridSize.y),
                    static_cast<std::size_t>(t.z * gridSize.z)]++;
    }
  }

  numberOfSamples++;
}

void PropertyDensityGrid::writeOutput(std::size_t systemId, [[maybe_unused]] const SimulationBox &simulationBox,
                                      const ForceField &forceField, const std::optional<Framework> &framework,
                                      const std::vector<Component> &components, std::size_t currentCycle)
{
  if (currentCycle % writeEvery != 0uz) return;

  std::filesystem::create_directory("density_grids");

  std::mdspan<double, std::dextents<std::size_t, 4>> data_cell(
      grid_cell.data(), numberOfComponents, numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);
  std::mdspan<double, std::dextents<std::size_t, 4>> data_unitcell(
      grid_unitcell.data(), numberOfComponents, numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);

  for (std::size_t i = 0; i < components.size(); ++i)
  {
    std::ofstream ostream(std::format("density_grids/grid_component_{}.s{}.cube", components[i].name, systemId));
    const double3x3 cell = simulationBox.cell;

    std::vector<Atom> frameworkAtoms = framework.has_value() ? framework->atoms : std::vector<Atom>{};

    std::print(ostream, "Cube density file\n");
    std::print(ostream, "Written by RASPA-3\n");
    std::print(ostream, "{} {} {} {}\n", frameworkAtoms.size(), 0.0, 0.0, 0.0);

    std::print(ostream, "{} {} {} {}\n", -numberOfGridPoints.x, cell.ax / gridSize.x, cell.ay / gridSize.x,
               cell.az / gridSize.x);
    std::print(ostream, "{} {} {} {}\n", -numberOfGridPoints.y, cell.bx / gridSize.y, cell.by / gridSize.y,
               cell.bz / gridSize.y);
    std::print(ostream, "{} {} {} {}\n", -numberOfGridPoints.z, cell.cx / gridSize.z, cell.cy / gridSize.z,
               cell.cz / gridSize.z);
    for (std::vector<Atom>::iterator it = frameworkAtoms.begin(); it != frameworkAtoms.end(); ++it)
    {
      double3 pos = it->position;
      std::size_t type = static_cast<std::size_t>(it->type);
      double charge = it->charge;
      std::size_t atomicNumber = forceField.pseudoAtoms[type].atomicNumber;

      std::print(ostream, "{} {} {} {} {}\n", atomicNumber, charge, pos.x, pos.y, pos.z);
    }

    std::vector<double>::iterator it_begin =
        grid_cell.begin() + std::distance(&data_cell[0, 0, 0, 0], &data_cell[i, 0, 0, 0]);
    std::vector<double>::iterator it_end = it_begin + static_cast<std::vector<double>::difference_type>(totalGridSize);

    //[[maybe_unused]]std::vector<double>::iterator maximum = std::max_element(it_begin, it_end);
    double normalization{1.0};
    switch (normType)
    {
      case Normalization::Max:
      {
        std::vector<double>::iterator maximum = std::max_element(it_begin, it_end);
        if (maximum != it_end)
        {
          normalization = static_cast<double>(1.0 / *maximum);
        }
        break;
      }
      case Normalization::NumberDensity:
      {
        normalization =
            (gridSize.x * gridSize.y * gridSize.z) / (cell.determinant() * static_cast<double>(numberOfSamples));
        break;
      }
    }

    for (std::vector<double>::iterator it = it_begin; it != it_end; ++it)
    {
      std::print(ostream, "{}\n", *it * normalization);
    }
  }

  if (framework.has_value())
  {
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      std::ofstream ostream(
          std::format("density_grids/grid_unitcell_component_{}.s{}.cube", components[i].name, systemId));

      std::vector<Atom> frameworkAtoms = framework->unitCellAtoms;
      double3x3 unitCell = framework->simulationBox.cell;

      std::print(ostream, "Cube density file\n");
      std::print(ostream, "Written by RASPA-3\n");
      std::print(ostream, "{} {} {} {}\n", frameworkAtoms.size(), 0.0, 0.0, 0.0);

      std::print(ostream, "{} {} {} {}\n", -numberOfGridPoints.x, unitCell.ax / gridSize.x, unitCell.ay / gridSize.x,
                 unitCell.az / gridSize.x);
      std::print(ostream, "{} {} {} {}\n", -numberOfGridPoints.y, unitCell.bx / gridSize.y, unitCell.by / gridSize.y,
                 unitCell.bz / gridSize.y);
      std::print(ostream, "{} {} {} {}\n", -numberOfGridPoints.z, unitCell.cx / gridSize.z, unitCell.cy / gridSize.z,
                 unitCell.cz / gridSize.z);
      for (std::vector<Atom>::iterator it = frameworkAtoms.begin(); it != frameworkAtoms.end(); ++it)
      {
        double3 pos = it->position;
        std::size_t type = static_cast<std::size_t>(it->type);
        double charge = it->charge;
        std::size_t atomicNumber = forceField.pseudoAtoms[type].atomicNumber;

        std::print(ostream, "{} {} {} {} {}\n", atomicNumber, charge, pos.x, pos.y, pos.z);
      }

      std::vector<double>::iterator it_begin =
          grid_unitcell.begin() + std::distance(&data_unitcell[0, 0, 0, 0], &data_unitcell[i, 0, 0, 0]);
      std::vector<double>::iterator it_end =
          it_begin + static_cast<std::vector<double>::difference_type>(totalGridSize);

      double normalization{1.0};
      switch (normType)
      {
        case Normalization::Max:
        {
          std::vector<double>::iterator maximum = std::max_element(it_begin, it_end);
          if (maximum != it_end)
          {
            normalization = static_cast<double>(1.0 / *maximum);
          }
          break;
        }
        case Normalization::NumberDensity:
        {
          double3 numberOfUnitCells = simulationBox.lengths() / framework->simulationBox.lengths();
          normalization = (gridSize.x * gridSize.y * gridSize.z) /
                          (unitCell.determinant() * static_cast<double>(numberOfSamples) * numberOfUnitCells.x *
                           numberOfUnitCells.y * numberOfUnitCells.z);
          break;
        }
      }

      for (std::vector<double>::iterator it = it_begin; it != it_end; ++it)
      {
        std::print(ostream, "{}\n", *it * normalization);
      }
    }
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyDensityGrid &temp)
{
  archive << temp.versionNumber;

  archive << temp.numberOfFrameworks;
  archive << temp.numberOfComponents;
  archive << temp.grid_cell;
  archive << temp.grid_unitcell;
  archive << temp.totalGridSize;
  archive << temp.numberOfGridPoints;
  archive << temp.gridSize;

  archive << temp.sampleEvery;
  archive << temp.writeEvery;

  archive << temp.normType;
  archive << temp.numberOfSamples;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyDensityGrid &temp)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > temp.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyDensityGrid' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> temp.numberOfFrameworks;
  archive >> temp.numberOfComponents;
  archive >> temp.grid_cell;
  archive >> temp.grid_unitcell;
  archive >> temp.totalGridSize;
  archive >> temp.numberOfGridPoints;
  archive >> temp.gridSize;

  archive >> temp.sampleEvery;
  archive >> temp.writeEvery;

  archive >> temp.normType;
  archive >> temp.numberOfSamples;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyDensityGrid: Error in binary restart\n"));
  }
#endif

  return archive;
}
