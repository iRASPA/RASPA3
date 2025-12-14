module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#include "mdspanwrapper.h"
#endif

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
#include "mdspanwrapper.h"
#endif

module property_density_grid;

#ifdef USE_STD_IMPORT
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
//import mdspan;
#endif

// Gaussian cube file are stored row-order (std::layout_right)
// The grid is arranged with the x axis as the outer loop and the z axis as the inner loop

namespace
{
struct AxisStencil
{
  std::size_t left;
  std::size_t right;
  double wLeft;
  double wRight;
};

inline AxisStencil equitableAxis(double f, std::size_t n)
{
  f -= std::floor(f);

  const double g = static_cast<double>(n) * f;
  const long idx = static_cast<long>(std::ceil(g));
  const double part = g - std::floor(g);

  AxisStencil out{};

  if (part > 0.5)
  {
    if (idx >= static_cast<long>(n))
    {
      out.left = n - 1;
      out.right = 0;
    }
    else
    {
      out.left = static_cast<std::size_t>(idx - 1);
      out.right = static_cast<std::size_t>(idx);
    }
    out.wLeft = 1.5 - part;
    out.wRight = part - 0.5;
  }
  else
  {
    if (idx <= 1)
    {
      out.left = n - 1;
      out.right = 0;
    }
    else
    {
      out.left = static_cast<std::size_t>(idx - 2);
      out.right = static_cast<std::size_t>(idx - 1);
    }
    out.wLeft = 0.5 - part;
    out.wRight = 0.5 + part;
  }

  out.wLeft = std::max(0.0, out.wLeft);
  out.wRight = std::max(0.0, out.wRight);

  return out;
}

inline void depositEquitable(std::mdspan<double, std::dextents<std::size_t, 5>> grid,
                             std::size_t comp, std::size_t channel, double3 f, const double3 gridSize)
{
  const std::size_t Nx = static_cast<std::size_t>(gridSize.x);
  const std::size_t Ny = static_cast<std::size_t>(gridSize.y);
  const std::size_t Nz = static_cast<std::size_t>(gridSize.z);

  const AxisStencil A = equitableAxis(f.x, Nx);
  const AxisStencil B = equitableAxis(f.y, Ny);
  const AxisStencil C = equitableAxis(f.z, Nz);

  grid[comp, channel, A.left,  B.left,  C.left]  += A.wLeft  * B.wLeft  * C.wLeft;
  grid[comp, channel, A.right, B.left,  C.left]  += A.wRight * B.wLeft  * C.wLeft;
  grid[comp, channel, A.left,  B.right, C.left]  += A.wLeft  * B.wRight * C.wLeft;
  grid[comp, channel, A.left,  B.left,  C.right] += A.wLeft  * B.wLeft  * C.wRight;
  grid[comp, channel, A.right, B.right, C.left]  += A.wRight * B.wRight * C.wLeft;
  grid[comp, channel, A.right, B.left,  C.right] += A.wRight * B.wLeft  * C.wRight;
  grid[comp, channel, A.left,  B.right, C.right] += A.wLeft  * B.wRight * C.wRight;
  grid[comp, channel, A.right, B.right, C.right] += A.wRight * B.wRight * C.wRight;
}
}

void PropertyDensityGrid::sample(const std::optional<Framework> &framework, const SimulationBox &simulationBox,
                                 std::span<const Atom> moleculeAtoms, std::size_t currentCycle)
{
  if (currentCycle % sampleEvery != 0uz) return;

  const std::size_t nChannels = numberOfChannels;

  std::mdspan<double, std::dextents<std::size_t, 5>> data_cell(
      grid_cell.data(), numberOfComponents, nChannels, numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);
  std::mdspan<double, std::dextents<std::size_t, 5>> data_unitcell(
      grid_unitcell.data(), numberOfComponents, nChannels, numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);

  for (std::span<const Atom>::iterator it = moleculeAtoms.begin(); it != moleculeAtoms.end(); ++it)
  {
    if (it->isFractional) continue;

    const std::size_t comp = static_cast<std::size_t>(it->componentId);

    // Channel selection:
    // - If DensityGridPseudoAtomsList is empty: accumulate ALL atoms for the component into channel 0.
    // - If it is non-empty: accumulate only atoms whose pseudo atom type is in the list, into the
    //   corresponding channel.
    std::size_t channel = 0uz;
    if (!densityGridPseudoAtomsList.empty())
    {
      auto found = std::find(densityGridPseudoAtomsList.begin(), densityGridPseudoAtomsList.end(),
                             static_cast<std::size_t>(it->type));
      if (found == densityGridPseudoAtomsList.end()) continue;
      channel = static_cast<std::size_t>(std::distance(densityGridPseudoAtomsList.begin(), found));
    }

    const double3 pos = it->position;

    const double3 s = (simulationBox.inverseCell * pos).fract();
    if (binningMode == Binning::Standard)
    {
      data_cell[comp, channel, static_cast<std::size_t>(s.x * gridSize.x), static_cast<std::size_t>(s.y * gridSize.y),
                static_cast<std::size_t>(s.z * gridSize.z)]++;
    }
    else
    {
      depositEquitable(data_cell, comp, channel, s, gridSize);
    }

    if (framework.has_value())
    {
      const double3 t = (framework->simulationBox.inverseCell * pos).fract();
      if (binningMode == Binning::Standard)
      {
        data_unitcell[comp, channel, static_cast<std::size_t>(t.x * gridSize.x), static_cast<std::size_t>(t.y * gridSize.y),
                      static_cast<std::size_t>(t.z * gridSize.z)]++;
      }
      else
      {
        depositEquitable(data_unitcell, comp, channel, t, gridSize);
      }
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

  const std::size_t nChannels = numberOfChannels;

  std::mdspan<double, std::dextents<std::size_t, 5>> data_cell(
      grid_cell.data(), numberOfComponents, nChannels, numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);
  std::mdspan<double, std::dextents<std::size_t, 5>> data_unitcell(
      grid_unitcell.data(), numberOfComponents, nChannels, numberOfGridPoints.x, numberOfGridPoints.y, numberOfGridPoints.z);

  for (std::size_t i = 0; i < components.size(); ++i)
  {
    for (std::size_t c = 0; c < nChannels; ++c)
    {
      std::string channelName;
      if (!densityGridPseudoAtomsList.empty())
      {
        const std::size_t pseudoType = densityGridPseudoAtomsList[c];
        channelName = forceField.pseudoAtoms[pseudoType].name;
      }

      std::vector<double>::iterator it_begin =
          grid_cell.begin() + static_cast<std::vector<double>::difference_type>(((i * nChannels) + c) * totalGridSize);
      std::vector<double>::iterator it_end = it_begin + static_cast<std::vector<double>::difference_type>(totalGridSize);
      std::vector<double>::iterator maximum = std::max_element(it_begin, it_end);
      if (maximum == it_end || *maximum <= 0.0)
      {
        continue;
      }
      std::ofstream ostream(
          channelName.empty()
              ? std::format("density_grids/grid_cell_component_{}.s{}.cube", components[i].name, systemId)
              : std::format("density_grids/grid_cell_component_{}_{}.s{}.cube", components[i].name, channelName, systemId));
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

    //[[maybe_unused]]std::vector<double>::iterator maximum = std::max_element(it_begin, it_end);
      double normalization{1.0};
    switch (normType)
    {
      case Normalization::Max:
      {
        normalization = (*maximum > 0.0) ? (1.0 / *maximum) : 1.0;
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
  }

  if (framework.has_value())
  {
    for (std::size_t i = 0; i < components.size(); ++i)
    {
      for (std::size_t c = 0; c < nChannels; ++c)
      {
        std::string channelName;
        if (!densityGridPseudoAtomsList.empty())
        {
          const std::size_t pseudoType = densityGridPseudoAtomsList[c];
          channelName = forceField.pseudoAtoms[pseudoType].name;
        }

        std::vector<double>::iterator it_begin =
            grid_unitcell.begin() + static_cast<std::vector<double>::difference_type>(((i * nChannels) + c) * totalGridSize);
        std::vector<double>::iterator it_end =
            it_begin + static_cast<std::vector<double>::difference_type>(totalGridSize);
        std::vector<double>::iterator maximum = std::max_element(it_begin, it_end);
        if (maximum == it_end || *maximum <= 0.0)
        {
          continue;
        }
        std::ofstream ostream(
            channelName.empty()
                ? std::format("density_grids/grid_unitcell_component_{}.s{}.cube", components[i].name, systemId)
                : std::format("density_grids/grid_unitcell_component_{}_{}.s{}.cube", components[i].name, channelName, systemId));

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


        double normalization{1.0};
      switch (normType)
      {
        case Normalization::Max:
        {
        normalization = (*maximum > 0.0) ? (1.0 / *maximum) : 1.0;
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

  archive << temp.densityGridPseudoAtomsList;

  archive << temp.normType;
  archive << temp.binningMode;
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

  if (versionNumber >= 3)
  {
    archive >> temp.densityGridPseudoAtomsList;
  }
  else
  {
    temp.densityGridPseudoAtomsList.clear();
  }

  temp.numberOfChannels =
      temp.densityGridPseudoAtomsList.empty() ? 1uz : temp.densityGridPseudoAtomsList.size();

  archive >> temp.normType;
  if (versionNumber >= 4)
  {
    archive >> temp.binningMode;
  }
  else
  {
    temp.binningMode = PropertyDensityGrid::Binning::Standard;
  }
  archive >> temp.numberOfSamples;

  #if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(
        std::format("PropertyDensityGrid: Error in binary restart\n"));
  }
  #endif

  return archive;
}