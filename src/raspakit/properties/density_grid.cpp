module;

module property_density_grid;

import std;

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
  // (the fractional part of a coordinate just below zero rounds to exactly one: that is the origin)
  if (!(f < 1.0)) f = 0.0;

  // The voxel that holds the point, counted from one: floor(g) + 1. It was ceil(g), which is the same number except
  // for a point on a voxel boundary (g an integer, as at the special positions 1/4, 1/2, 3/4 of a grid whose size is
  // a multiple of four): there it is one less, and the two halves of the point went to the voxels one to the left of
  // the two that share the boundary; at g = 1 one half went to the last voxel, across the whole box.
  const double g = std::min(static_cast<double>(n) * f, std::nextafter(static_cast<double>(n), 0.0));
  const long idx = static_cast<long>(std::floor(g)) + 1L;
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
                             std::size_t comp, std::size_t channel, double3 f, const double3 gridSize,
                             double weight)
{
  const std::size_t Nx = static_cast<std::size_t>(gridSize.x);
  const std::size_t Ny = static_cast<std::size_t>(gridSize.y);
  const std::size_t Nz = static_cast<std::size_t>(gridSize.z);

  const AxisStencil A = equitableAxis(f.x, Nx);
  const AxisStencil B = equitableAxis(f.y, Ny);
  const AxisStencil C = equitableAxis(f.z, Nz);

  grid[comp, channel, A.left,  B.left,  C.left]  += weight * A.wLeft  * B.wLeft  * C.wLeft;
  grid[comp, channel, A.right, B.left,  C.left]  += weight * A.wRight * B.wLeft  * C.wLeft;
  grid[comp, channel, A.left,  B.right, C.left]  += weight * A.wLeft  * B.wRight * C.wLeft;
  grid[comp, channel, A.left,  B.left,  C.right] += weight * A.wLeft  * B.wLeft  * C.wRight;
  grid[comp, channel, A.right, B.right, C.left]  += weight * A.wRight * B.wRight * C.wLeft;
  grid[comp, channel, A.right, B.left,  C.right] += weight * A.wRight * B.wLeft  * C.wRight;
  grid[comp, channel, A.left,  B.right, C.right] += weight * A.wLeft  * B.wRight * C.wRight;
  grid[comp, channel, A.right, B.right, C.right] += weight * A.wRight * B.wRight * C.wRight;
}
}

namespace
{
// The grid point of a fractional coordinate. The fractional part s - floor(s) of a coordinate just below zero
// rounds to exactly one, and one times the number of points is the index past the last point: the count was
// written outside the grid. That coordinate is the last point.
inline std::size_t gridIndex(double fraction, double size, std::size_t numberOfPoints)
{
  if (numberOfPoints == 0uz) return 0uz;
  const double scaled = fraction * size;
  if (!(scaled > 0.0)) return 0uz;
  return std::min(static_cast<std::size_t>(scaled), numberOfPoints - 1uz);
}
}  // namespace

void PropertyDensityGrid::sample(const std::optional<Framework> &framework, const SimulationBox &simulationBox,
                                 std::span<const Atom> moleculeAtoms, std::size_t currentCycle, double weight)
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
      data_cell[comp, channel, gridIndex(s.x, gridSize.x, numberOfGridPoints.x),
                gridIndex(s.y, gridSize.y, numberOfGridPoints.y), gridIndex(s.z, gridSize.z, numberOfGridPoints.z)] += weight;
    }
    else
    {
      depositEquitable(data_cell, comp, channel, s, gridSize, weight);
    }

    if (framework.has_value())
    {
      // the unit cell of the current (possibly deformed) box: the supercell fraction times the number of cells
      const double3 t = double3(s.x * static_cast<double>(framework->numberOfUnitCells.x),
                                s.y * static_cast<double>(framework->numberOfUnitCells.y),
                                s.z * static_cast<double>(framework->numberOfUnitCells.z))
                            .fract();
      if (binningMode == Binning::Standard)
      {
        data_unitcell[comp, channel, gridIndex(t.x, gridSize.x, numberOfGridPoints.x),
                      gridIndex(t.y, gridSize.y, numberOfGridPoints.y),
                      gridIndex(t.z, gridSize.z, numberOfGridPoints.z)] += weight;
      }
      else
      {
        depositEquitable(data_unitcell, comp, channel, t, gridSize, weight);
      }
    }
  }

  // the grids, the volume and the number of samples carry the weight of the configuration (a biased walk visits
  // the loadings with other frequencies than the ensemble: the integrated density was that of the biased walk)
  numberOfSamples++;
  weightSum += weight;
  volumeSum += weight * simulationBox.volume;
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
        // the average volume of the samples: the current volume differs from it in a variable-cell run
        const double averageVolume =
            volumeSum > 0.0 ? volumeSum / weightSum : cell.determinant();
        normalization = (gridSize.x * gridSize.y * gridSize.z) / (averageVolume * weightSum);
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
          // every unit cell of the box folds onto the grid: the average unit-cell volume times the number of cells
          // is the average volume of the box
          const double averageVolume =
              volumeSum > 0.0 ? volumeSum / weightSum : simulationBox.volume;
          normalization =
              (gridSize.x * gridSize.y * gridSize.z) / (averageVolume * weightSum);
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
  archive << temp.volumeSum;
  archive << temp.weightSum;

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

  archive >> temp.densityGridPseudoAtomsList;

  temp.numberOfChannels =
      temp.densityGridPseudoAtomsList.empty() ? 1uz : temp.densityGridPseudoAtomsList.size();

  archive >> temp.normType;
  archive >> temp.binningMode;
  archive >> temp.numberOfSamples;
  if (versionNumber >= 2) archive >> temp.volumeSum;
  // (before version 3 every sample had the weight one)
  temp.weightSum = static_cast<double>(temp.numberOfSamples);
  if (versionNumber >= 3) archive >> temp.weightSum;

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
