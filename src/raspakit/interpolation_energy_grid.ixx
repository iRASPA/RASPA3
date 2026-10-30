module;

export module interpolation_energy_grid;

import std;

import archive;
import int3;
import uint3;
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
  uint3 numberOfGridPoints;
  uint3 numberOfCells;
  ForceFieldSettings::InterpolationScheme order;

  // omit data from binary restart (the grid will recreated on the fly)
  std::vector<double> data;

  InterpolationEnergyGrid() {}

  InterpolationEnergyGrid(const SimulationBox unitCellBox, double3 origin, uint3 numberOfGridPoints, ForceFieldSettings::InterpolationScheme order);

  constexpr static std::make_signed_t<std::size_t> num_points_interpolation{6};

  void makeExternalFieldInterpolationGrid(std::ostream& stream, const ForceField& forceField, const SimulationBox &simulationBox);

  void makeFrameworkInterpolationGrid(std::ostream &stream, ForceFieldSettings::InterpolationGridType interpolationGridType,
                             const ForceField &forceField, const Framework &framework, double cutOff,
                             std::size_t pseudo_atom_index);

  double interpolate(double3 pos) const;
  std::pair<double, double3> interpolateGradient(double3 pos) const;
  std::tuple<double, double3, double3x3> interpolateHessian(double3 pos) const;

  /// The energies of this grid and of 'other' (typically the Van der Waals and the Coulomb grid of the same
  /// framework) at 'pos'. The position is located once and the corner data of both grids are fetched together, so
  /// the two memory stalls overlap; equivalent to calling 'interpolate' on each grid.
  std::pair<double, double> interpolate(const InterpolationEnergyGrid& other, double3 pos) const;
  std::pair<std::pair<double, double3>, std::pair<double, double3>> interpolateGradient(
      const InterpolationEnergyGrid& other, double3 pos) const;

  /// Starts fetching the corner data of the cell that contains 'pos' into the cache. The grid is far larger
  /// than the caches and a probe position is effectively random, so an evaluation otherwise waits for memory;
  /// issuing the prefetches for a batch of positions (the trial positions of a CBMC step, the atoms of a
  /// molecule) before evaluating them overlaps those waits.
  void prefetch(double3 pos) const;

  /// The cell of a position: the lower-corner node indices and the position inside the cell (each in [0, 1)).
  struct Cell
  {
    std::size_t x0, y0, z0;
    double3 t;
  };
  Cell locateCell(double3 pos) const;

  /// Derivatives with respect to the cell coordinate converted to Cartesian space.
  double3 cartesianGradient(const double3& cellGradient) const;
  double3x3 cartesianHessian(const double3x3& cellHessian) const;

  void writeOutput(std::size_t systemId, const SimulationBox &simulationBox, const ForceField &forceField);

  static const uint3 parseExternalFieldGridDimensions(const std::string& filename);

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
