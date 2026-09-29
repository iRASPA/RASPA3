module;

export module spatial_decomposition_settings;

import std;

import int3;
import archive;

/**
 * \brief User-facing settings of the spatial-decomposition MD force engine.
 *
 * Read from the input file by the input reader for 'SimulationType' : 'MolecularDynamicsSpatialDecomposition'
 * and handed to the driver, which constructs one engine per system from them. The engine itself holds only
 * derived state (cells, neighbour lists, mesh, FFT plans) and is rebuilt from these settings on a restart.
 */
export struct SpatialDecompositionSettings
{
  std::uint64_t versionNumber{1};

  /// Number of worker threads (1: every phase runs on the calling thread through the same code path).
  std::size_t numberOfThreads{1};

  /// Verlet skin [Angstrom] added to the cutoff for the neighbour lists; the lists are rebuilt when an atom
  /// has moved more than half the skin since the last build. Zero rebuilds every step.
  double verletSkin{2.0};

  /// Target mesh spacing [Angstrom] of the particle-mesh Ewald sum; the mesh dimensions are the smallest
  /// FFT-friendly sizes (2^a 3^b 5^c) at or below this spacing.
  double meshSpacing{1.0};

  /// Order of the cardinal B-spline charge assignment (3 to 7; the standard choice is 5).
  std::size_t interpolationOrder{5};

  /// Optional explicit sub-domain grid (nx, ny, nz); the product must equal the number of threads. When not
  /// given the thread count is factorized into the grid with the smallest total sub-domain surface.
  std::optional<int3> domainGrid{};

  bool operator==(const SpatialDecompositionSettings&) const = default;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const SpatialDecompositionSettings& s);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, SpatialDecompositionSettings& s);
};

export Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const SpatialDecompositionSettings& s)
{
  archive << s.versionNumber;
  archive << s.numberOfThreads;
  archive << s.verletSkin;
  archive << s.meshSpacing;
  archive << s.interpolationOrder;
  archive << s.domainGrid.has_value();
  if (s.domainGrid.has_value())
  {
    archive << s.domainGrid->x;
    archive << s.domainGrid->y;
    archive << s.domainGrid->z;
  }
  return archive;
}

export Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, SpatialDecompositionSettings& s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > s.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'SpatialDecompositionSettings' at line {} in file {}\n", location.line(),
                    location.file_name()));
  }
  archive >> s.numberOfThreads;
  archive >> s.verletSkin;
  archive >> s.meshSpacing;
  archive >> s.interpolationOrder;
  bool hasGrid;
  archive >> hasGrid;
  if (hasGrid)
  {
    int3 grid;
    archive >> grid.x;
    archive >> grid.y;
    archive >> grid.z;
    s.domainGrid = grid;
  }
  else
  {
    s.domainGrid.reset();
  }
  return archive;
}
