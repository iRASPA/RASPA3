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
/// Arithmetic precision of the specialised pair kernel. `Double` evaluates the pairs in double (the scalar
/// kernel over the half lists); `Mixed` runs the SIMD cluster kernel: the per-pair arithmetic (positions relative
/// to the sub-domain, distances, Lennard-Jones and tabulated Ewald terms, the forces within one cluster block) in
/// single precision, the per-atom forces, energies and the virial accumulated in double. The lists, the mesh,
/// the bonded terms and the integration are always double.
export enum class PairPrecision : std::uint8_t
{
  Double = 0,
  Mixed = 1
};

export inline std::string pairPrecisionName(PairPrecision precision)
{
  return precision == PairPrecision::Mixed ? "Mixed" : "Double";
}

/// Where the specialised pair kernel runs. `CPU`: the SIMD cluster kernel (or the scalar kernel) on the worker
/// threads. `OpenCL` / `Metal`: the cluster pair kernel on the device (GPU) in single precision, overlapped with
/// the mesh, bonded and exclusion work of the worker threads; requires a device of the kind (Metal: macOS builds)
/// and the specialised kernel (Lennard-Jones with Ewald or no electrostatics, fully coupled atoms).
export enum class PairDevice : std::uint8_t
{
  CPU = 0,
  OpenCL = 1,
  Metal = 2
};

export inline std::string pairDeviceName(PairDevice device)
{
  switch (device)
  {
    case PairDevice::OpenCL:
      return "OpenCL";
    case PairDevice::Metal:
      return "Metal";
    case PairDevice::CPU:
      break;
  }
  return "CPU";
}
/// Whether the device is a GPU backend (the device step), not the CPU.
export inline bool isDevice(PairDevice device) { return device != PairDevice::CPU; }

export struct SpatialDecompositionSettings
{
  std::uint64_t versionNumber{5};

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

  /// Precision of the specialised pair kernel (see PairPrecision).
  PairPrecision pairPrecision{PairPrecision::Double};

  /// Skin [Angstrom] of the pruned (inner) pair list of the cluster kernel: every few steps the blocks of the
  /// Verlet list are pruned to the pairs within cutoff + pruneSkin, and the pruned list is reused while no atom
  /// has moved more than half this skin. Must be below the Verlet skin; 0 disables pruning.
  double pruneSkin{0.5};

  /// Runs the double-precision instantiation of the cluster kernel for `PairPrecision::Double` instead of the
  /// scalar kernel. Not an input option: the two give the same numbers to rounding (used by the tests to validate
  /// the cluster code path), and on 128-bit SIMD (NEON) the scalar kernel is the faster of the two.
  bool clusterKernelForDouble{false};

  /// Device of the specialised pair kernel (see PairDevice).
  PairDevice pairDevice{PairDevice::CPU};

  /// With a device pair kernel: the particle-mesh Ewald sum runs on the device as well (charge spreading, FFTs,
  /// influence function, interpolation). Not an input option (the tests and benchmarks switch it off to compare
  /// the device pairs with the host mesh).
  bool deviceMesh{true};

  /// With a device pair kernel: the bonded terms and the self / exclusion corrections run on the device when
  /// the device kernels cover the intramolecular potentials of the system (else they stay on the host). Not an
  /// input option.
  bool deviceBonded{true};

  /// With a device pair kernel, mesh and bonded terms: the MD state (positions, velocities, molecule records)
  /// stays on the device and the velocity-Verlet integrator runs there in double-float arithmetic (hi + lo
  /// floats, "emulated double"); the host keeps the thermostat and the list rebuilds and downloads the state
  /// only when it samples properties or writes a restart file. Falls back to the host integrator when the system
  /// is not covered (semi-flexible molecules, barostats, bonded terms on the host). Input option 'Resident'.
  bool resident{true};

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
  archive << static_cast<std::uint8_t>(s.pairPrecision);
  archive << s.pruneSkin;
  archive << s.clusterKernelForDouble;
  archive << static_cast<std::uint8_t>(s.pairDevice);
  archive << s.deviceMesh;
  archive << s.deviceBonded;
  archive << s.resident;
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
  s.pairPrecision = PairPrecision::Double;
  s.pruneSkin = 0.5;
  s.clusterKernelForDouble = false;
  if (versionNumber >= 2)
  {
    std::uint8_t precision;
    archive >> precision;
    s.pairPrecision = precision == 1 ? PairPrecision::Mixed : PairPrecision::Double;
    archive >> s.pruneSkin;
    archive >> s.clusterKernelForDouble;
  }
  s.pairDevice = PairDevice::CPU;
  if (versionNumber >= 3)
  {
    std::uint8_t device;
    archive >> device;
    s.pairDevice = device == 1 ? PairDevice::OpenCL : device == 2 ? PairDevice::Metal : PairDevice::CPU;
  }
  s.deviceMesh = true;
  s.deviceBonded = true;
  if (versionNumber >= 4)
  {
    archive >> s.deviceMesh;
    archive >> s.deviceBonded;
  }
  s.resident = true;
  if (versionNumber >= 5)
  {
    archive >> s.resident;
  }
  return archive;
}
