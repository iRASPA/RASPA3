module;

export module spatial_decomposition_device_bonded_topology;

import std;

import system;

/**
 * \brief The per-molecule terms of a system in the flattened form the device bonded kernels read, built on the
 * host once per system and uploaded by the backend.
 *
 * The terms of every component (bonds, bends, torsions, improper torsions, intramolecular Lennard-Jones and
 * Coulomb pairs) are flattened and grouped by kind, with the gradient slots per component atom in CSR form; the
 * molecules refer to the block of their component and get consecutive term instances and gradient slots. The
 * struct layouts mirror the Term and MoleculeInfo structs of bonded_kernel_source.cpp.
 *
 * `supports` reports the components that use other intramolecular terms (Urey-Bradley, inversion bends, cross
 * terms), molecules with more than 256 atoms and lambda-group atoms, for which the engine keeps the bonded work
 * on the host.
 */
export struct BondedTopology
{
  /// Mirrors the Term struct of the kernel source.
  struct Term
  {
    std::uint32_t kind{0};  ///< 0 bond, 1 bend, 2 torsion, 3 improper torsion, 4 intramolecular LJ, 5 Coulomb
    std::uint32_t type{0};  ///< the BondType / BendType / TorsionType value
    std::uint32_t atoms[4]{};
    float parameters[6]{};
  };
  static_assert(sizeof(Term) == 48);

  /// Mirrors the MoleculeInfo struct of the kernel source.
  struct MoleculeInfo
  {
    std::uint32_t firstAtom{0}, numberOfAtoms{0}, atomOffset{0}, termOffset{0}, instanceBase{0}, gradientBase{0};
    std::uint32_t padding[2]{};
  };
  static_assert(sizeof(MoleculeInfo) == 32);

  static constexpr std::size_t maximumAtomsPerMolecule = 256;
  static constexpr std::uint32_t noAtom = std::numeric_limits<std::uint32_t>::max();

  /// Whether the device kernels cover the intramolecular terms of the system; on false, `reason` names the first
  /// unsupported feature.
  static bool supports(const System& system, std::string& reason);

  /// Flattens the terms, the molecule table and the masses of the system; throws when the 32-bit instance or
  /// gradient indices would overflow.
  void build(const System& system);

  /// The slot layout of a list build: per slot (molecule << 8) | index in the molecule (noAtom for a dummy
  /// slot), and per sorted atom the sorted index of the first atom of its molecule (the reference of the
  /// relative positions the host packs).
  void layout(std::span<const std::uint32_t> slotOfSorted, std::span<const std::uint32_t> originalToSorted,
              std::size_t slots, std::vector<std::uint32_t>& slotMolecule,
              std::vector<std::uint32_t>& referenceOfSorted) const;

  std::vector<Term> terms{};                        ///< the terms of all components, by component then kind
  std::vector<std::uint32_t> gradientOffset{};      ///< per term: offset of its gradient slots in the component block
  std::vector<std::uint32_t> atomGradientStart{};   ///< CSR offsets of the gradient slots per component atom
  std::vector<std::uint32_t> atomGradients{};       ///< gradient slot (within the molecule's block)
  std::vector<std::uint32_t> instanceMolecule{};    ///< molecule of every term instance
  std::vector<MoleculeInfo> molecules{};
  std::vector<float> massOfAtom{};                  ///< per atom in the system order
  std::size_t numberOfInstances{0};                 ///< term instances over all molecules
  std::size_t numberOfGradients{0};                 ///< gradient slots over all molecules
  double chargeSquaredSum{0.0};                     ///< sum q^2 over the atoms (the self energy on the host)
};
