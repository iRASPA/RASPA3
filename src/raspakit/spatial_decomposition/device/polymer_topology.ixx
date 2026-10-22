module;

export module spatial_decomposition_device_polymer_topology;

import std;

import system;
import polymer_links;

/**
 * \brief The polymer-bond overlay of the Polymerization driver (PolymerLinkTable's derived term lists) in the
 * flattened form the device overlay kernels read (polymer_kernel_source.cpp), built on the host after every
 * reaction and uploaded by DevicePolymerOverlay.
 *
 * The terms are the device transcription of what Interactions::computePolymerLinkGradient evaluates: the link
 * bonds, bends and torsions (kinds 0..2), the overrides of the retyped monomer terms (3..5, 'replacement -
 * original'), the pair corrections (6, 'value(wanted) - value(counted)' with the Lennard-Jones parameters of the
 * two treatments baked in) and the stereocentre restraints (7). Every term carries the index of its ramp (one
 * per formation step of a bond or override), the ramp factors themselves are written per step by `lambdas`.
 * The gradient slots per participating atom are a CSR table, the inter-molecular slots of an atom first (the
 * virial correction applies to them only).
 *
 * One addition with no host counterpart: the device bonded kernels evaluate the scaled (1-4) pairs of a molecule
 * with the Lennard-Jones parameters of the component's pseudo-atom types, while the host bonded work follows the
 * current types of the atoms. For the molecules with a retyped atom the overlay therefore adds, as pair
 * corrections of the intramolecular kind, 'scaled(current types) - scaled(component types)' for the scaled pairs
 * that contain a retyped atom, so that the device total equals the host total.
 */
export struct PolymerOverlayTopology
{
  /// Mirrors the PolymerTerm struct of the kernel source.
  struct Term
  {
    std::uint32_t kind{0};
    std::uint32_t type{0};
    std::uint32_t type2{0};
    std::uint32_t ramp{0};
    std::uint32_t atoms[4]{};
    std::uint32_t flags{0};
    std::uint32_t padding[3]{};
    float parameters[6]{};
    float parameters2[6]{};
  };
  static_assert(sizeof(Term) == 96);

  static constexpr std::uint32_t interMolecular = 1u;
  static constexpr std::uint32_t intramolecularPair = 2u;
  /// The step of the ramp that is always fully on (the retype corrections).
  static constexpr std::uint64_t alwaysOn = std::numeric_limits<std::uint64_t>::max();

  /// Flattens the derived term lists of the table (which must be valid: PolymerLinkTable::rebuild) for the
  /// current atoms and force field of the system.
  void build(const System& system, const PolymerLinkTable& table);

  /// The ramp factor of every ramp at the table's current step (PolymerLinkTable::rampFactor), into `out`.
  void lambdas(const PolymerLinkTable& table, std::vector<float>& out) const;

  std::vector<Term> terms{};
  std::vector<std::uint32_t> atoms{};       ///< the participating atoms (system order indices), ascending
  std::vector<std::uint32_t> refStart{};    ///< CSR offsets of the gradient slots per participating atom
  std::vector<std::uint32_t> refLinkEnd{};  ///< per participating atom: end of its inter-molecular slots
  std::vector<std::uint32_t> refs{};        ///< gradient slots (4 * term + role)
  std::vector<std::uint64_t> rampSteps{};   ///< per ramp: the step the ramp started at
  std::size_t numberOfRetypeCorrections{0};
};
