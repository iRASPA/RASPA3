module;

export module cbmc_external_energy;

import std;

import atom;
import molecule;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import framework;
import component;
import double3;
import double3x3;
import forcefield;
import simulationbox;
import units;
import cbmc_intermolecular;
import cbmc_framework_molecule;
import interpolation_energy_grid;
import cbmc_growth_context;

export namespace CBMC
{
/// External energy of a single first-bead trial position.
struct FirstBeadTrial
{
  Atom position;
  RunningEnergy energy;
};

/// Whether any of 'molecule_atoms' lies inside one of the component's blocking pockets (spheres in
/// fractional framework coordinates, radius scaled by the atom's 'scalingVDW'). Every trial set of a
/// CBMC grow is filtered with this, so a grown molecule never lies in a pocket; 'System' delegates its
/// own check (for non-CBMC placements and rescaled fractional molecules) to this function.
bool insideBlockedPockets(const std::optional<Framework> &framework, const Component &component,
                          std::span<const Atom> molecule_atoms);

// 'skipBackgroundMolecule', where present, is the molecule id whose atoms in the context's background
// are ignored: the molecule being regrown (reinsertion, identity change), which is still present in
// the background but must not interact with its own trial positions. std::nullopt skips nothing.

/// External energies of single first-bead trial positions; positions in a blocked pocket or with an
/// overlap are dropped from the result.
[[nodiscard]] std::vector<FirstBeadTrial> computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::span<const Atom> trialPositions,
    std::optional<std::size_t> skipBackgroundMolecule = std::nullopt) noexcept;

/// External (external-field, framework-molecule, inter-molecular) energy of one trial set of bead
/// positions at the context's cut-offs; std::nullopt when the set lies in a blocked pocket or
/// overlaps. The single evaluation everything else is built on: the first-bead overload above, the
/// dual cut-off correction, and both chain schemes (one trial set at a time, no intermediate
/// containers).
[[nodiscard]] std::optional<RunningEnergy> computeExternalNonOverlappingEnergy(
    const GrowContext &context, const Component &component, std::span<const Atom> trialPositionSet,
    std::optional<std::size_t> skipBackgroundMolecule = std::nullopt) noexcept;

/// Dual cut-off correction of a grown or retraced configuration: the external (external-field,
/// framework-molecule, and inter-molecular) energy of 'trialPositionSet' evaluated at the full
/// cut-offs minus the same energy evaluated at the inner cut-off used during growth. The passed
/// context supplies the background (framework and molecule atoms); its cut-offs are ignored.
/// Folding the correction into the growth or retrace results (energies += correction,
/// multiplyRosenbluthWeight(-beta * correction)) makes the configuration behave as if it had been
/// grown at the full cut-offs. Returns std::nullopt when the configuration overlaps at the full
/// cut-offs.
[[nodiscard]] std::optional<RunningEnergy> computeDualCutOffCorrection(
    const GrowContext &context, const Component &component, std::span<const Atom> trialPositionSet,
    std::optional<std::size_t> skipBackgroundMolecule = std::nullopt) noexcept;
}  // namespace CBMC
