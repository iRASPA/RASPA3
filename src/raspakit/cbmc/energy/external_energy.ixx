module;

export module cbmc_external_energy;

import std;

import atom;
import running_energy;
import framework;
import component;
import forcefield;
import simulationbox;
import interpolation_energy_grid;
import cbmc_grow_context;
import cbmc_results;

// The external energy of a CBMC trial set: external field, framework-molecule, and inter-molecular,
// each with an overlap short-circuit, combined by 'computeExternalNonOverlappingEnergy', plus the dual
// cut-off correction built on it. The pair evaluations are the shared kernel of the 'interactions'
// module ('Interactions::evaluatePair'), so the energies a molecule is grown with are by construction
// those the system accounts for; only the explicit cut-offs (dual cut-off scheme), the overlap
// short-circuit, the background exclusions, and the interpolation-grid path are CBMC's own.
export namespace CBMC
{
/// Whether any of 'molecule_atoms' lies inside one of the component's blocking pockets (spheres in
/// fractional framework coordinates, radius scaled by the atom's 'scalingVDW'). Every trial set of a
/// CBMC grow is filtered with this, so a grown molecule never lies in a pocket; 'System' delegates its
/// own check (for non-CBMC placements and rescaled fractional molecules) to this function.
bool insideBlockedPockets(const std::optional<Framework> &framework, const Component &component,
                          std::span<const Atom> molecule_atoms);

/// External-field energy of 'atoms' (the force field's potential energy surface or its interpolation
/// grid); std::nullopt when an atom lies outside a confining geometry.
[[nodiscard]] std::optional<RunningEnergy> computeExternalFieldEnergy(
    bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid, double cutOffVDW,
    double cutOffCoulomb, std::span<const Atom> atoms) noexcept;

/// Framework-molecule energy of 'atoms' at the given cut-offs, from the interpolation grids where one
/// exists for the atom type (never for a fractional atom) and by explicit pair summation otherwise;
/// std::nullopt on a hard overlap.
[[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms, double cutOffVDW,
    double cutOffCoulomb, std::span<const Atom> atoms) noexcept;

/// Inter-molecular energy of the trial atoms 'atoms' against the background 'moleculeAtoms' (pairs
/// within the same molecule id are skipped, as are all background atoms of 'skipBackgroundMolecule'
/// when given -- the molecule being regrown). std::nullopt on a hard overlap.
[[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    double cutOffVDW, double cutOffCoulomb, std::span<const Atom> atoms,
    std::optional<std::size_t> skipBackgroundMolecule = std::nullopt) noexcept;

// The background a trial set is evaluated against is the context's: its molecule atoms minus those
// with the trial atoms' own molecule id and minus 'context.skipBackgroundMolecule' (see GrowContext).

/// External energies of single first-bead trial positions; positions in a blocked pocket or with an
/// overlap are dropped from the result.
[[nodiscard]] std::vector<FirstBeadTrial> computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::span<const Atom> trialPositions) noexcept;

/// External (external-field, framework-molecule, inter-molecular) energy of one trial set of bead
/// positions at the context's cut-offs; std::nullopt when the set lies in a blocked pocket or
/// overlaps. The single evaluation everything else is built on: the first-bead overload above, the
/// dual cut-off correction, and both chain schemes (one trial set at a time, no intermediate
/// containers).
[[nodiscard]] std::optional<RunningEnergy> computeExternalNonOverlappingEnergy(
    const GrowContext &context, const Component &component, std::span<const Atom> trialPositionSet) noexcept;

/// Dual cut-off correction of a grown or retraced configuration: the external (external-field,
/// framework-molecule, and inter-molecular) energy of 'trialPositionSet' evaluated at the full
/// cut-offs minus the same energy evaluated at the inner cut-off used during growth. The passed
/// context supplies the background (framework and molecule atoms); its cut-offs are ignored. The
/// entry points of the 'cbmc' module fold it into their results (energies += correction,
/// multiplyRosenbluthWeight(-beta * correction)), so a configuration behaves as if it had been grown
/// at the full cut-offs. Returns std::nullopt when the configuration overlaps at the full cut-offs.
[[nodiscard]] std::optional<RunningEnergy> computeDualCutOffCorrection(
    const GrowContext &context, const Component &component, std::span<const Atom> trialPositionSet) noexcept;
}  // namespace CBMC
