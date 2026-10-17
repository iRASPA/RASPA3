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
import mc_cell_list;
import interactions_neighbourhood;
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

/// Framework-molecule energy of a trial set in the context's environment: the atoms with a grid read it as
/// above; the others are evaluated against the context's framework neighbourhood when it covers the set, or
/// against the framework atoms gathered around the set through the context's framework cell list (see
/// GrowContext::frameworkCellList), and by the brute-force loop when the context has no such list.
[[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy(const GrowContext& context,
                                                                          std::span<const Atom> atoms) noexcept;

/// Inter-molecular energy of the trial atoms 'atoms' against the background 'moleculeAtoms' (pairs
/// within the same molecule id are skipped, as are all background atoms of 'skipBackgroundMolecule'
/// when given -- the molecule being regrown). std::nullopt on a hard overlap.
[[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    double cutOffVDW, double cutOffCoulomb, std::span<const Atom> atoms,
    std::optional<std::size_t> skipBackgroundMolecule = std::nullopt) noexcept;

/// Cell-list variant of 'computeInterMolecularEnergy': the same result (to the summation order), gathering
/// the background within the cut-off of the trial set's bounding sphere through the cell list into a
/// thread-local Interactions::Neighbourhood and evaluating with its dense kernel. 'cellList' must be
/// built over 'moleculeAtoms'; with a disabled or stale list the neighbourhood is gathered by a pass over
/// all atoms instead.
[[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, const MCCellList &cellList,
    std::span<const Atom> moleculeAtoms, double cutOffVDW, double cutOffCoulomb, std::span<const Atom> atoms,
    std::optional<std::size_t> skipBackgroundMolecule = std::nullopt) noexcept;

/// Inter-molecular energy of a trial set against a pre-gathered neighbourhood (which must cover the set
/// and the cut-offs, see 'Neighbourhood::covers'); the same result as the overloads above.
[[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergy(const ForceField& forceField,
                                                                       const SimulationBox& simulationBox,
                                                                       const Interactions::Neighbourhood& neighbourhood,
                                                                       double cutOffVDW, double cutOffCoulomb,
                                                                       std::span<const Atom> atoms) noexcept;

/// Gathers the background of a batch of trial sets (the context's molecule atoms minus its skipped
/// molecule, within the force field's full cut-offs of the batch's bounding sphere) into 'neighbourhood',
/// through the context's cell list when it has one. The batch is then evaluated with
/// 'context.withNeighbourhood(&neighbourhood)'. 'extraRadius' widens the bounding sphere, so that the
/// neighbourhood also covers trial sets near the batch that are generated later (see Neighbourhood::covers).
void gatherNeighbourhood(const GrowContext& context, Interactions::Neighbourhood& neighbourhood,
                         std::span<const std::span<const Atom>> trialSets, double extraRadius = 0.0) noexcept;

/// The framework counterpart of 'gatherNeighbourhood': the context's framework atoms within the full cut-offs
/// of the batch's bounding sphere, through its framework cell list. Returns false (and gathers nothing) when
/// the batch needs no framework pairs -- no framework cell list in the context, no framework atoms, or every
/// trial atom reads an interpolation grid. On true the batch is evaluated with
/// 'context.withFrameworkNeighbourhood(&neighbourhood)'.
[[nodiscard]] bool gatherFrameworkNeighbourhood(const GrowContext& context, Interactions::Neighbourhood& neighbourhood,
                                                std::span<const std::span<const Atom>> trialSets,
                                                double extraRadius = 0.0) noexcept;

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
