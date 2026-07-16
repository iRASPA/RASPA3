module;

export module cbmc_interactions;

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
import cbmc_interactions_intermolecular;
import cbmc_interactions_framework_molecule;
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

/// External energy of a trial set of bead positions.
struct ChainTrial
{
  std::vector<Atom> positions;
  RunningEnergy energy;
};

/// External energy of a trial set of bead positions, together with the torsion Rosenbluth weight
/// accumulated while generating that orientation.
struct ChainTrialTorsion
{
  std::vector<Atom> positions;
  RunningEnergy energy;
  double torsionWeight;
};

/// External energy of a trial rigid-molecule orientation (molecule record + its atoms).
struct MoleculeTrial
{
  Molecule molecule;
  std::vector<Atom> positions;
  RunningEnergy energy;
};

bool insideBlockedPockets(const std::optional<Framework> &frameworks, const Component &component,
                          std::span<const Atom> molecule_atoms);

[[nodiscard]] std::vector<FirstBeadTrial> computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::vector<Atom> &trialPositions,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

[[nodiscard]] std::vector<ChainTrial> computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::vector<std::vector<Atom>> &trialPositionSets,
    std::make_signed_t<std::size_t> skip = -1, std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

[[nodiscard]] std::vector<ChainTrialTorsion> computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::vector<std::vector<Atom>> &trialPositionSets,
    const std::vector<double> &RosenbluthWeightsTorsion, std::make_signed_t<std::size_t> skip = -1,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

[[nodiscard]] std::vector<MoleculeTrial> computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component,
    std::vector<std::pair<Molecule, std::vector<Atom>>> &trialPositionSets, std::make_signed_t<std::size_t> skip = -1,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

[[nodiscard]] std::optional<RunningEnergy> computeExternalNonOverlappingEnergyDualCutOff(
    const GrowContext &context, const Component &component, std::vector<Atom> &trialPositionSet,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

/// Dual cut-off correction of a grown or retraced configuration: the external (external-field,
/// framework-molecule, and inter-molecular) energy of 'trialPositionSet' evaluated at the full
/// cut-offs minus the same energy evaluated at the inner cut-off used during growth. The passed
/// context supplies the background (framework and molecule atoms); its cut-offs are ignored.
/// Folding the correction into the growth or retrace results (energies += correction,
/// RosenbluthWeight *= exp(-beta * correction)) makes the configuration behave as if it had been
/// grown at the full cut-offs. Returns std::nullopt when the configuration overlaps at the full
/// cut-offs.
[[nodiscard]] std::optional<RunningEnergy> computeDualCutOffCorrection(
    const GrowContext &context, const Component &component, std::vector<Atom> &trialPositionSet,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;
}  // namespace CBMC
