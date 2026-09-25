module;

module cbmc_external_energy;

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
import cbmc_external_field;
import cbmc_framework_molecule;
import cbmc_intermolecular;
import cbmc_growth_context;

bool CBMC::insideBlockedPockets(const std::optional<Framework> &framework, const Component &component,
                                std::span<const Atom> molecule_atoms)
{
  if (framework.has_value())
  {
    for (std::size_t i = 0; i != component.blockingPockets.size(); ++i)
    {
      double radius_squared = component.blockingPockets[i].w * component.blockingPockets[i].w;
      double3 pos =
          framework->simulationBox.cell *
          double3(component.blockingPockets[i].x, component.blockingPockets[i].y, component.blockingPockets[i].z);
      for (const Atom &atom : molecule_atoms)
      {
        double3 dr = atom.position - pos;

        // compute the periodic boundary conditions with the single unit cell of the framework
        dr = framework->simulationBox.applyPeriodicBoundaryConditions(dr);

        double vdwScaling = atom.scalingVDW;
        if (dr.length_squared() < vdwScaling * radius_squared)
        {
          return true;
        }
      }
    }
  }
  return false;
}

[[nodiscard]] std::vector<CBMC::FirstBeadTrial> CBMC::computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::span<const Atom> trialPositions,
    std::optional<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<CBMC::FirstBeadTrial> energies{};
  energies.reserve(trialPositions.size());

  // Each first-bead trial is a one-atom trial set.
  for (const Atom &trialPosition : trialPositions)
  {
    std::optional<RunningEnergy> energy =
        computeExternalNonOverlappingEnergy(context, component, {&trialPosition, 1}, skipBackgroundMolecule);
    if (!energy.has_value()) continue;
    energies.push_back({trialPosition, energy.value()});
  }
  return energies;
}

std::optional<RunningEnergy> CBMC::computeExternalNonOverlappingEnergy(
    const GrowContext &context, const Component &component, std::span<const Atom> trialPositionSet,
    std::optional<std::size_t> skipBackgroundMolecule) noexcept
{
  if (CBMC::insideBlockedPockets(context.framework, component, trialPositionSet))
  {
    return std::nullopt;
  }

  std::optional<RunningEnergy> externalFieldEnergy =
      CBMC::computeExternalFieldEnergy(context.hasExternalField, context.forceField, context.simulationBox,
                                       context.externalFieldInterpolationGrid, context.cutOffFrameworkVDW,
                                       context.cutOffCoulomb, trialPositionSet);
  if (!externalFieldEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkEnergy = CBMC::computeFrameworkMoleculeEnergy(
      context.forceField, context.simulationBox, context.interpolationGrids, context.framework,
      context.frameworkAtoms, context.cutOffFrameworkVDW, context.cutOffCoulomb, trialPositionSet);
  if (!frameworkEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> interEnergy =
      CBMC::computeInterMolecularEnergy(context.forceField, context.simulationBox, context.moleculeAtoms,
                                        context.cutOffMoleculeVDW, context.cutOffCoulomb, trialPositionSet,
                                        skipBackgroundMolecule);
  if (!interEnergy.has_value()) return std::nullopt;

  return externalFieldEnergy.value() + interEnergy.value() + frameworkEnergy.value();
}

std::vector<CBMC::ChainTrialTorsion> CBMC::computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::vector<std::vector<Atom>> &trialPositionSets,
    const std::vector<double> &RosenbluthWeightsTorsion, std::optional<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<CBMC::ChainTrialTorsion> energies{};
  energies.reserve(trialPositionSets.size());

  for (std::size_t i = 0; i != trialPositionSets.size(); ++i)
  {
    std::optional<RunningEnergy> energy =
        computeExternalNonOverlappingEnergy(context, component, trialPositionSets[i], skipBackgroundMolecule);
    if (!energy.has_value()) continue;
    energies.push_back({std::move(trialPositionSets[i]), energy.value(), RosenbluthWeightsTorsion[i]});
  }
  return energies;
}

std::optional<RunningEnergy> CBMC::computeDualCutOffCorrection(const GrowContext &context, const Component &component,
                                                               std::span<const Atom> trialPositionSet,
                                                               std::optional<std::size_t> skipBackgroundMolecule) noexcept
{
  // The same configuration and background evaluated at the full and at the inner cut-offs.
  const GrowContext fullCutOffContext = context.withFullCutOffs();
  const GrowContext innerCutOffContext = context.withInnerCutOffs();

  std::optional<RunningEnergy> fullCutOffEnergy =
      CBMC::computeExternalNonOverlappingEnergy(fullCutOffContext, component, trialPositionSet, skipBackgroundMolecule);
  if (!fullCutOffEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> innerCutOffEnergy = CBMC::computeExternalNonOverlappingEnergy(
      innerCutOffContext, component, trialPositionSet, skipBackgroundMolecule);
  if (!innerCutOffEnergy.has_value()) return std::nullopt;

  return fullCutOffEnergy.value() - innerCutOffEnergy.value();
}
