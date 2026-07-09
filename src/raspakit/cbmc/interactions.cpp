module;

module cbmc_interactions;

import std;

import atom;
import molecule;
import energy_factor;
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
import cbmc_interactions_external_field;
import cbmc_interactions_framework_molecule;
import cbmc_interactions_intermolecular;
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
    const GrowContext &context, const Component &component, std::vector<Atom> &trialPositions,
    std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<CBMC::FirstBeadTrial> energies{};
  energies.reserve(trialPositions.size());

  // loop over the trial-positions and compute the external energy of each trial position '{it, 1}'
  for (auto it = trialPositions.begin(); it != trialPositions.end(); ++it)
  {
    if (CBMC::insideBlockedPockets(context.framework, component, {it, 1}))
    {
      continue;
    }

    std::optional<RunningEnergy> externalFieldEnergy =
        CBMC::computeExternalFieldEnergy(context.hasExternalField, context.forceField, context.simulationBox,
                                         context.externalFieldInterpolationGrid, context.cutOffFrameworkVDW,
                                         context.cutOffCoulomb, {it, 1});

    // skip trial-positions that have an overlap in external-field energy
    if (!externalFieldEnergy.has_value()) continue;

    std::optional<RunningEnergy> frameworkEnergy = CBMC::computeFrameworkMoleculeEnergy(
        context.forceField, context.simulationBox, context.interpolationGrids, context.framework,
        context.frameworkAtoms, context.cutOffFrameworkVDW, context.cutOffCoulomb, {it, 1});

    // skip trial-positions that have an overlap in framework-molecule energy
    if (!frameworkEnergy.has_value()) continue;

    std::optional<RunningEnergy> interEnergy =
        CBMC::computeInterMolecularEnergy(context.forceField, context.simulationBox, context.moleculeAtoms,
                                          context.cutOffMoleculeVDW, context.cutOffCoulomb, {it, 1}, -1,
                                          skipBackgroundMolecule);

    // skip trial-positions that have an overlap in inter-molecular energy
    if (!interEnergy.has_value()) continue;

    // store position and energy
    energies.push_back({*it, externalFieldEnergy.value() + interEnergy.value() + frameworkEnergy.value()});
  }
  return energies;
}

std::vector<CBMC::ChainTrial> CBMC::computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::vector<std::vector<Atom>> &trialPositionSets,
    std::make_signed_t<std::size_t> skip, std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<CBMC::ChainTrial> energies{};
  energies.reserve(trialPositionSets.size());

  for (std::vector<Atom> trialPositionSet : trialPositionSets)
  {
    if (CBMC::insideBlockedPockets(context.framework, component, trialPositionSet))
    {
      continue;
    }

    std::optional<RunningEnergy> externalFieldEnergy =
        CBMC::computeExternalFieldEnergy(context.hasExternalField, context.forceField, context.simulationBox,
                                         context.externalFieldInterpolationGrid, context.cutOffFrameworkVDW,
                                         context.cutOffCoulomb, trialPositionSet, skip);
    if (!externalFieldEnergy.has_value()) continue;

    std::optional<RunningEnergy> frameworkEnergy = CBMC::computeFrameworkMoleculeEnergy(
        context.forceField, context.simulationBox, context.interpolationGrids, context.framework,
        context.frameworkAtoms, context.cutOffFrameworkVDW, context.cutOffCoulomb, trialPositionSet, skip);
    if (!frameworkEnergy.has_value()) continue;

    std::optional<RunningEnergy> interEnergy = CBMC::computeInterMolecularEnergy(
        context.forceField, context.simulationBox, context.moleculeAtoms, context.cutOffMoleculeVDW,
        context.cutOffCoulomb, trialPositionSet, skip, skipBackgroundMolecule);
    if (!interEnergy.has_value()) continue;

    energies.push_back(
        {trialPositionSet, externalFieldEnergy.value() + interEnergy.value() + frameworkEnergy.value()});
  }
  return energies;
}

std::vector<CBMC::ChainTrialTorsion> CBMC::computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component, std::vector<std::vector<Atom>> &trialPositionSets,
    const std::vector<double> &RosenbluthWeightsTorsion, std::make_signed_t<std::size_t> skip,
    std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<CBMC::ChainTrialTorsion> energies{};
  energies.reserve(trialPositionSets.size());

  for (std::size_t i = 0; i != trialPositionSets.size(); ++i)
  {
    if (CBMC::insideBlockedPockets(context.framework, component, trialPositionSets[i]))
    {
      continue;
    }

    std::optional<RunningEnergy> externalFieldEnergy =
        CBMC::computeExternalFieldEnergy(context.hasExternalField, context.forceField, context.simulationBox,
                                         context.externalFieldInterpolationGrid, context.cutOffFrameworkVDW,
                                         context.cutOffCoulomb, trialPositionSets[i], skip);
    if (!externalFieldEnergy.has_value())
    {
      continue;
    }

    std::optional<RunningEnergy> frameworkEnergy = CBMC::computeFrameworkMoleculeEnergy(
        context.forceField, context.simulationBox, context.interpolationGrids, context.framework,
        context.frameworkAtoms, context.cutOffFrameworkVDW, context.cutOffCoulomb, trialPositionSets[i], skip);
    if (!frameworkEnergy.has_value())
    {
      continue;
    }

    std::optional<RunningEnergy> interEnergy = CBMC::computeInterMolecularEnergy(
        context.forceField, context.simulationBox, context.moleculeAtoms, context.cutOffMoleculeVDW,
        context.cutOffCoulomb, trialPositionSets[i], skip, skipBackgroundMolecule);
    if (!interEnergy.has_value())
    {
      continue;
    }

    energies.push_back({trialPositionSets[i],
                        externalFieldEnergy.value() + frameworkEnergy.value() + interEnergy.value(),
                        RosenbluthWeightsTorsion[i]});
  }
  return energies;
}

std::vector<CBMC::MoleculeTrial> CBMC::computeExternalNonOverlappingEnergies(
    const GrowContext &context, const Component &component,
    std::vector<std::pair<Molecule, std::vector<Atom>>> &trialPositionSets, std::make_signed_t<std::size_t> skip,
    std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  std::vector<CBMC::MoleculeTrial> energies{};
  energies.reserve(trialPositionSets.size());

  for (auto &[molecule, trialPositionSet] : trialPositionSets)
  {
    if (CBMC::insideBlockedPockets(context.framework, component, trialPositionSet))
    {
      continue;
    }

    std::optional<RunningEnergy> externalFieldEnergy =
        CBMC::computeExternalFieldEnergy(context.hasExternalField, context.forceField, context.simulationBox,
                                         context.externalFieldInterpolationGrid, context.cutOffFrameworkVDW,
                                         context.cutOffCoulomb, trialPositionSet, skip);
    if (!externalFieldEnergy.has_value()) continue;

    std::optional<RunningEnergy> frameworkEnergy = CBMC::computeFrameworkMoleculeEnergy(
        context.forceField, context.simulationBox, context.interpolationGrids, context.framework,
        context.frameworkAtoms, context.cutOffFrameworkVDW, context.cutOffCoulomb, trialPositionSet, skip);
    if (!frameworkEnergy.has_value()) continue;

    std::optional<RunningEnergy> interEnergy = CBMC::computeInterMolecularEnergy(
        context.forceField, context.simulationBox, context.moleculeAtoms, context.cutOffMoleculeVDW,
        context.cutOffCoulomb, trialPositionSet, skip, skipBackgroundMolecule);
    if (!interEnergy.has_value()) continue;

    energies.push_back(
        {molecule, trialPositionSet, externalFieldEnergy.value() + interEnergy.value() + frameworkEnergy.value()});
  }
  return energies;
}

std::optional<RunningEnergy> CBMC::computeExternalNonOverlappingEnergyDualCutOff(
    const GrowContext &context, const Component &component, std::vector<Atom> &trialPositionSet,
    std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  if (CBMC::insideBlockedPockets(context.framework, component, trialPositionSet))
  {
    return std::nullopt;
  }

  std::optional<RunningEnergy> externalFieldEnergy =
      CBMC::computeExternalFieldEnergy(context.hasExternalField, context.forceField, context.simulationBox,
                                       context.externalFieldInterpolationGrid, context.cutOffFrameworkVDW,
                                       context.cutOffCoulomb, trialPositionSet, -1);
  if (!externalFieldEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkEnergy = CBMC::computeFrameworkMoleculeEnergy(
      context.forceField, context.simulationBox, context.interpolationGrids, context.framework,
      context.frameworkAtoms, context.cutOffFrameworkVDW, context.cutOffCoulomb, trialPositionSet, -1);
  if (!frameworkEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> interEnergy =
      CBMC::computeInterMolecularEnergy(context.forceField, context.simulationBox, context.moleculeAtoms,
                                        context.cutOffMoleculeVDW, context.cutOffCoulomb, trialPositionSet, -1,
                                        skipBackgroundMolecule);
  if (!interEnergy.has_value()) return std::nullopt;

  return externalFieldEnergy.value() + interEnergy.value() + frameworkEnergy.value();
}
