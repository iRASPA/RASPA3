module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <utility>
#include <vector>
#endif

module cbmc;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import randomnumbers;
import component;
import atom;
import molecule;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import energy_status;
import forcefield;
import energy_factor;
import cbmc_first_bead_data;
import cbmc_multiple_first_bead;
import cbmc_rigid_insertion;
import cbmc_rigid_deletion;
import cbmc_flexible_insertion;
import cbmc_flexible_deletion;
import framework;
import component;
import cbmc_chain_data;
import interpolation_energy_grid;

// Insertion:
// Insertion means growing a new molecule, and therefore the atrributes of the atoms are
// taken from 'component.atoms'. The parameters 'scaling', 'groupId',
// and 'isFractional' are passed to the function and set on the atoms.

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculeSwapInsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::size_t selectedMolecule, double scaling, bool groupId,
    bool isFractional) noexcept
{
  std::size_t startingBead = component.startingBead;
  Atom firstBead = component.atoms[startingBead];
  firstBead.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
  firstBead.groupId = groupId;
  firstBead.isFractional = isFractional;
  firstBead.setScaling(scaling);

  std::optional<FirstBeadData> const firstBeadData = CBMC::growMoleculeMultipleFirstBeadSwapInsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtomData,
      moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, firstBead);

  if (!firstBeadData) return std::nullopt;

  if (component.atoms.size() == 1)
  {
    return ChainGrowData(Molecule(double3(firstBeadData->atom.position), simd_quatd(0.0, 0.0, 0.0, 1.0),
                                  component.totalMass, component.componentId, component.definedAtoms.size()),
                         {firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  // place the molecule centered around the first bead at 'firstBeadData->atom.position'
  std::vector<Atom> molecule_atoms = component.atoms;
  for (Atom &atom : molecule_atoms)
  {
    atom.position += firstBeadData->atom.position - component.atoms[startingBead].position;
    atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
    atom.groupId = groupId;
    atom.isFractional = isFractional;
    atom.setScaling(scaling);
  }

  std::optional<ChainGrowData> chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::growRigidMoleculeChainInsertion(random, component, hasExternalField, forceField, simulationBox,
                                                        interpolationGrids, framework, frameworkAtomData,
                                                        moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW,
                                                        cutOffCoulomb, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = CBMC::growFlexibleMoleculeChainInsertion(
          random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
          molecule_atoms, {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  if (!chainData) return std::nullopt;

  return ChainGrowData(chainData->molecule, chainData->atom, firstBeadData->energies + chainData->energies,
                       firstBeadData->RosenbluthWeight * chainData->RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculeSwapDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::span<Atom> molecule_atoms) noexcept
{
  std::size_t startingBead = component.startingBead;

  const FirstBeadData firstBeadData = CBMC::retraceMultipleFirstBeadSwapDeletion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtomData,
      moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, molecule_atoms[startingBead]);

  if (molecule_atoms.size() == 1)
  {
    return ChainRetraceData(firstBeadData.energies, firstBeadData.RosenbluthWeight, 0.0);
  }

  ChainRetraceData chainData;
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::retraceRigidMoleculeChainDeletion(random, component, hasExternalField, forceField,
                                                          simulationBox, interpolationGrids, framework,
                                                          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW,
                                                          cutOffMoleculeVDW, cutOffCoulomb, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = CBMC::retraceFlexibleMoleculeChainDeletion(
          random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
          molecule_atoms, {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  return ChainRetraceData(firstBeadData.energies + chainData.energies,
                          firstBeadData.RosenbluthWeight * chainData.RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculeReinsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, Molecule &molecule, std::span<Atom> molecule_atoms) noexcept
{
  std::size_t startingBead = component.startingBead;

  std::optional<FirstBeadData> const firstBeadData = CBMC::growMultipleFirstBeadReinsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtomData,
      moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, molecule_atoms[startingBead]);

  if (!firstBeadData) return std::nullopt;

  // place the molecule centered around the first bead at 'firstBeadData->atom.position'
  std::vector<Atom> atoms = component.atoms;
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    atoms[i].position += firstBeadData->atom.position - component.atoms[startingBead].position;
    atoms[i].charge = molecule_atoms[i].charge;
    atoms[i].scalingVDW = molecule_atoms[i].scalingVDW;
    atoms[i].scalingCoulomb = molecule_atoms[i].scalingCoulomb;
    atoms[i].moleculeId = molecule_atoms[i].moleculeId;
    atoms[i].componentId = molecule_atoms[i].componentId;
    atoms[i].groupId = molecule_atoms[i].groupId;
    atoms[i].isFractional = molecule_atoms[i].isFractional;
  }

  if (molecule_atoms.size() == 1)
  {
    Molecule firstBeadMolecule = Molecule(firstBeadData->atom.position, simd_quatd(), component.totalMass,
                                          component.componentId, component.definedAtoms.size());
    firstBeadMolecule.atomIndex = molecule.atomIndex;
    firstBeadMolecule.numberOfAtoms = molecule.numberOfAtoms;

    return ChainGrowData(firstBeadMolecule, {firstBeadData->atom}, firstBeadData->energies,
                         firstBeadData->RosenbluthWeight, firstBeadData->storedR);
  }

  std::optional<ChainGrowData> chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::growRigidMoleculeChainInsertion(
          random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = CBMC::growFlexibleMoleculeChainInsertion(
          random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, atoms,
          {component.startingBead});
      break;
    default:
      std::unreachable();
  }
  if (!chainData) return std::nullopt;

  // Copy data over from old molecule, since we are reinserting the same molecule
  chainData->molecule.atomIndex = molecule.atomIndex;
  chainData->molecule.numberOfAtoms = molecule.numberOfAtoms;

  return ChainGrowData(chainData->molecule, chainData->atom, firstBeadData->energies + chainData->energies,
                       firstBeadData->RosenbluthWeight * chainData->RosenbluthWeight, firstBeadData->storedR);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, [[maybe_unused]] Molecule &molecule, std::span<Atom> molecule_atoms,
    double storedR) noexcept
{
  std::size_t startingBead = component.startingBead;

  const FirstBeadData firstBeadData = CBMC::retraceMultipleFirstBeadReinsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtomData,
      moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, molecule_atoms[startingBead],
      storedR);

  if (molecule_atoms.size() == 1)
  {
    return ChainRetraceData(firstBeadData.energies, firstBeadData.RosenbluthWeight, 0.0);
  }

  ChainRetraceData chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::retraceRigidMoleculeChainDeletion(random, component, hasExternalField, forceField,
                                                          simulationBox, interpolationGrids, framework,
                                                          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW,
                                                          cutOffMoleculeVDW, cutOffCoulomb, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = CBMC::retraceFlexibleMoleculeChainDeletion(
          random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
          molecule_atoms, {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  return ChainRetraceData(firstBeadData.energies + chainData.energies,
                          firstBeadData.RosenbluthWeight * chainData.RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculePartialReinsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, Molecule &molecule, std::span<Atom> moleculeAtoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced) noexcept
{
  std::optional<ChainGrowData> chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::growRigidMoleculeChainInsertion(random, component, hasExternalField, forceField, simulationBox,
                                                        interpolationGrids, framework, frameworkAtomData,
                                                        moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW,
                                                        cutOffCoulomb, moleculeAtoms);
      break;
    case Component::GrowType::Flexible:
      chainData = CBMC::growFlexibleMoleculeChainInsertion(
          random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
          moleculeAtoms, beadsAlreadyPlaced);
      break;
    default:
      std::unreachable();
  }

  if (!chainData) return std::nullopt;

  // Copy data over from old molecule, since we are reinserting the same molecule
  chainData->molecule.atomIndex = molecule.atomIndex;
  chainData->molecule.numberOfAtoms = molecule.numberOfAtoms;

  return ChainGrowData(chainData->molecule, chainData->atom, chainData->energies, chainData->RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculePartialReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, [[maybe_unused]] Molecule &molecule, std::span<Atom> moleculeAtoms,
    const std::vector<std::size_t> &beadsAlreadyPlaced) noexcept
{
  ChainRetraceData chainData;
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::retraceRigidMoleculeChainDeletion(random, component, hasExternalField, forceField,
                                                          simulationBox, interpolationGrids, framework,
                                                          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW,
                                                          cutOffMoleculeVDW, cutOffCoulomb, moleculeAtoms);
      break;
    case Component::GrowType::Flexible:
      chainData = CBMC::retraceFlexibleMoleculeChainDeletion(
          random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
          moleculeAtoms, beadsAlreadyPlaced);
      break;
    default:
      std::unreachable();
  }

  return ChainRetraceData(chainData.energies, chainData.RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<ChainGrowData> CBMC::growMoleculeIdentityChangeInsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::size_t selectedMolecule, double scaling, bool groupId,
    bool isFractional) noexcept
{
  std::size_t startingBead = component.startingBead;
  Atom firstBead = component.atoms[startingBead];
  firstBead.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
  firstBead.groupId = groupId;
  firstBead.isFractional = isFractional;
  firstBead.setScaling(scaling);

  std::optional<FirstBeadData> const firstBeadData = CBMC::growMoleculeMultipleFirstBeadSwapInsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtomData,
      moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, firstBead);

  if (!firstBeadData) return std::nullopt;

  if (component.atoms.size() == 1)
  {
    // update atom index
    return ChainGrowData(Molecule(double3(firstBeadData->atom.position), simd_quatd(0.0, 0.0, 0.0, 1.0),
                                  component.totalMass, component.componentId, component.definedAtoms.size()),
                         {firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  // place the molecule centered around the first bead at 'firstBeadData->atom.position'
  std::vector<Atom> molecule_atoms = component.atoms;
  for (Atom &atom : molecule_atoms)
  {
    atom.position += firstBeadData->atom.position - component.atoms[startingBead].position;
    atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
    atom.groupId = groupId;
    atom.isFractional = isFractional;
    atom.setScaling(scaling);
  }

  std::optional<ChainGrowData> chainData{};
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::growRigidMoleculeChainInsertion(random, component, hasExternalField, forceField, simulationBox,
                                                        interpolationGrids, framework, frameworkAtomData,
                                                        moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW,
                                                        cutOffCoulomb, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = CBMC::growFlexibleMoleculeChainInsertion(
          random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
          molecule_atoms, {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  if (!chainData) return std::nullopt;
  // update atom index

  return ChainGrowData(chainData->molecule, chainData->atom, firstBeadData->energies + chainData->energies,
                       firstBeadData->RosenbluthWeight * chainData->RosenbluthWeight, 0.0);
}

[[nodiscard]] ChainRetraceData CBMC::retraceMoleculeIdentityChangeDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtomData,
    std::span<const Atom> moleculeAtomData, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::span<Atom> molecule_atoms) noexcept
{
  std::size_t startingBead = component.startingBead;

  const FirstBeadData firstBeadData = CBMC::retraceMultipleFirstBeadSwapDeletion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtomData,
      moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, molecule_atoms[startingBead]);

  if (molecule_atoms.size() == 1)
  {
    // update atom index
    return ChainRetraceData(firstBeadData.energies, firstBeadData.RosenbluthWeight, 0.0);
  }

  ChainRetraceData chainData;
  switch (growType)
  {
    case Component::GrowType::Rigid:
      chainData = CBMC::retraceRigidMoleculeChainDeletion(random, component, hasExternalField, forceField,
                                                          simulationBox, interpolationGrids, framework,
                                                          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW,
                                                          cutOffMoleculeVDW, cutOffCoulomb, molecule_atoms);
      break;
    case Component::GrowType::Flexible:
      chainData = CBMC::retraceFlexibleMoleculeChainDeletion(
          random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtomData, moleculeAtomData, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
          molecule_atoms, {component.startingBead});
      break;
    default:
      std::unreachable();
  }

  // update atom index
  return ChainRetraceData(firstBeadData.energies + chainData.energies,
                          firstBeadData.RosenbluthWeight * chainData.RosenbluthWeight, 0.0);
}
