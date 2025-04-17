module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module cbmc;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <tuple>;
import <optional>;
import <span>;
import <iostream>;
import <algorithm>;
import <numeric>;
#endif

import randomnumbers;
import component;
import atom;
import molecule;
import double3;
import double3x3;
import simulationbox;
import energy_status;
import cbmc_growing_status;
import forcefield;
import energy_factor;
import cbmc_rigid_insertion;
import cbmc_rigid_deletion;
import cbmc_rigid_reinsertion;
// import cbmc_flexible_insertion;
import framework;
import component;
import cbmc_chain_data;

[[nodiscard]] std::optional<ChainData> CBMC::growMoleculeSwapInsertion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const std::vector<Component> &components, const ForceField &forceField, const SimulationBox &simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta,
    Component::GrowType growType, double cutOffFrameworkVDW, double cutOffMoleculeVDW, double cutOffCoulomb,
    size_t selectedComponent, size_t selectedMolecule, double scaling, size_t groupId,
    size_t numberOfTrialDirections) noexcept
{
  switch (growType)
  {
    default:
      return CBMC::growRigidMoleculeSwapInsertion(
          random, framework, component, hasExternalField, components, forceField, simulationBox, frameworkAtoms,
          moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent,
          selectedMolecule, scaling, groupId, numberOfTrialDirections);
  }
}

[[nodiscard]] std::optional<ChainData> CBMC::growMoleculeReinsertion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const std::vector<Component> &components, const ForceField &forceField, const SimulationBox &simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule,
    Molecule &molecule, std::span<Atom> molecule_atoms, size_t numberOfTrialDirections) noexcept
{
  return CBMC::growRigidMoleculeReinsertion(random, framework, component, hasExternalField, components, forceField,
                                            simulationBox, frameworkAtoms, moleculeAtoms, beta, cutOffFrameworkVDW,
                                            cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, selectedMolecule,
                                            molecule, molecule_atoms, numberOfTrialDirections);
}

[[nodiscard]] ChainData CBMC::retraceMoleculeReinsertion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const std::vector<Component> &components, const ForceField &forceField, const SimulationBox &simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule,
    Molecule &molecule, std::span<Atom> molecule_atoms, double storedR, size_t numberOfTrialDirections) noexcept
{
  return CBMC::retraceRigidMoleculeReinsertion(random, framework, component, hasExternalField, components, forceField,
                                               simulationBox, frameworkAtoms, moleculeAtoms, beta, cutOffFrameworkVDW,
                                               cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, selectedMolecule,
                                               molecule, molecule_atoms, storedR, numberOfTrialDirections);
}

[[nodiscard]] ChainData CBMC::retraceMoleculeSwapDeletion(
    RandomNumber &random, const Framework &framework, const Component &component, bool hasExternalField,
    const std::vector<Component> &components, const ForceField &forceField, const SimulationBox &simulationBox,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule,
    std::span<Atom> molecule, double scaling, size_t numberOfTrialDirections) noexcept
{
  return CBMC::retraceRigidMoleculeSwapDeletion(random, framework, component, hasExternalField, components, forceField,
                                                simulationBox, frameworkAtoms, moleculeAtoms, beta, cutOffFrameworkVDW,
                                                cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, selectedMolecule,
                                                molecule, scaling, numberOfTrialDirections);
}
