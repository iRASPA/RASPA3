module;

module cbmc;

import <vector>;
import <tuple>;
import <optional>;
import <span>;
import <iostream>;
import <algorithm>;
import <numeric>;

import randomnumbers;
import component;
import atom;
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
import cbmc_flexible_insertion;
import component;


[[nodiscard]] std::optional<ChainData> 
CBMC::growMoleculeSwapInsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                                const ForceField &forceField, const SimulationBox &simulationBox, 
                                std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, 
                                Component::GrowType growType, double cutOff, double cutOffCoulomb, 
                                size_t selectedComponent, size_t selectedMolecule, double scaling, 
                                std::vector<Atom> atoms, size_t numberOfTrialDirections) noexcept
{
  switch(growType)
  {
    default:
    return CBMC::growRigidMoleculeSwapInsertion(random, hasExternalField, components, forceField, simulationBox, frameworkAtoms, 
                                                moleculeAtoms, beta, cutOff, cutOffCoulomb, selectedComponent, 
                                                selectedMolecule, scaling, atoms, numberOfTrialDirections);
  }
}

[[nodiscard]] std::optional<ChainData> 
CBMC::growMoleculeReinsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                              const ForceField &forceField, const SimulationBox &simulationBox, 
                              std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, 
                              double cutOff, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule, 
                              std::span<Atom> molecule, size_t numberOfTrialDirections) noexcept 
{
  return CBMC::growRigidMoleculeReinsertion(random, hasExternalField, components, forceField, simulationBox, frameworkAtoms, 
                                            moleculeAtoms, beta, cutOff, cutOffCoulomb, selectedComponent, 
                                            selectedMolecule, molecule, numberOfTrialDirections);
}

[[nodiscard]] ChainData 
CBMC::retraceMoleculeReinsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                                 const ForceField &forceField, const SimulationBox &simulationBox, 
                                 std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                 double beta, double cutOff, double cutOffCoulomb, size_t selectedComponent,
                                 size_t selectedMolecule, std::span<Atom> molecule, double storedR, 
                                 size_t numberOfTrialDirections) noexcept
{
  return CBMC::retraceRigidMoleculeReinsertion(random, hasExternalField, components, forceField, simulationBox, frameworkAtoms, 
                                               moleculeAtoms, beta, cutOff, cutOffCoulomb, selectedComponent, 
                                               selectedMolecule, molecule, storedR, numberOfTrialDirections);
}

[[nodiscard]] ChainData 
CBMC::retraceMoleculeSwapDeletion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                                  const ForceField &forceField, const SimulationBox &simulationBox, 
                                  std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                  double beta, double cutOff, double cutOffCoulomb, size_t selectedComponent,
                                  size_t selectedMolecule, std::span<Atom> molecule, double scaling, double storedR, 
                                  size_t numberOfTrialDirections) noexcept
{
  return CBMC::retraceRigidMoleculeSwapDeletion(random, hasExternalField, components, forceField, simulationBox, frameworkAtoms, 
                                                moleculeAtoms, beta, cutOff,cutOffCoulomb, selectedComponent, 
                                                selectedMolecule, molecule, scaling, storedR, numberOfTrialDirections);
}

