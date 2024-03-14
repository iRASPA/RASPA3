export module cbmc_rigid_reinsertion;

import <vector>;
import <optional>;
import <span>;

import atom;
import double3x3;
import double3;
import randomnumbers;
import forcefield;
import simulationbox;
import cbmc_chain_data;
import cbmc_interactions;
import component;


export namespace CBMC                                                                                                   
{   
  [[nodiscard]] std::optional<ChainData> 
  growRigidMoleculeReinsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                               const ForceField &forceField, const SimulationBox &simulationBox, 
                               std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                               double beta, double cutOff, double cutOffCoulomb, size_t selectedComponent, 
                               size_t selectedMolecule, std::span<Atom> molecule, size_t numberOfTrialDirections) 
                               noexcept;

  [[nodiscard]] ChainData 
  retraceRigidMoleculeReinsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                                  const ForceField &forceField, const SimulationBox &simulationBox, 
                                  std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                  double beta, double cutOff, double cutOffCoulomb, size_t selectedComponent, 
                                  size_t selectedMolecule, std::span<Atom> molecule, double storedR, 
                                  size_t numberOfTrialDirections);
}
