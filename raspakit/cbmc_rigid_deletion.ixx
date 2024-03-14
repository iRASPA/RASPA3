export module cbmc_rigid_deletion;

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
  [[nodiscard]] ChainData
  retraceRigidMoleculeSwapDeletion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                                   const ForceField &forcefield, const SimulationBox &simulationBox, 
                                   std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                   double beta, double cutOff, double cutOffCoulomb, size_t selectedComponent, 
                                   size_t selectedMolecule, std::span<Atom> molecule, double scaling, double storedR, 
                                   size_t numberOfTrialDirections) noexcept;
}
