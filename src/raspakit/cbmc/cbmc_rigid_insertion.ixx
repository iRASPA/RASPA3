module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <optional>
#include <span>
#endif

export module cbmc_rigid_insertion;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <optional>;
import <span>;
#endif

import atom;
import double3x3;
import double3;
import randomnumbers;
import forcefield;
import simulationbox;
import cbmc_chain_data;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import component;


export namespace CBMC                                                                                                   
{
  [[nodiscard]] std::optional<ChainData> 
  growRigidMoleculeSwapInsertion(RandomNumber &random, bool hasExternalField,  const std::vector<Component> &components, 
                                 const ForceField &forceField, const SimulationBox &simulationBox, 
                                 std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                 double beta, double cutOff, double cutOffCoulomb, size_t selectedComponent, 
                                 size_t selectedMolecule, double scaling, std::vector<Atom> atoms, 
                                 size_t numberOfTrialDirections) noexcept;

  [[nodiscard]] std::optional<ChainData>                                                                                  
  growRigidMoleculeChain(RandomNumber &random, bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,    
                         std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta,    
                         double cutOff, double cutOffCoulomb, size_t startingBead,                                  
                         std::vector<Atom> molecule, size_t numberOfTrialDirections) noexcept;
}
