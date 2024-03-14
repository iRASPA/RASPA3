export module cbmc;

import <vector>;
import <optional>;
import <span>;

import atom;
import double3x3;
import double3;
import randomnumbers;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import cbmc_chain_data;
import component;
import forcefield;
import simulationbox;


export namespace CBMC
{
  // insertion
  [[nodiscard]] std::optional<ChainData>
  growMoleculeSwapInsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                            const ForceField &forceField, const SimulationBox &simulationBox, 
                            std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, 
                            Component::GrowType growType, double cutOff, double cutOffCoulomb, 
                            size_t selectedComponent, size_t selectedMolecule, double scaling, 
                            std::vector<Atom> atoms, size_t numberOfTrialDirections) noexcept;
  
  // deletion
  [[nodiscard]] ChainData
  retraceMoleculeSwapDeletion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                              const ForceField &forceField, const SimulationBox &simulationBox, 
                              std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, 
                              double cutOff, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule,
                              std::span<Atom> molecule, double scaling, double storedR, 
                              size_t numberOfTrialDirections) noexcept; 
  
  // reinsertion grow
  [[nodiscard]] std::optional<ChainData>
  growMoleculeReinsertion(RandomNumber &random, bool hasExternalField, 
                          const std::vector<Component> &components, const ForceField &forceField, 
                          const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms, 
                          std::span<const Atom> moleculeAtoms, double beta, double cutOff, double cutOffCoulomb, 
                          size_t selectedComponent, size_t selectedMolecule, std::span<Atom> molecule, 
                          size_t numberOfTrialDirections) noexcept;
  
  // reinsertion retrace
  [[nodiscard]] ChainData
  retraceMoleculeReinsertion(RandomNumber &random, bool hasExternalField, const std::vector<Component> &components, 
                             const ForceField &forceField, const SimulationBox &simulationBox, 
                             std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, double beta, 
                             double cutOff, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule, 
                             std::span<Atom> molecule, double storedR, size_t numberOfTrialDirections) noexcept;
}


