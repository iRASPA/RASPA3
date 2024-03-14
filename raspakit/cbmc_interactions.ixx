export module cbmc_interactions;

import <vector>;
import <tuple>;
import <type_traits>;
import <span>;
import <optional>;

import atom;
import energy_factor;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import component;
import double3;
import double3x3;
import forcefield;
import simulationbox;
import units;
import cbmc_interactions_intermolecular;
import cbmc_interactions_framework_molecule;


export namespace CBMC                                                                                                   
{ 
  [[nodiscard]] const std::vector<std::pair<Atom, RunningEnergy>> 
  computeExternalNonOverlappingEnergies(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox, 
                                        std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                        double cutOffVDW, double cutOffCoulomb,  std::vector<Atom>& trialPositions) 
                                        noexcept;
  
  const std::vector<std::pair<std::vector<Atom>,RunningEnergy>> 
  computeExternalNonOverlappingEnergies(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox, 
                                        std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                        double cutOffVDW, double cutOffCoulomb, 
                                        std::vector<std::vector<Atom>>& trialPositionSets, 
                                        std::make_signed_t<std::size_t> skip = -1) noexcept;
  
  const std::optional<RunningEnergy> 
  computeExternalNonOverlappingEnergyDualCutOff(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox, 
                                                std::span<const Atom> frameworkAtoms, 
                                                std::span<const Atom> moleculeAtoms, 
                                                double cutOffVDW, double cutOffCoulomb, 
                                                std::vector<Atom>& trialPositionSet) noexcept;
}
