module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <tuple>
#include <type_traits>
#include <span>
#include <optional>
#endif

export module cbmc_interactions;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <tuple>;
import <type_traits>;
import <span>;
import <optional>;
#endif

import atom;
import molecule;
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

  const std::vector<std::tuple<Molecule, std::vector<Atom>,RunningEnergy>> 
  computeExternalNonOverlappingEnergies(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox, 
                                        std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                        double cutOffVDW, double cutOffCoulomb, 
                                        std::vector<std::pair<Molecule, std::vector<Atom>>>& trialPositionSets, 
                                        std::make_signed_t<std::size_t> skip = -1) noexcept;
  
  const std::optional<RunningEnergy> 
  computeExternalNonOverlappingEnergyDualCutOff(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox, 
                                                std::span<const Atom> frameworkAtoms, 
                                                std::span<const Atom> moleculeAtoms, 
                                                double cutOffVDW, double cutOffCoulomb, 
                                                std::vector<Atom>& trialPositionSet) noexcept;
}
