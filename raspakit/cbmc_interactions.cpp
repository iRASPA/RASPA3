module;

module cbmc_interactions;

import <iomanip>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <optional>;
import <cmath>;
import <tuple>;
import <type_traits>;

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
import cbmc_interactions_external_field;
import cbmc_interactions_framework_molecule;
import cbmc_interactions_intermolecular;


inline std::pair<EnergyStatus, double3x3> 
pair_acc(const std::pair<EnergyStatus, double3x3> &lhs, const std::pair<EnergyStatus, double3x3> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

[[nodiscard]] const std::vector<std::pair<Atom, RunningEnergy>> 
CBMC::computeExternalNonOverlappingEnergies(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox, 
                                            std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
                                            [[maybe_unused]] double cutOffVDW, [[maybe_unused]] double cutOffCoulomb,  
                                            [[maybe_unused]] std::vector<Atom>& trialPositions) noexcept
{
    std::vector<std::pair<Atom, RunningEnergy>> energies{};
    for (auto it = trialPositions.begin(); it != trialPositions.end(); ++it)
    {
        // skip trial-positions that have an overlap in inter-molecular energy
        std::optional<RunningEnergy> externalFieldEnergy = 
          CBMC::computeExternalFieldEnergy(hasExternalField, forceField, simulationBox,
                                           cutOffVDW, cutOffCoulomb, { it,1 });
        if(!externalFieldEnergy.has_value()) continue;

        // skip trial-positions that have an overlap in inter-molecular energy
        std::optional<RunningEnergy> interEnergy = 
          CBMC::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtoms, 
                                            cutOffVDW, cutOffCoulomb, { it,1 });
        if(!interEnergy.has_value()) continue;

        // skip trial-positions that have an overlap in framework-molecule energy
        std::optional<RunningEnergy> frameworkEnergy = 
          CBMC::computeFrameworkMoleculeEnergy(forceField, simulationBox, frameworkAtoms, 
                                               cutOffVDW, cutOffCoulomb, { it,1 });
        if(!frameworkEnergy.has_value()) continue;

        energies.push_back(std::make_pair(*it, 
              externalFieldEnergy.value() + interEnergy.value() + frameworkEnergy.value()));
    }
    return energies;
}

const std::vector<std::pair<std::vector<Atom>,RunningEnergy>> 
CBMC::computeExternalNonOverlappingEnergies(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox, 
                                            std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, 
                                            [[maybe_unused]] double cutOffVDW, [[maybe_unused]] double cutOffCoulomb, 
                                            std::vector<std::vector<Atom>>& trialPositionSets, 
                                            std::make_signed_t<std::size_t> skip) noexcept
{
    std::vector<std::pair<std::vector<Atom>,RunningEnergy>> energies{};

    for (std::vector<Atom> trialPositionSet : trialPositionSets)
    {
        std::optional<RunningEnergy> eternalFieldEnergy = 
          CBMC::computeExternalFieldEnergy(hasExternalField, forceField, simulationBox,
                                           cutOffVDW, cutOffCoulomb, trialPositionSet);
        if(!eternalFieldEnergy.has_value()) continue;

        std::optional<RunningEnergy> interEnergy = 
          CBMC::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtoms, cutOffVDW, cutOffCoulomb, 
                                            trialPositionSet, skip);
        if(!interEnergy.has_value()) continue;

        std::optional<RunningEnergy> frameworkEnergy = 
          CBMC::computeFrameworkMoleculeEnergy(forceField, simulationBox, frameworkAtoms, cutOffVDW, cutOffCoulomb, 
                                               trialPositionSet, skip);
        if(!frameworkEnergy.has_value()) continue;

        energies.push_back(std::make_pair(trialPositionSet, 
              eternalFieldEnergy.value() + interEnergy.value() + frameworkEnergy.value()));
    }
    return energies;
}

const std::optional<RunningEnergy> 
CBMC::computeExternalNonOverlappingEnergyDualCutOff(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox, 
                                                    std::span<const Atom> frameworkAtoms, 
                                                    std::span<const Atom> moleculeAtoms, 
                                                    double cutOffVDW, double cutOffCoulomb, 
                                                    std::vector<Atom>& trialPositionSet) noexcept
{
  std::pair<std::vector<Atom>,RunningEnergy> energies;

  std::optional<RunningEnergy> externalFieldEnergy = 
    CBMC::computeExternalFieldEnergy(hasExternalField, forceField, simulationBox, cutOffVDW, cutOffCoulomb, 
                                     trialPositionSet);
  if(!externalFieldEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> interEnergy = 
    CBMC::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtoms, cutOffVDW, cutOffCoulomb, 
                                      trialPositionSet, -1);
  if(!interEnergy.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkEnergy = 
    CBMC::computeFrameworkMoleculeEnergy(forceField, simulationBox, frameworkAtoms, cutOffVDW, cutOffCoulomb, 
                                         trialPositionSet, -1);
  if(!frameworkEnergy.has_value()) return std::nullopt;

  return externalFieldEnergy.value() + interEnergy.value() + frameworkEnergy.value();
}
