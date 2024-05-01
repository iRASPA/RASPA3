module;

#ifdef USE_LEGACY_HEADERS
#include <span>
#include <optional>
#include <tuple>
#include <vector>
#endif

export module interactions_intermolecular;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <optional>;
import <tuple>;
import <vector>;
#endif

import double3x3;
import atom;
import running_energy;
import energy_status;
import simulationbox;
import force_factor;
import forcefield;
import component;

export namespace Interactions
{
  void computeInterMolecularEnergy(const ForceField &forceField, const SimulationBox &simulationBox, 
                                   std::span<const Atom> moleculeAtoms, RunningEnergy &energyStatus) noexcept;

  void computeInterMolecularTailEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                       std::span<const Atom> moleculeAtoms, 
                                       RunningEnergy &energyStatus) noexcept;

  [[nodiscard]] std::optional<RunningEnergy> 
  computeInterMolecularEnergyDifference(const ForceField &forceField, const SimulationBox &simulationBox, 
                                        std::span<const Atom> moleculeAtoms, 
                                        std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept;

  [[nodiscard]] RunningEnergy 
  computeInterMolecularTailEnergyDifference(const ForceField &forceField, const SimulationBox &simulationBox, 
                                            std::span<const Atom> moleculeAtoms, 
                                            std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept;

  std::pair<ForceFactor, ForceFactor> computeInterMolecularGradient(const ForceField &forceField, const SimulationBox &simulationBox, 
                                                                    std::span<Atom> moleculeAtoms) noexcept;

  std::pair<EnergyStatus, double3x3>                                                                                      
  computeInterMolecularEnergyStrainDerivative(const ForceField &forceField,
                                              const std::vector<Component> &components,
                                              const SimulationBox &simulationBox,
                                              std::span<Atom> moleculeAtoms) noexcept;
};
