export module interactions_intermolecular;

import <span>;
import <optional>;
import <tuple>;

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

  ForceFactor computeInterMolecularGradient(const ForceField &forceField, const SimulationBox &simulationBox, 
                                            std::span<Atom> moleculeAtoms) noexcept;

  std::pair<EnergyStatus, double3x3>                                                                                      
  computeInterMolecularEnergyStrainDerivative(const ForceField &forceField,
                                              const std::vector<Component> &components,
                                              const SimulationBox &simulationBox,
                                              std::span<Atom> moleculeAtoms) noexcept;
};
