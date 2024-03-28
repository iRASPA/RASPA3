export module interactions_framework_molecule;

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
import framework;
import component;

export namespace Interactions
{
  void computeFrameworkMoleculeEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                      std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
                                      RunningEnergy &energyStatus) noexcept;

  void computeFrameworkMoleculeTailEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                          std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
                                          RunningEnergy &energyStatus) noexcept;

  [[nodiscard]] std::optional<RunningEnergy>
  computeFrameworkMoleculeEnergyDifference(const ForceField &forceField, const SimulationBox &simulationBox,
                                           std::span<const Atom> frameworkAtoms, std::span<const Atom> newatoms,
                                           std::span<const Atom> oldatoms) noexcept;

  [[nodiscard]] RunningEnergy
  computeFrameworkMoleculeTailEnergyDifference(const ForceField &forceField, const SimulationBox &simulationBox,
                                               std::span<const Atom> frameworkAtoms, std::span<const Atom> newatoms,
                                               std::span<const Atom> oldatoms) noexcept;

  ForceFactor computeFrameworkMoleculeGradient(const ForceField &forceField, const SimulationBox &simulationBox,
                                               std::span<Atom> frameworkAtoms,
                                               std::span<Atom> moleculeAtoms) noexcept;

  [[nodiscard]] std::pair<EnergyStatus, double3x3>
  computeFrameworkMoleculeEnergyStrainDerivative(const ForceField &forceField, 
                                                 const std::vector<Framework> &frameworkComponents,
                                                 const std::vector<Component> &components,
                                                 const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms,
                                                 std::span<Atom> moleculeAtoms) noexcept;
};
