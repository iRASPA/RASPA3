module;

#ifdef USE_LEGACY_HEADERS
#include <span>
#include <optional>
#include <tuple>
#include <vector>
#endif

export module interactions_framework_molecule;

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
import framework;
import component;

export namespace Interactions
{
  RunningEnergy computeFrameworkMoleculeEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                      std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms) noexcept;

  RunningEnergy computeFrameworkMoleculeTailEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                                   std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms) noexcept;

  [[nodiscard]] std::optional<RunningEnergy>
  computeFrameworkMoleculeEnergyDifference(const ForceField &forceField, const SimulationBox &simulationBox,
                                           std::span<const Atom> frameworkAtoms, std::span<const Atom> newatoms,
                                           std::span<const Atom> oldatoms) noexcept;

  [[nodiscard]] RunningEnergy
  computeFrameworkMoleculeTailEnergyDifference(const ForceField &forceField, const SimulationBox &simulationBox,
                                               std::span<const Atom> frameworkAtoms, std::span<const Atom> newatoms,
                                               std::span<const Atom> oldatoms) noexcept;

  RunningEnergy computeFrameworkMoleculeGradient(const ForceField &forceField, const SimulationBox &simulationBox,
                                                 std::span<Atom> frameworkAtoms,
                                                 std::span<Atom> moleculeAtoms) noexcept;

  [[nodiscard]] std::pair<EnergyStatus, double3x3>
  computeFrameworkMoleculeEnergyStrainDerivative(const ForceField &forceField, 
                                                 const std::vector<Framework> &frameworkComponents,
                                                 const std::vector<Component> &components,
                                                 const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms,
                                                 std::span<Atom> moleculeAtoms) noexcept;
};
