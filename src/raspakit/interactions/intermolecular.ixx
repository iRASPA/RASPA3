module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
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

import double3;
import double3x3;
import atom;
import running_energy;
import energy_status;
import simulationbox;
import energy_factor;
import force_factor;
import forcefield;
import component;

export namespace Interactions
{
RunningEnergy computeInterMolecularEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                          std::span<const Atom> moleculeAtoms) noexcept;

RunningEnergy computeInterMolecularTailEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                              std::span<const Atom> moleculeAtoms) noexcept;

[[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergyDifference(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept;

[[nodiscard]] RunningEnergy computeInterMolecularTailEnergyDifference(const ForceField &forceField,
                                                                      const SimulationBox &simulationBox,
                                                                      std::span<const Atom> moleculeAtoms,
                                                                      std::span<const Atom> newatoms,
                                                                      std::span<const Atom> oldatoms) noexcept;

RunningEnergy computeInterMolecularGradient(const ForceField &forceField, const SimulationBox &simulationBox,
                                            std::span<Atom> moleculeAtoms) noexcept;

std::pair<EnergyStatus, double3x3> computeInterMolecularEnergyStrainDerivative(const ForceField &forceField,
                                                                               const std::vector<Component> &components,
                                                                               const SimulationBox &simulationBox,
                                                                               std::span<Atom> moleculeAtoms) noexcept;

void computeInterMolecularElectricPotential(const ForceField &forceField, const SimulationBox &box,
                                            std::span<double> electricPotentialMolecules,
                                            std::span<const Atom> moleculeAtoms) noexcept;

EnergyFactor computeInterMolecularPolarizationEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                                     std::span<const Atom> moleculeAtoms) noexcept;

RunningEnergy computeInterMolecularElectricField(const ForceField &forceField, const SimulationBox &box,
                                                 std::span<double3> electricFieldMolecules,
                                                 std::span<const Atom> moleculeAtoms) noexcept;

std::optional<RunningEnergy> computeInterMolecularElectricFieldDifference(
    const ForceField &forceField, const SimulationBox &box, std::span<double3> electricFieldMolecules,
    std::span<double3> electricFieldMolecule, std::span<const Atom> moleculeAtoms, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept;
};  // namespace Interactions
