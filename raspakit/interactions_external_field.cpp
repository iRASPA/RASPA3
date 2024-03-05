module;

module interactions_external_field;

import <complex>;
import <span>;
import <numbers>;
import <cmath>;
import <vector>;
import <iostream>;
import <algorithm>;
import <type_traits>;

import int3;
import double3;
import double3x3;
import atom;
import simulationbox;
import energy_status;
import energy_status_inter;
import units;
import energy_factor;
import force_factor;
import running_energy;
import component;
import forcefield;

void Interactions::computeExternalFieldEnergy([[maybe_unused]] const ForceField &forceField, [[maybe_unused]] const SimulationBox &simulationBox,
                                              [[maybe_unused]] std::span<const Atom> moleculeAtoms, [[maybe_unused]] RunningEnergy &energyStatus) noexcept
{
}

void Interactions::computeExternalFieldTailEnergy([[maybe_unused]] const ForceField &forceField, [[maybe_unused]] const SimulationBox &simulationBox,
                                                  [[maybe_unused]] std::span<const Atom> moleculeAtoms, [[maybe_unused]] RunningEnergy &energyStatus) noexcept
{
}

[[nodiscard]] std::optional<RunningEnergy>
Interactions::computeExternalFieldEnergyDifference([[maybe_unused]] const ForceField &forceField, [[maybe_unused]] const SimulationBox &simulationBox,
                                                   [[maybe_unused]] std::span<const Atom> newatoms, [[maybe_unused]] std::span<const Atom> oldatoms) noexcept
{
  return RunningEnergy();
}
