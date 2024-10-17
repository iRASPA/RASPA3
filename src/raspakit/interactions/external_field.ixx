module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module interactions_external_field;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <optional>;
import <tuple>;
import <complex>;
import <vector>;
#endif

import double3;
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
void computeExternalFieldEnergy(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
                                std::span<const Atom> moleculeAtoms, RunningEnergy &energyStatus) noexcept;

void computeExternalFieldTailEnergy(bool hasExternalField, const ForceField &forceField,
                                    const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
                                    RunningEnergy &energyStatus) noexcept;

[[nodiscard]] std::optional<RunningEnergy> computeExternalFieldEnergyDifference(
    bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept;
}  // namespace Interactions
