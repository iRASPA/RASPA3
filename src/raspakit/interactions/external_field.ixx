module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module interactions_external_field;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import double3x3;
import atom;
import running_energy;
import energy_status;
import simulationbox;
import gradient_factor;
import forcefield;
import component;

export namespace Interactions
{
/**
 * \brief Computes the external field energy contribution for the system.
 *
 * Calculates the energy contribution from an external field acting on a set of atoms,
 * updating the running energy status accordingly.
 *
 * \param hasExternalField Flag indicating if the external field is active.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms to compute the energy for.
 * \param energyStatus The running total of energies to be updated.
 */
void computeExternalFieldEnergy(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
                                std::span<const Atom> moleculeAtoms, RunningEnergy &energyStatus) noexcept;

/**
 * \brief Computes the tail correction for the external field energy.
 *
 * Calculates the tail corrections for the external field energy, accounting for
 * interactions beyond the cutoff distance.
 *
 * \param hasExternalField Flag indicating if the external field is active.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms to compute the tail energy for.
 * \param energyStatus The running total of energies to be updated.
 */
void computeExternalFieldTailEnergy(bool hasExternalField, const ForceField &forceField,
                                    const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
                                    RunningEnergy &energyStatus) noexcept;

/**
 * \brief Computes the difference in external field energy between two states.
 *
 * Calculates the energy difference due to an external field when atoms have moved or changed,
 * which is useful for evaluating Monte Carlo moves. Returns an optional RunningEnergy,
 * which is `std::nullopt` if the energy exceeds the overlap criteria.
 *
 * \param hasExternalField Flag indicating if the external field is active.
 * \param forceField The force field parameters.
 * \param simulationBox The simulation box containing the atoms.
 * \param newatoms A span of new atom positions and properties.
 * \param oldatoms A span of old atom positions and properties.
 * \return An optional RunningEnergy containing the energy difference, or `std::nullopt` if overlap occurs.
 */
[[nodiscard]] std::optional<RunningEnergy> computeExternalFieldEnergyDifference(
    bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept;

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>>
calculateThirdDerivativeAtPositionExternalField(const ForceField &forceField, const SimulationBox &simulationBox,
                                                double3 posA);

std::array<double, 8> calculateTricubicFractionalAtPositionExternalField(const ForceField &forceField,
                                                                         const SimulationBox &simulationBox,
                                                                         double3 posA);

std::tuple<double, std::array<double, 3>, std::array<std::array<double, 3>, 3>,
           std::array<std::array<std::array<double, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>,
           std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>>
calculateSixthOrderDerivativeAtPositionExternalField(const ForceField &forceField, double3 pos);

std::array<double, 27> calculateTriquinticFractionalAtPositionExternalField(const ForceField &forceField,
                                                                            const SimulationBox &simulationBox,
                                                                            double3 posA);

}  // namespace Interactions
