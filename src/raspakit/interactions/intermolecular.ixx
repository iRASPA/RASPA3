module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
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
import gradient_factor;
import forcefield;
import component;

export namespace Interactions
{
/**
 * \brief Computes the inter-molecular energy between atoms.
 *
 * Calculates the van der Waals and Coulombic energy contributions between all pairs of atoms
 * in \p moleculeAtoms, excluding interactions within the same molecule.
 *
 * \param forceField The force field parameters used for the energy calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms for which to compute inter-molecular energies.
 * \return The total inter-molecular energy contributions.
 */
RunningEnergy computeInterMolecularEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                          std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the tail correction for inter-molecular van der Waals energy.
 *
 * Calculates the long-range correction to the van der Waals energy between all pairs of atoms
 * in \p moleculeAtoms, excluding interactions within the same molecule.
 *
 * \param forceField The force field parameters used for the energy calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms for which to compute the tail energy correction.
 * \return The tail energy contributions to the total energy.
 */
RunningEnergy computeInterMolecularTailEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                              std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the difference in inter-molecular energy due to atom changes.
 *
 * Calculates the energy difference resulting from replacing \p oldatoms with \p newatoms
 * in the system, considering interactions with \p moleculeAtoms. Excludes interactions within
 * the same molecule. Returns std::nullopt if an overlap is detected based on the force field's overlap criteria.
 *
 * \param forceField The force field parameters used for the energy calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of existing atoms in the system.
 * \param newatoms A span of new atoms to be added to the system.
 * \param oldatoms A span of atoms to be removed from the system.
 * \return The energy difference due to the atom changes, or std::nullopt if an overlap occurs.
 */
[[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergyDifference(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Computes the difference in inter-molecular tail energy due to atom changes.
 *
 * Calculates the change in tail correction to the van der Waals energy resulting from replacing
 * \p oldatoms with \p newatoms in the system, considering interactions with \p moleculeAtoms.
 * Excludes interactions within the same molecule.
 *
 * \param forceField The force field parameters used for the energy calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of existing atoms in the system.
 * \param newatoms A span of new atoms to be added to the system.
 * \param oldatoms A span of atoms to be removed from the system.
 * \return The change in tail energy contributions due to the atom changes.
 */
[[nodiscard]] RunningEnergy computeInterMolecularTailEnergyDifference(const ForceField &forceField,
                                                                      const SimulationBox &simulationBox,
                                                                      std::span<const Atom> moleculeAtoms,
                                                                      std::span<const Atom> newatoms,
                                                                      std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Computes the inter-molecular forces and energy.
 *
 * Calculates the van der Waals and Coulombic forces between all pairs of atoms
 * in \p moleculeAtoms, excluding interactions within the same molecule. Updates the
 * gradient (force) field of each Atom accordingly.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms for which to compute inter-molecular forces and energies.
 * \return The total inter-molecular energy contributions.
 */
RunningEnergy computeInterMolecularGradient(const ForceField &forceField, const SimulationBox &simulationBox,
                                            std::span<Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes inter-molecular energy, forces, and strain derivative tensor.
 *
 * Calculates the van der Waals and Coulombic energies and forces between all pairs of atoms
 * in \p moleculeAtoms, excluding interactions within the same molecule. Updates the gradient
 * (force) field of each Atom and computes the strain derivative tensor for stress calculations.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param components The list of components in the system.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms for which to compute energies, forces, and strain derivatives.
 * \return A pair containing the energy status and the strain derivative tensor.
 */
std::pair<EnergyStatus, double3x3> computeInterMolecularEnergyStrainDerivative(const ForceField &forceField,
                                                                               const std::vector<Component> &components,
                                                                               const SimulationBox &simulationBox,
                                                                               std::span<Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the inter-molecular electric potential for each atom.
 *
 * Calculates the electric potential at each atom in \p moleculeAtoms due to the charges of
 * other atoms, excluding interactions within the same molecule. The results are stored in
 * \p electricPotentialMolecules.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param box The simulation box containing the atoms.
 * \param electricPotentialMolecules A span to store the computed electric potentials.
 * \param moleculeAtoms A span of atoms for which to compute the electric potentials.
 */
void computeInterMolecularElectricPotential(const ForceField &forceField, const SimulationBox &box,
                                            std::span<double> electricPotentialMolecules,
                                            std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the inter-molecular polarization energy.
 *
 * Calculates the energy contribution from polarization interactions between atoms in \p moleculeAtoms,
 * excluding interactions within the same molecule.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms for which to compute the polarization energy.
 * \return The total inter-molecular polarization energy.
 */
Potentials::EnergyFactor computeInterMolecularPolarizationEnergy(const ForceField &forceField,
                                                                 const SimulationBox &simulationBox,
                                                                 std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the inter-molecular electric field for each atom.
 *
 * Calculates the electric field at each atom in \p moleculeAtoms due to the charges of
 * other atoms, excluding interactions within the same molecule. The results are stored in
 * \p electricFieldMolecules. Also computes inter-molecular energy contributions.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param box The simulation box containing the atoms.
 * \param electricFieldMolecules A span to store the computed electric fields.
 * \param moleculeAtoms A span of atoms for which to compute the electric fields.
 * \return The total inter-molecular energy contributions.
 */
RunningEnergy computeInterMolecularElectricField(const ForceField &forceField, const SimulationBox &box,
                                                 std::span<double3> electricFieldMolecules,
                                                 std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the difference in inter-molecular electric field due to atom changes.
 *
 * Calculates the change in electric field at each atom in \p moleculeAtoms resulting from replacing
 * \p oldatoms with \p newatoms in the system, considering interactions with \p moleculeAtoms.
 * Excludes interactions within the same molecule. Updates \p electricFieldMolecules and \p electricFieldMolecule
 * accordingly. Returns std::nullopt if an overlap is detected based on the force field's overlap criteria.
 *
 * \param forceField The force field parameters used for the calculations.
 * \param box The simulation box containing the atoms.
 * \param electricFieldMolecules A span to store the updated electric fields of existing atoms.
 * \param electricFieldMolecule A span to store the updated electric fields of new atoms.
 * \param moleculeAtoms A span of existing atoms in the system.
 * \param newatoms A span of new atoms to be added to the system.
 * \param oldatoms A span of atoms to be removed from the system.
 * \return The energy difference due to the atom changes, or std::nullopt if an overlap occurs.
 */
std::optional<RunningEnergy> computeInterMolecularElectricFieldDifference(
    const ForceField &forceField, const SimulationBox &box, std::span<double3> electricFieldMolecules,
    std::span<double3> electricFieldMolecule, std::span<const Atom> moleculeAtoms, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept;
};  // namespace Interactions
