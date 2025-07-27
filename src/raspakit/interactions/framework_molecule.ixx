module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module interactions_framework_molecule;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import double3x3;
import double3x3x3;
import atom;
import running_energy;
import energy_status;
import simulationbox;
import energy_factor;
import gradient_factor;
import hessian_factor;
import forcefield;
import framework;
import component;
import interpolation_energy_grid;

export namespace Interactions
{
/**
 * \brief Computes the interaction energy between the framework and molecule atoms.
 *
 * Calculates the van der Waals and Coulombic interaction energy between atoms in the framework
 * and the molecule. It sums the contributions from all atom pairs within the specified cut-off distances.
 *
 * \param forceField The force field parameters for the simulation.
 * \param simulationBox The simulation box containing periodic boundary conditions.
 * \param frameworkAtoms A span of atoms representing the framework.
 * \param moleculeAtoms A span of atoms representing the molecule.
 * \return A RunningEnergy object containing the total interaction energy.
 */
RunningEnergy computeFrameworkMoleculeEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the tail correction energy between the framework and molecule atoms.
 *
 * Calculates the van der Waals tail correction energy for interactions between atoms in the framework
 * and the molecule. This correction accounts for the long-range interactions beyond the cut-off distance.
 *
 * \param forceField The force field parameters for the simulation.
 * \param simulationBox The simulation box containing periodic boundary conditions.
 * \param frameworkAtoms A span of atoms representing the framework.
 * \param moleculeAtoms A span of atoms representing the molecule.
 * \return A RunningEnergy object containing the tail correction energy.
 */
RunningEnergy computeFrameworkMoleculeTailEnergy(const ForceField &forceField, const SimulationBox &simulationBox,
                                                 std::span<const Atom> frameworkAtoms,
                                                 std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the difference in interaction energy between the framework and molecule atoms.
 *
 * Calculates the difference in van der Waals and Coulombic interaction energy between the framework
 * and the molecule atoms, due to changes from oldatoms to newatoms. If an overlap is detected (energy
 * exceeds overlap criteria), returns std::nullopt.
 *
 * \param forceField The force field parameters for the simulation.
 * \param simulationBox The simulation box containing periodic boundary conditions.
 * \param frameworkAtoms A span of atoms representing the framework.
 * \param newatoms A span of new atom positions representing the molecule.
 * \param oldatoms A span of old atom positions representing the molecule.
 * \return An optional RunningEnergy object containing the energy difference, or std::nullopt if overlap occurs.
 */

[[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergyDifference(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> framework, std::span<const Atom> frameworkAtoms, std::span<const Atom> newatoms,
    std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Computes the difference in interaction energy and electric-field between the framework and molecule atoms.
 *
 * Calculates the difference in van der Waals and Coulombic interaction energy between the framework
 * and the molecule atoms, due to changes from oldatoms to newatoms. If an overlap is detected (energy
 * exceeds overlap criteria), returns std::nullopt.
 *
 * \param forceField The force field parameters for the simulation.
 * \param simulationBox The simulation box containing periodic boundary conditions.
 * \param frameworkAtoms A span of atoms representing the framework.
 * \param newatoms A span of new atom positions representing the molecule.
 * \param oldatoms A span of old atom positions representing the molecule.
 * \param electricFieldMoleculeNew A span of new electric-field at the molecule atom positions
 * \param electricFieldMoleculeOld A span of old electric-field at the molecule atom positions
 * \return An optional RunningEnergy object containing the energy difference, or std::nullopt if overlap occurs.
 */

[[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergyDifference(
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> framework, std::span<const Atom> frameworkAtoms,
    std::span<double3> electricFieldMoleculeNew, std::span<double3> electricFieldMoleculeOld,
    std::span<const Atom> newatoms, std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Computes the difference in electric field due to changes in molecule atoms.
 *
 * Calculates the difference in Coulombic electric field at molecule atom positions due to changes from
 * oldatoms to newatoms, resulting from interactions with the framework atoms. Updates the electric field vector
 * accordingly. Also computes the difference in van der Waals and Coulombic interaction energy.
 * Returns std::nullopt if an overlap occurs.
 *
 * \param forceField The force field parameters for the simulation.
 * \param simulationBox The simulation box containing periodic boundary conditions.
 * \param frameworkAtoms A span of atoms representing the framework.
 * \param electricFieldMoleculeNew A span of new electric field
 * \param electricFieldMoleculeOld A span of old electric field
 * \param newatoms A span of new atom positions representing the molecule.
 * \param oldatoms A span of old atom positions representing the molecule.
 * \return An optional RunningEnergy object containing the energy difference, or std::nullopt if overlap occurs.
 */

void computeFrameworkMoleculeElectricFieldDifference(const ForceField &forceField, const SimulationBox &simulationBox,
                                                     std::span<const Atom> frameworkAtoms,
                                                     std::span<double3> electricFieldMoleculeNew,
                                                     std::span<double3> electricFieldMoleculeOld,
                                                     std::span<const Atom> newatoms,
                                                     std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Computes the difference in tail correction energy between the framework and molecule atoms.
 *
 * Calculates the difference in van der Waals tail correction energy due to changes in the molecule atoms
 * from oldatoms to newatoms. This accounts for the long-range interactions beyond the cut-off distance.
 *
 * \param forceField The force field parameters for the simulation.
 * \param simulationBox The simulation box containing periodic boundary conditions.
 * \param frameworkAtoms A span of atoms representing the framework.
 * \param newatoms A span of new atom positions representing the molecule.
 * \param oldatoms A span of old atom positions representing the molecule.
 * \return A RunningEnergy object containing the tail energy difference.
 */
[[nodiscard]] RunningEnergy computeFrameworkMoleculeTailEnergyDifference(const ForceField &forceField,
                                                                         const SimulationBox &simulationBox,
                                                                         std::span<const Atom> frameworkAtoms,
                                                                         std::span<const Atom> newatoms,
                                                                         std::span<const Atom> oldatoms) noexcept;

/**
 * \brief Computes the interaction energy and gradients between the framework and molecule atoms.
 *
 * Calculates the van der Waals and Coulombic interaction energy and updates the gradient vectors
 * (forces) for atoms in the framework and the molecule. The gradients are accumulated in the gradient
 * field of the Atom structures.
 *
 * \param forceField The force field parameters for the simulation.
 * \param simulationBox The simulation box containing periodic boundary conditions.
 * \param frameworkAtoms A span of atoms representing the framework; their gradients will be updated.
 * \param moleculeAtoms A span of atoms representing the molecule; their gradients will be updated.
 * \return A RunningEnergy object containing the total interaction energy.
 */
RunningEnergy computeFrameworkMoleculeGradient(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms,
    std::span<Atom> moleculeAtoms,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids) noexcept;

/**
 * \brief Computes the interaction energy, gradients, and strain derivative between the framework and molecule atoms.
 *
 * Calculates the van der Waals and Coulombic interaction energy, updates the gradient vectors (forces) for atoms
 * in the framework and the molecule, and computes the strain derivative tensor for the stress calculation.
 *
 * \param forceField The force field parameters for the simulation.
 * \param frameworkComponents A vector of frameworks in the simulation.
 * \param components A vector of components (molecules) in the simulation.
 * \param simulationBox The simulation box containing periodic boundary conditions.
 * \param frameworkAtoms A span of atoms representing the framework; their gradients will be updated.
 * \param moleculeAtoms A span of atoms representing the molecule; their gradients will be updated.
 * \return A pair consisting of EnergyStatus containing the interaction energies, and a double3x3 tensor representing
 * the strain derivative.
 */
[[nodiscard]] std::pair<EnergyStatus, double3x3> computeFrameworkMoleculeEnergyStrainDerivative(
    const ForceField &forceField, const std::optional<Framework> &framework,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::vector<Component> &components, const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms,
    std::span<Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the electric potential at molecule atom positions due to the framework atoms.
 *
 * Calculates the Coulombic electric potential at each molecule atom position, resulting from interactions
 * with the framework atoms, and accumulates the values in the provided electric potential vector.
 *
 * \param forceField The force field parameters for the simulation.
 * \param simulationBox The simulation box containing periodic boundary conditions.
 * \param electricPotentialMolecules A span of doubles where the electric potentials will be accumulated.
 * \param frameworkAtoms A span of atoms representing the framework.
 * \param moleculeAtoms A span of atoms representing the molecule.
 */
void computeFrameworkMoleculeElectrostaticPotential(const ForceField &forceField, const SimulationBox &simulationBox,
                                                    std::span<double> electricPotentialMolecules,
                                                    std::span<const Atom> frameworkAtoms,
                                                    std::span<const Atom> moleculeAtoms) noexcept;

/**
 * \brief Computes the electric field at molecule atom positions due to the framework atoms.
 *
 * Calculates the Coulombic electric field at each molecule atom position, resulting from interactions
 * with the framework atoms, and accumulates the values in the provided electric field vector.
 * Also computes the van der Waals and Coulombic interaction energy between the framework and molecule atoms.
 *
 * \param forceField The force field parameters for the simulation.
 * \param simulationBox The simulation box containing periodic boundary conditions.
 * \param electricField A span of double3 where the electric field vectors will be accumulated.
 * \param frameworkAtoms A span of atoms representing the framework.
 * \param moleculeAtoms A span of atoms representing the molecule.
 * \return A RunningEnergy object containing the total interaction energy.
 */
RunningEnergy computeFrameworkMoleculeElectricField(const ForceField &forceField, const SimulationBox &simulationBox,
                                                    std::span<double3> electricField,
                                                    std::span<const Atom> frameworkAtoms,
                                                    std::span<const Atom> moleculeAtoms) noexcept;

std::tuple<double, double3, double3x3> calculateHessianAtPositionVDW(const ForceField &forceField,
                                                                     const SimulationBox &simulationBox, double3 posA,
                                                                     std::size_t typeA,
                                                                     std::span<const Atom> frameworkAtoms);

std::tuple<double, double3, double3x3> calculateHessianAtPositionCoulomb(const ForceField &forceField,
                                                                         const SimulationBox &simulationBox,
                                                                         double3 posA, double chargeA,
                                                                         std::span<const Atom> frameworkAtoms);

};  // namespace Interactions
