module;

export module integrators_compute;

import std;

import molecule;
import atom;
import atom_dynamics;
import double3;
import component;
import framework;
import forcefield;

export namespace Integrators
{
/**
 * \brief Computes the total translational kinetic energy of molecules.
 *
 * Calculates the sum of translational kinetic energies for all movable atoms and molecules in the system.
 * For rigid molecules this uses the center-of-mass velocity; for flexible molecules the
 * atomic velocities carry all kinetic energy. Flexible-framework atomic kinetic energy is
 * included when framework state and force-field masses are supplied.
 *
 * \param moleculeData A span of molecules for which to compute the kinetic energy.
 * \param moleculeAtomPositions A span of atoms corresponding to the molecules.
 * \param components A vector of components containing rigidity and mass information.
 * \param framework Optional framework definition.
 * \param frameworkAtomPositions Framework atoms used to obtain pseudo-atom masses.
 * \param frameworkDynamics Per-framework-atom velocities.
 * \param forceField Force field containing framework pseudo-atom masses.
 * \return The total translational kinetic energy.
 */
double computeTranslationalKineticEnergy(std::span<const Molecule> moleculeData,
                                         std::span<const Atom> moleculeAtomPositions,
                                         std::span<const AtomDynamics> moleculeDynamics,
                                         const std::vector<Component>& components,
                                         const std::optional<Framework>& framework = std::nullopt,
                                         std::span<const Atom> frameworkAtomPositions = {},
                                         std::span<const AtomDynamics> frameworkDynamics = {},
                                         const ForceField* forceField = nullptr,
                                         std::span<const GroupState> groupData = {});

/**
 * \brief Computes the total rotational kinetic energy of molecules.
 *
 * Calculates the sum of rotational kinetic energies for all molecules in the system,
 * based on their angular velocities and moments of inertia.
 *
 * \param moleculeData A span of molecules for which to compute the rotational kinetic energy.
 * \param components A vector of components containing inertia information.
 * \return The total rotational kinetic energy.
 */
double computeRotationalKineticEnergy(std::span<const Molecule> moleculeData, const std::vector<Component> components,
                                      std::span<const GroupState> groupData = {});

/**
 * \brief Computes the center of mass position of the system.
 *
 * Calculates the weighted average position of all molecules in the system.
 *
 * \param moleculeData A span of molecules for which to compute the center of mass.
 * \return The center of mass position as a double3 vector.
 */
double3 computeCenterOfMass(std::span<const Molecule> moleculeData);

/**
 * \brief Computes the velocity of the center of mass of the system.
 *
 * The velocity of the center of mass is the average velocity of all objects in the system weighted by their masses.
 *
 * \param moleculeData A span of molecules for which to compute the center of mass velocity.
 * \return The center of mass velocity as a double3 vector.
 */
double3 computeCenterOfMassVelocity(std::span<const Molecule> moleculeData);

/**
 * \brief Computes the total linear momentum of the system.
 *
 * Calculates the sum of linear momenta for all molecules in the system.
 *
 * \param moleculeData A span of molecules for which to compute the linear momentum.
 * \return The total linear momentum as a double3 vector.
 */
double3 computeLinearMomentum(std::span<const Molecule> moleculeData);
}  // namespace Integrators
