module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>
#endif

export module integrators_compute;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <optional>;
#endif

import molecule;
import atom;
import double3;
import component;

export namespace Integrators
{
/**
 * \brief Computes the total translational kinetic energy of molecules.
 *
 * Calculates the sum of translational kinetic energies for all molecules in the system.
 *
 * \param moleculePositions A span of molecules for which to compute the kinetic energy.
 * \return The total translational kinetic energy.
 */
double computeTranslationalKineticEnergy(std::span<const Molecule> moleculePositions);

/**
 * \brief Computes the total rotational kinetic energy of molecules.
 *
 * Calculates the sum of rotational kinetic energies for all molecules in the system,
 * based on their angular velocities and moments of inertia.
 *
 * \param moleculePositions A span of molecules for which to compute the rotational kinetic energy.
 * \param components A vector of components containing inertia information.
 * \return The total rotational kinetic energy.
 */
double computeRotationalKineticEnergy(std::span<const Molecule> moleculePositions,
                                      const std::vector<Component> components);

/**
 * \brief Computes the center of mass position of the system.
 *
 * Calculates the weighted average position of all molecules in the system.
 *
 * \param moleculePositions A span of molecules for which to compute the center of mass.
 * \return The center of mass position as a double3 vector.
 */
double3 computeCenterOfMass(std::span<const Molecule> moleculePositions);

/**
 * \brief Computes the velocity of the center of mass of the system.
 *
 * The velocity of the center of mass is the average velocity of all objects in the system weighted by their masses.
 *
 * \param moleculePositions A span of molecules for which to compute the center of mass velocity.
 * \return The center of mass velocity as a double3 vector.
 */
double3 computeCenterOfMassVelocity(std::span<const Molecule> moleculePositions);

/**
 * \brief Computes the total linear momentum of the system.
 *
 * Calculates the sum of linear momenta for all molecules in the system.
 *
 * \param moleculePositions A span of molecules for which to compute the linear momentum.
 * \return The total linear momentum as a double3 vector.
 */
double3 computeLinearMomentum(std::span<const Molecule> moleculePositions);
}  // namespace Integrators
