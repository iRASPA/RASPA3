module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <print>
#include <tuple>
#endif

export module rigid;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <tuple>;
import <print>;
#endif

import archive;
import double3;
import simd_quatd;
import stringutils;

export namespace Rigid
{
/**
 * \brief Integrates the rotational motion of a rigid body using a second-order splitting method.
 *
 * This function updates the orientation and angular momentum quaternions of a rigid body over a time step \p dt,
 * using a second-order integrator based on operator splitting. It applies rotations about the principal axes
 * in a specific sequence to achieve second-order accuracy in time.
 *
 * \param dt The total time step over which to integrate.
 * \param q A pair of quaternions \p (p, q) representing the current orientation (\p p) and angular momentum (\p q) of
 * the rigid body. \param inverseInertiaVector The inverse of the rigid body's principal moments of inertia, as a \p
 * double3 vector. \return A pair of quaternions representing the updated orientation and angular momentum after
 * integration.
 */
std::pair<simd_quatd, simd_quatd> NoSquishFreeRotorOrderTwo(double dt, std::pair<simd_quatd, simd_quatd> q,
                                                            double3 inverseInertiaVector);
}  // namespace Rigid

namespace Rigid
{
/**
 * \brief Applies a rotation about a specified principal axis to a rigid body's orientation and angular momentum
 * quaternions.
 *
 * This function performs a rotation of the rigid body's orientation (\p p) and angular momentum (\p q) quaternions
 * about one of the principal axes, specified by \p k (1 for x-axis, 2 for y-axis, 3 for z-axis). The rotation is
 * computed using the angular velocity derived from the rigid body's angular momentum and inverse inertia tensor.
 *
 * \param k The index of the principal axis (1: x-axis, 2: y-axis, 3: z-axis) about which to rotate.
 * \param dt The time increment over which to perform the rotation.
 * \param q A pair of quaternions \p (p, q) representing the current orientation (\p p) and angular momentum (\p q) of
 * the rigid body. \param inverseInertiaVector The inverse of the rigid body's principal moments of inertia, as a \p
 * double3 vector. \return A pair of quaternions representing the updated orientation and angular momentum after the
 * rotation.
 */
std::pair<simd_quatd, simd_quatd> NoSquishRotate(size_t k, double dt, std::pair<simd_quatd, simd_quatd> q,
                                                 double3 inverseInertiaVector);
}  // namespace Rigid
