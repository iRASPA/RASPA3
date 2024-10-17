module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <optional>
#include <span>
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
double computeTranslationalKineticEnergy(std::span<const Molecule> moleculePositions);
double computeRotationalKineticEnergy(std::span<const Molecule> moleculePositions,
                                      const std::vector<Component> components);
double3 computeCenterOfMass(std::span<const Molecule> moleculePositions);
double3 computeCenterOfMassVelocity(std::span<const Molecule> moleculePositions);
double3 computeLinearMomentum(std::span<const Molecule> moleculePositions);
}  // namespace Integrators