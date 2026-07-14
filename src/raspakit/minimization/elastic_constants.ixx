module;

export module elastic_constants;

import std;

import system;
import double3x3;

export struct ElasticConstantsResult
{
  // Voigt order: xx, yy, zz, yz, xz, xy. Values use internal pressure units.
  std::array<double, 36> born{};
  std::array<double, 36> relaxation{};
  std::array<double, 36> pressureCorrection{};
  std::array<double, 36> stiffness{};
  std::array<double, 36> compliance{};
  std::array<double, 6> stabilityEigenvalues{};
  std::array<double, 3> youngModuli{};
  // poissonRatios[loadingDirection * 3 + transverseDirection].
  std::array<double, 9> poissonRatios{};
  double bulkModulusVoigt{};
  double shearModulusVoigt{};
  double bulkModulusReuss{};
  double shearModulusReuss{};
  double bulkModulusHill{};
  double shearModulusHill{};
  std::size_t discardedInternalModes{};
  bool complianceAvailable{};
};

/**
 * Compute the relaxed, static elastic tensor at the current (minimized) structure.
 *
 * The Born and relaxation terms are evaluated from the analytic generalized
 * Hessian. The cell Hessian is converted from logarithmic to infinitesimal
 * symmetric strain before applying the hydrostatic-pressure correction.
 */
export ElasticConstantsResult computeElasticConstants(const System& system,
                                                      double relativeEigenvalueTolerance = 1.0e-8);

/**
 * Compute the instantaneous affine (Born) tensor at the current configuration.
 *
 * No internal-coordinate relaxation or external-pressure correction is applied.
 */
export std::array<double, 36> computeAffineBornTensor(const System& system);

/**
 * Return the instantaneous molecular kinetic virial, sum(m v outer v), before volume normalization.
 */
export double3x3 computeMolecularKineticVirial(const System& system);

export std::string writeElasticConstants(const ElasticConstantsResult& result);
