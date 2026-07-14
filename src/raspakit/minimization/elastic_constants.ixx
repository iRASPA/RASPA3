module;

export module elastic_constants;

import std;

import system;

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

export std::string writeElasticConstants(const ElasticConstantsResult& result);
