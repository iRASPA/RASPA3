module;

export module elastic_constants;

import std;

import system;
import double3x3;
import archive;

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

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ElasticConstantsResult& r);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ElasticConstantsResult& r);
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
 * Return the instantaneous kinetic virial of the barostat's coupling points, sum(m v outer v), before
 * volume normalization: rigid molecules and rigid groups contribute their center-of-mass momentum flux,
 * flexible atoms contribute individually. This is the kinetic partner of 'computeBarostatVirial'.
 */
export double3x3 computeMolecularKineticVirial(const System& system);

/**
 * Convert the molecular virial (sum over molecules of R_com outer F_molecule, the quantity behind the
 * Monte Carlo and the reported pressure) into the virial conjugate to the coupling the molecular-dynamics
 * barostat actually applies: the cell drives the centers of mass of rigid molecules and rigid groups and
 * every flexible atom individually ('propagateCell'). For a coupled point k of a molecule at R_k with
 * total force F_k (non-bonded and bonded) the virial is sum_k R_k outer F_k, so
 *
 *     V_barostat = V_molecular + sum_molecules sum_k (R_k - R_com) outer F_k .
 *
 * The extra term vanishes for a rigid molecule (one coupled point, its center of mass) and is what makes
 * the flexible atoms' kinetic term N_atoms k T the right partner of the virial; feeding the molecular
 * virial with the atomic kinetic energy over-estimates the pressure by (N_atoms - N_molecules) k T / V,
 * which for a liquid of flexible molecules is of the order of a thousand bar and blows the cell up.
 * Framework atoms are left as they enter the molecular virial (atomically); a flexible framework with
 * rigid groups is not corrected. Uses the current atom positions and the stored total atom gradients.
 */
export double3x3 computeBarostatVirial(const System& system, const double3x3& molecularVirial);

export std::string writeElasticConstants(const ElasticConstantsResult& result);
