module;

export module phonon_dynamical_matrix;

import std;

import double3;
import system;
import phonon_force_constants;

/**
 * k-dependent dynamical matrix and phonon frequencies built from the image-resolved force constants.
 *
 * The dynamical matrix is the mass-weighted lattice Fourier transform of the force constants:
 *
 *   D_{i alpha, j beta}(k) = (1/sqrt(m_i m_j)) sum_R Phi_{i alpha, j beta}(R) exp(2 pi i k . R),
 *
 * where k is expressed in fractional (crystallographic) reciprocal-lattice coordinates and R runs over
 * the integer lattice vectors. Diagonalizing the (Hermitian) matrix yields the squared angular
 * frequencies omega^2; the eigenvalues use the same internal units as the Gamma-point normal-mode
 * analysis, so at k = 0 they reproduce computeNormalModes for a fully flexible, neutral system.
 */
export struct PhononModes
{
  double3 kFractional{};
  std::vector<double> eigenvalues;  ///< Squared angular frequencies (ascending), internal units.
};

/** Per-atom inverse square-root masses in force-constant site order (flexible framework atoms, then molecules). */
export std::vector<double> phononInverseSqrtMasses(const System& system);

/** Assemble the complex Hermitian dynamical matrix D(k) (row-major, (3N)x(3N)) from the real-space force constants. */
export std::vector<std::complex<double>> computeDynamicalMatrix(const ForceConstants& forceConstants,
                                                                std::span<const double> inverseSqrtMass,
                                                                double3 kFractional);

/**
 * Reciprocal-space (Ewald Fourier) contribution to the mass-weighted dynamical matrix D(k).
 *
 * Returns the long-ranged Coulomb term that cannot be represented by finite-range real-space force
 * constants: a Bloch-summed lattice Fourier series (pair term) plus the k-independent self term that
 * enforces the acoustic sum rule. The site ordering matches computeRealSpaceForceConstants (flexible
 * framework atoms first, then molecule atoms). When the force field does not use Ewald Fourier sums,
 * the result is all zeros.
 *
 * Only the reciprocal part is produced here; the short-ranged erfc real-space Coulomb and the
 * intramolecular exclusion correction (molecule pairs and flexible-framework bonded pairs) are already
 * part of the real-space force constants. Charges of a rigid framework are fixed, so they add no matrix
 * rows/columns but enter the molecule self term as a static structure factor.
 */
export std::vector<std::complex<double>> computeEwaldReciprocalDynamicalMatrix(const System& system,
                                                                               std::span<const double> inverseSqrtMass,
                                                                               double3 kFractional);

/** Assemble and diagonalize the real-space-only D(k), returning the squared angular frequencies. */
export PhononModes computePhononModes(const ForceConstants& forceConstants, std::span<const double> inverseSqrtMass,
                                      double3 kFractional);

/** Assemble and diagonalize D(k) including the Ewald reciprocal contribution for `system`. */
export PhononModes computePhononModes(const System& system, const ForceConstants& forceConstants,
                                      std::span<const double> inverseSqrtMass, double3 kFractional);

/**
 * Convenience: assemble the force constants and masses once, then evaluate the modes along a k-path.
 *
 * For fully flexible systems the Cartesian dynamical matrix is used, including the Ewald reciprocal
 * contribution when the force field uses charges. When the system contains rigid molecules the modes are
 * built in generalized center-of-mass + orientation coordinates: the Cartesian force constants are
 * projected with the rigid-body kinematic Jacobian, the on-site gradient-curvature term is added, and the
 * matrix is mass-weighted with the inertia-tensor metric (so at k = 0 the result matches computeNormalModes).
 * Charged rigid molecules are supported: the reciprocal-space Ewald matrix is added to the Cartesian force
 * constants (as an ordinary position-dependent potential) before the same projection and mass weighting.
 */
export std::vector<PhononModes> computePhononDispersion(const System& system, std::span<const double3> kPath);

/** Convert squared angular frequencies to wavenumbers [cm^-1]; imaginary modes are returned as negatives. */
export std::vector<double> phononFrequenciesWavenumber(std::span<const double> eigenvalues);
