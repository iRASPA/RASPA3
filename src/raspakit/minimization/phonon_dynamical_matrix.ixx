module;

export module phonon_dynamical_matrix;

import std;

import int3;
import double3;
import system;
import phonon_force_constants;
import phonon_kpath;

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

/** Phonon band structure along a labeled high-symmetry path: sampled k-points and their modes (parallel). */
export struct PhononDispersionResult
{
  std::vector<PhononKPoint> path;   ///< Sampled k-points (fractional, path coordinate, labels).
  std::vector<PhononModes> modes;   ///< Modes at each path point, parallel to `path`.
};

/**
 * Sample the phonon dispersion along the band-structure path defined by `nodes`.
 *
 * The nodes are expanded into a dense path with `pointsPerSegment` steps per segment (see
 * buildPhononKPath), then computePhononDispersion is evaluated at every sampled k-point. The k-point
 * metadata (path coordinate, labels) is retained alongside the squared frequencies for plotting.
 */
export PhononDispersionResult computePhononDispersionAlongPath(const System& system,
                                                               std::span<const PhononPathNode> nodes,
                                                               std::size_t pointsPerSegment);

/** Human-readable band-structure table (path coordinate, k-point, and per-band wavenumbers). */
export std::string writePhononDispersion(const PhononDispersionResult& result);

/** Phonon density of states obtained by integrating the modes over a uniform Brillouin-zone q-mesh. */
export struct PhononDensityOfStates
{
  int3 mesh{};                     ///< Gamma-centered Monkhorst-Pack mesh dimensions used for the integration.
  std::size_t numberOfQPoints{};   ///< Total number of sampled q-points (mesh.x * mesh.y * mesh.z).
  std::vector<double> frequency;   ///< Bin-center frequencies [cm^-1] (or reduced units).
  std::vector<double> dos;         ///< Density of states; integral over frequency equals the number of branches 3N.
};

/**
 * Compute the phonon density of states on a Gamma-centered uniform q-mesh over the Brillouin zone.
 *
 * The dynamical matrix is diagonalized at every q-point of the `mesh` (using the same machinery as the
 * dispersion), the resulting frequencies are converted to wavenumbers, and each frequency is accumulated
 * onto a linear frequency grid of `numberOfBins` bins with a normalized Gaussian of width `broadening`
 * (same units as the frequency axis). The DOS is normalized by the number of q-points so that its integral
 * over frequency equals the number of branches (3N). Because a uniform grid over the reciprocal unit cell is
 * an unbiased sampling of the Brillouin zone, this is the quantity to compare with an inelastic-neutron
 * generalized phonon DOS. Intended to be evaluated on the primitive-cell system (unfolded spectrum).
 *
 * \param system System carrying the (minimized) framework/molecules and force field.
 * \param mesh Gamma-centered q-mesh dimensions; each component is clamped to at least 1.
 * \param numberOfBins Number of frequency bins in the histogram.
 * \param broadening Gaussian broadening (standard deviation) on the frequency axis; if <= 0 the bin width is used.
 */
export PhononDensityOfStates computePhononDensityOfStates(const System& system, int3 mesh, std::size_t numberOfBins,
                                                          double broadening);

/** Human-readable two-column density-of-states table (frequency, DOS). */
export std::string writePhononDensityOfStates(const PhononDensityOfStates& result);
