module;

export module normal_modes;

import std;

import system;

export struct NormalModesResult
{
  std::size_t numberOfModes{};
  /// Eigenvalues omega^2 of the mass-weighted Hessian in internal units, ascending.
  std::vector<double> eigenvalues;
  /// Column-major mass-weighted eigenvectors: eigenvectors[dof + numberOfModes * mode].
  std::vector<double> eigenvectors;
  std::size_t negativeModes{};
  std::size_t zeroModes{};
  /// Rigid-molecule orientation DOFs with (near-)zero moment of inertia, excluded from mass weighting.
  std::size_t discardedRotationalDofs{};
};

/**
 * Compute Gamma-point normal modes at the current (minimized) configuration.
 *
 * The generalized position-position Hessian is transformed with the inverse-square-root
 * mass metric: atomic masses for Cartesian DOFs, the molecular mass for rigid-body
 * center-of-mass DOFs, and the space-frame inertia tensor for rigid-body orientation
 * DOFs. The eigenvalues of the resulting dynamical matrix are the squared angular
 * frequencies omega^2; negative eigenvalues correspond to imaginary modes.
 *
 * Orientation directions with vanishing moment of inertia (linear or single-bead rigid
 * molecules) are excluded via a pseudo-inverse and appear as zero modes.
 */
export NormalModesResult computeNormalModes(const System& system, double relativeZeroTolerance = 1.0e-8);

/**
 * Signed mode frequencies: wavenumbers in cm^-1 (RASPA units) or angular frequency in
 * reduced units. Negative values denote imaginary modes.
 */
export std::vector<double> normalModeFrequencies(const NormalModesResult& result);

export std::string writeNormalModes(const NormalModesResult& result);

/**
 * Write one PDB movie per normal mode into `directory` (created if needed).
 *
 * Each atom oscillates as x(t) = x0 + A * e_a * sin(2 pi t), where e_a is the
 * mode's Cartesian displacement pattern and A scales the pattern so that the
 * largest atomic displacement equals `amplitude` (in Angstrom). Rigid-molecule
 * modes are animated through their center-of-mass and orientation degrees of
 * freedom. Files are named `mode_XXXX.s{systemIndex}.pdb`.
 *
 * \param numberOfPeriods    number of full oscillation periods per movie.
 * \param pointsPerPeriod    number of frames sampled within one period.
 * \param amplitude          maximum atomic displacement in Angstrom.
 */
export void writeNormalModeMovies(const System& system, const NormalModesResult& result, std::size_t systemIndex,
                                  std::size_t numberOfPeriods, std::size_t pointsPerPeriod, double amplitude,
                                  const std::filesystem::path& directory = "normal_modes");
