module;

export module phonon_force_constants;

import std;

import int3;
import system;

/**
 * Image-resolved interatomic force constants in Cartesian atomic coordinates.
 *
 * Stores blocks Phi_{i,j}(R) = d^2 U / du_i(0) du_j(R), where atom i sits in the
 * home cell and atom j sits in the periodic image displaced by the lattice
 * vector R (integer multiples of the cell vectors). The blocks are the seed for
 * the dynamical matrix D(k) = sum_R Phi(R) exp(i k.R); summing over R recovers
 * the Gamma-point (folded) Cartesian Hessian.
 *
 * The container is dense per image: each key R owns a full (3N)x(3N) matrix in
 * row-major order. This is memory-heavy for large systems and is intended as the
 * correctness-first foundation for phonon dispersion.
 */
export class ForceConstants
{
 public:
  ForceConstants() = default;
  explicit ForceConstants(std::size_t numberOfAtoms);

  std::size_t numberOfAtoms() const noexcept { return _numberOfAtoms; }
  std::size_t dimension() const noexcept { return _dimension; }
  std::size_t numberOfImages() const noexcept { return _blocks.size(); }

  /** Accumulate a 3x3 Cartesian block (row-major) at atom pair (i, j) and lattice vector R. */
  void addBlock(std::size_t i, std::size_t j, int3 latticeVector, const std::array<double, 9>& block);

  /** Accumulate a full dense (3N)x(3N) row-major matrix into the block for lattice vector R. */
  void addToImage(int3 latticeVector, std::span<const double> matrix);

  using BlockMap = std::unordered_map<int3, std::vector<double>, int3::hashFunction>;

  const BlockMap& blocks() const noexcept { return _blocks; }

  /** Sum over all lattice vectors R, yielding the folded (Gamma-point) Cartesian Hessian, row-major (3N)x(3N). */
  std::vector<double> foldToGamma() const;

 private:
  std::size_t _numberOfAtoms{};
  std::size_t _dimension{};
  BlockMap _blocks;
};

/**
 * Assemble the real-space (van der Waals + real-space Coulomb) force constants for the
 * flexible atoms of `system`, resolved by periodic image.
 *
 * Supported contributions:
 *  - molecule-molecule intermolecular pairs,
 *  - framework-molecule pairs (rigid framework: only the molecule-side home block; flexible
 *    framework: the full framework/molecule self and cross blocks),
 *  - framework-framework intramolecular terms for a flexible framework (bonds, bends, torsions,
 *    improper torsions and the scaled/nonbonded intra van der Waals and Coulomb pairs),
 *  - molecule intramolecular terms (bonds, bends, urey-bradley, torsions, cross terms and the intra
 *    van der Waals / Coulomb pairs), which are image-free and land in the home cell Phi(0),
 *  - the Ewald intramolecular exclusion correction (-C q_i q_j erf(alpha r)/r for every intramolecular
 *    pair), whose short-ranged Hessian is folded into the image-resolved blocks when the force field
 *    uses Ewald charges.
 *
 * Global atom ordering matches the minimization DOF layout: flexible framework atoms first, then
 * molecule atoms. Rigid molecules and polarization are rejected with an exception so results are never
 * silently partial. The long-ranged reciprocal-space Ewald contribution is NOT part of the real-space
 * force constants; it is added separately at each k-point (see computeEwaldReciprocalDynamicalMatrix).
 *
 * Summed over images, the result reproduces the position-position block of the analytic generalized
 * Hessian for such a system (see the fold-consistency tests).
 */
export ForceConstants computeRealSpaceForceConstants(const System& system);
