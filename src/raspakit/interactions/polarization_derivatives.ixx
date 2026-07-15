module;

export module interactions_polarization_derivatives;

import std;

import int3;
import double3;
import double3x3;
import system;

/**
 * Analytic first and second derivatives of the (non-self-consistent) polarization energy
 *
 *   U_pol = -1/2 sum_A alpha_A |E_A|^2,
 *
 * where E_A is the electric field on molecule atom A produced by all permanent charges. The field is
 * built exactly as System::computeTotalElectricField:
 *   - real-space molecule-molecule Coulomb (different molecules only, cut off at cutOffCoulomb),
 *   - real-space framework-molecule Coulomb,
 *   - the Ewald reciprocal field of the fixed (rigid) framework structure factor (an on-site term that
 *     depends only on the field-point position).
 * The intramolecular Ewald exclusion does NOT contribute to the field in this model.
 *
 * The derivatives w.r.t. Cartesian atomic displacements are
 *   dU/dr_m         = - sum_A alpha_A E_A . (dE_A/dr_m),
 *   d^2U/dr_m dr_n  = - sum_A alpha_A [ (dE_A/dr_m)^T (dE_A/dr_n) + E_A . (d^2 E_A/dr_m dr_n) ].
 * The first (field-gradient outer product) term is many-body: it couples any two atoms that both feed
 * the field on a common atom A. The reciprocal contribution enters only the on-site (A,A) block.
 */
export namespace Interactions
{
/**
 * Cartesian polarization derivatives in global-atom order (framework atoms first, then molecule atoms).
 * Blocks are 3x3 row-major matrices d^2U/dr_i dr_j keyed by the atom pair (i, j); only pairs between
 * movable atoms (plus movable self-blocks) are populated.
 *
 * `hessianBlocks` is the minimum-image-folded Hessian (correct for the Gamma point / normal modes).
 * `imageHessianBlocks` additionally resolves the periodic image: it keys the same 3x3 blocks by the
 * lattice vector R connecting the two atoms (field point in the home cell, source in image R) inside an
 * outer map. Summing over R reproduces `hessianBlocks`; keeping R lets the phonon force constants carry
 * the Bloch phase e^{2 pi i k.R} so the dynamical matrix is correct at every k, not only at Gamma.
 *
 * Strain derivatives (populated only when `strainBases` is non-empty) follow the method of homogeneous
 * deformation with the log-strain convention used by the cell degrees of freedom: for strain amplitudes
 * s_a along the (symmetric) generator matrices B_a the cell/positions deform as F = exp(sum_a s_a B_a),
 * reciprocal vectors k -> F^{-T} k and the volume as det(F). Rigid-molecule atoms deform via their center
 * of mass (y = F COM + sigma with the body offset sigma strain-independent), so separations use the COM-COM
 * arm and the reciprocal phase of a rigid field point carries the extra body-offset terms. The outputs are
 *   - `cellGradient[a]`       = dU/ds_a,
 *   - `positionStrain[i*nB+a]`= d(grad_i)/ds_a  (Cartesian *pure* force change, per global atom i),
 *   - `strainStrain[a*nB+b]`  = d^2U/ds_a ds_b.
 * `positionStrain` is the pure force change under strain (atoms moving via their COM); the caller projects
 * it onto the generalized DOFs and adds the kinematic term B_a . grad on the translational DOFs (which
 * vanishes for rigid-orientation DOFs). Both flexible and rigid molecules and frameworks are supported.
 */
struct PolarizationDerivatives
{
  double energy{};
  std::size_t numberOfFrameworkAtoms{};
  std::size_t numberOfMoleculeAtoms{};
  std::vector<double3> gradient;
  std::map<std::array<std::size_t, 2>, std::array<double, 9>> hessianBlocks;
  std::unordered_map<int3, std::map<std::array<std::size_t, 2>, std::array<double, 9>>, int3::hashFunction>
      imageHessianBlocks;

  std::size_t numberOfStrainBases{};
  std::vector<double> cellGradient;     // [a]
  std::vector<double3> positionStrain;  // [globalAtom * numberOfStrainBases + a]
  std::vector<double> strainStrain;     // [a * numberOfStrainBases + b]
};

/**
 * Compute the analytic polarization gradient and Cartesian Hessian for `system`.
 *
 * `movable` is indexed in global-atom order (framework atoms first, then molecule atoms) and flags the
 * atoms that carry degrees of freedom; blocks and gradient entries are produced only for movable atoms,
 * while the field itself is built from every charged atom (fixed framework included). Polarizable atoms
 * are required to be movable.
 *
 * When `strainBases` is non-empty the analytic cell-strain derivatives are computed as well (see above);
 * `strainBases` holds the (symmetric) strain generator matrices B_a of the active cell degrees of freedom.
 */
PolarizationDerivatives computePolarizationDerivatives(const System& system, std::span<const std::uint8_t> movable,
                                                       std::span<const double3x3> strainBases = {});
}  // namespace Interactions
