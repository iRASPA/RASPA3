module;

export module interactions_hessian_ewald;

import std;

import atom;
import atom_dynamics;
import molecule;
import component;
import framework;
import forcefield;
import simulationbox;
import generalized_hessian;
import minimization_dof_layout;
import minimization_cell_layout;
import running_energy;

export namespace Interactions
{
/**
 * \brief Computes the Ewald Fourier energy, gradient, and generalized Hessian blocks.
 *
 * Reciprocal-space analogue of RASPA2's CalculateEwaldFourierDerivatives. The Fourier energy
 * U = sum_k factor(k) |S(k)|^2 with S(k) = sum_j q_j exp(i k.r_j) is differentiated exactly with
 * respect to the minimization degrees of freedom (flexible atom positions, rigid-molecule
 * center-of-mass and orientation) and the single isotropic strain coordinate:
 *
 *   dU/dtheta_a          = 2 factor Re[ i dS/dtheta_a conj(S) ]
 *   d2U/dtheta_a dtheta_b = 2 factor Re[ (dS/dtheta_a) conj(dS/dtheta_b) ] +
 *                           2 factor Re[ conj(S) d2S/dtheta_a dtheta_b ]
 *
 * assembled per wave vector from complex per-DOF accumulators. The strain convention matches the
 * pair terms: the cell and all positions scale with exp(epsilon), so for flexible atoms the phases
 * k.r are invariant and only the k-dependent prefactor contributes, while rigid molecules pick up
 * corrections through their fixed internal offsets (pos - com).
 *
 * The intramolecular exclusion correction -q_i q_j erf(alpha r)/r is included as a radial pair
 * term for flexible molecules (energy only for rigid molecules, whose internal distances are
 * constant). The self energy and the Bogusz net-charge correction are included in the energy.
 * The self term is constant for fixed alpha; the Bogusz term is position-independent but its
 * strain gradient and isotropic exp(epsilon) strain curvature are included.
 *
 * Per-atom Cartesian gradients are accumulated into \p moleculeDynamics / \p frameworkDynamics;
 * the strain first derivative is accumulated into hessian.strainGradient().
 */
RunningEnergy computeEwaldFourierHessian(
    const ForceField& forceField, const SimulationBox& simulationBox, const std::optional<Framework>& framework,
    std::span<const std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> fixedFrameworkStoredEik,
    double netChargeFramework, std::span<const Molecule> moleculeData, std::span<const Component> components,
    std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms, const MinimizationDofLayout& layout,
    GeneralizedHessian& hessian, std::span<AtomDynamics> moleculeDynamics,
    std::span<AtomDynamics> frameworkDynamics = {}, const CellMinimizationLayout& cellLayout = {});

}  // namespace Interactions
