module;

export module interactions_pair_kernel;

import std;

import double3;
import double3x3;
import atom;
import simulationbox;
import forcefield;
export import potential_pair_derivatives;
import potential_pair_vdw;
import potential_pair_coulomb;

export namespace Interactions
{
/**
 * \brief Evaluates the VDW and Coulomb interaction of one atom pair at the requested derivative order.
 *
 * Shared pair kernel for the energy (Order 0), gradient (Order 1), and Hessian (Order 2) routines:
 * computes the minimum-image separation, applies the cutoffs, dispatches to the unified pair
 * potentials, and hands the resulting factors to the caller-supplied sinks. Each sink is invoked
 * with (const Potentials::PairDerivatives<Order>&, const double3& dr), where dr = posA - posB
 * after periodic boundary conditions; the Cartesian force on atom A is firstDerivativeFactor * dr.
 */
template <std::size_t Order, typename VDWSink, typename CoulombSink>
[[clang::always_inline]] inline void evaluatePair(const ForceField& forceField, const SimulationBox& simulationBox,
                                                  const Atom& atomA, const Atom& atomB, double cutOffVDWSquared,
                                                  double cutOffChargeSquared, bool useCharge, VDWSink&& vdwSink,
                                                  CoulombSink&& coulombSink)
{
  double3 dr = atomA.position - atomB.position;
  dr = simulationBox.applyPeriodicBoundaryConditions(dr);
  const double rr = double3::dot(dr, dr);

  if (rr < cutOffVDWSquared)
  {
    const Potentials::PairDerivatives<Order> factors =
        Potentials::potentialVDW<Order>(forceField, atomA.scalingVDW, atomB.scalingVDW, rr,
                                        static_cast<std::size_t>(atomA.type), static_cast<std::size_t>(atomB.type));
    vdwSink(factors, dr);
  }
  if (useCharge && rr < cutOffChargeSquared)
  {
    const double r = std::sqrt(rr);
    const Potentials::PairDerivatives<Order> factors = Potentials::potentialCoulomb<Order>(
        forceField, atomA.scalingCoulomb, atomB.scalingCoulomb, r, atomA.charge, atomB.charge);
    coulombSink(factors, dr);
  }
}

/**
 * \brief Loops over all intermolecular atom pairs and evaluates them at the requested derivative order.
 *
 * Triangular double loop with the molecule self-interaction skip (atoms sharing a molecule id do
 * not interact). Uses the molecule-molecule VDW cutoff. Sinks are invoked with
 * (std::size_t indexA, std::size_t indexB, const Atom& atomA, const Atom& atomB,
 *  const Potentials::PairDerivatives<Order>&, const double3& dr).
 */
template <std::size_t Order, typename VDWSink, typename CoulombSink>
inline void forEachMoleculeMoleculePair(const ForceField& forceField, const SimulationBox& simulationBox,
                                        std::span<const Atom> moleculeAtoms, VDWSink&& vdwSink,
                                        CoulombSink&& coulombSink)
{
  const bool useCharge = forceField.useCharge;
  const double cutOffVDWSquared = forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  if (moleculeAtoms.empty()) return;

  for (std::size_t indexA = 0; indexA + 1 < moleculeAtoms.size(); ++indexA)
  {
    const Atom& atomA = moleculeAtoms[indexA];
    for (std::size_t indexB = indexA + 1; indexB < moleculeAtoms.size(); ++indexB)
    {
      const Atom& atomB = moleculeAtoms[indexB];

      // skip interactions within the same molecule
      if (atomA.moleculeId == atomB.moleculeId) continue;

      evaluatePair<Order>(
          forceField, simulationBox, atomA, atomB, cutOffVDWSquared, cutOffChargeSquared, useCharge,
          [&](const Potentials::PairDerivatives<Order>& factors, const double3& dr)
          { vdwSink(indexA, indexB, atomA, atomB, factors, dr); },
          [&](const Potentials::PairDerivatives<Order>& factors, const double3& dr)
          { coulombSink(indexA, indexB, atomA, atomB, factors, dr); });
    }
  }
}

/**
 * \brief Loops over all framework atoms interacting with one molecule atom at the requested derivative order.
 *
 * Uses the framework-molecule VDW cutoff. dr = moleculeAtom.position - frameworkAtom.position, so
 * firstDerivativeFactor * dr is the force on the molecule atom. Sinks are invoked with
 * (std::size_t frameworkIndex, const Atom& frameworkAtom,
 *  const Potentials::PairDerivatives<Order>&, const double3& dr).
 */
template <std::size_t Order, typename VDWSink, typename CoulombSink>
inline void forEachFrameworkMoleculePair(const ForceField& forceField, const SimulationBox& simulationBox,
                                         const Atom& moleculeAtom, std::span<const Atom> frameworkAtoms,
                                         VDWSink&& vdwSink, CoulombSink&& coulombSink)
{
  const bool useCharge = forceField.useCharge;
  const double cutOffVDWSquared = forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::size_t index = 0; index < frameworkAtoms.size(); ++index)
  {
    const Atom& frameworkAtom = frameworkAtoms[index];
    evaluatePair<Order>(
        forceField, simulationBox, moleculeAtom, frameworkAtom, cutOffVDWSquared, cutOffChargeSquared, useCharge,
        [&](const Potentials::PairDerivatives<Order>& factors, const double3& dr)
        { vdwSink(index, frameworkAtom, factors, dr); },
        [&](const Potentials::PairDerivatives<Order>& factors, const double3& dr)
        { coulombSink(index, frameworkAtom, factors, dr); });
  }
}

/**
 * \brief Optional side-channel that lets the real-space Coulomb strain loops also gather the
 * polarization electric field and its cell-strain response for the molecular pressure.
 *
 * All spans are indexed by molecule atom (framework atoms act as field sources, never as field
 * points). For every real-space Coulomb pair the strain loops add the source contribution to the
 * field E_A and to the tensor
 *   fieldStrain_A[3*i + j][m] = dE_A[i] / dF[j][m],
 * the derivative of the field with respect to the cell deformation gradient F under molecular
 * center-of-mass scaling: separations deform with the COM-COM arm delta = d - sigma_A + sigma_source
 * (sigma = position - mass-weighted COM of the owning molecule; zero for framework sources). Atoms
 * with zero polarizability are skipped as field points.
 *
 * The caller completes the field with the Ewald reciprocal fixed-framework contribution and
 * contracts field and fieldStrain into the polarization energy -1/2 sum_A alpha_A |E_A|^2 and its
 * strain derivative -sum_A alpha_A E_A . dE_A/dF (Interactions::computePolarizationMolecularPressureStrain).
 */
struct PolarizationFieldStrain
{
  std::span<double3> field;                       ///< E_A per molecule atom.
  std::span<std::array<double3, 9>> fieldStrain;  ///< dE_A[i]/dF[j][m] stored at [3*i+j][m].
  std::span<const double3> centerOfMassOffset;    ///< sigma_A = r_A - COM(molecule of A), per molecule atom.
  std::span<const double> polarizability;         ///< Internal units; zero skips the field point.
};

/**
 * \brief Adds one real-space Coulomb source contribution to the gathered polarization field.
 *
 * d = fieldPoint - source (minimum image), delta = d - sigma_fieldPoint + sigma_source, and the
 * derivative factors are the *unit-charge* Coulomb factors of the pair; sourceCharge is the scaled
 * charge (scalingCoulomb * charge) of the source atom.
 */
inline void accumulatePolarizationFieldStrain(const PolarizationFieldStrain& gather, std::size_t index,
                                              double sourceCharge, const double3& d, const double3& delta,
                                              double firstDerivativeFactor, double secondDerivativeFactor)
{
  if (gather.polarizability[index] == 0.0) return;

  // E_A += -q_source * firstDerivativeFactor * d.
  gather.field[index] -= sourceCharge * firstDerivativeFactor * d;

  // M[i][j] = -q (f1 delta_ij + f2 d_i d_j) is dE_A/dd of this source; under strain the separation
  // moves with the COM arm, dd/dF[j][m] = e_j delta_m, so dE_A[i]/dF[j][m] += M[i][j] delta_m.
  std::array<double3, 9>& tensor = gather.fieldStrain[index];
  const std::array<double, 3> dComponents = {d.x, d.y, d.z};
  for (std::size_t i = 0; i < 3; ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      const double entry = -sourceCharge * (secondDerivativeFactor * dComponents[i] * dComponents[j] +
                                            (i == j ? firstDerivativeFactor : 0.0));
      tensor[3 * i + j] += entry * delta;
    }
  }
}

/**
 * \brief Accumulates the strain-derivative contribution gradient (x) displacement into the tensor.
 *
 * Column layout matches the existing strain conventions: tensor.{a,b,c}{x,y,z} +=
 * gradient.{x,y,z} * displacement.{x,y,z}.
 */
inline void accumulateStrainDerivative(double3x3& tensor, const double3& gradient, const double3& displacement)
{
  tensor.ax += gradient.x * displacement.x;
  tensor.bx += gradient.y * displacement.x;
  tensor.cx += gradient.z * displacement.x;

  tensor.ay += gradient.x * displacement.y;
  tensor.by += gradient.y * displacement.y;
  tensor.cy += gradient.z * displacement.y;

  tensor.az += gradient.x * displacement.z;
  tensor.bz += gradient.y * displacement.z;
  tensor.cz += gradient.z * displacement.z;
}
}  // namespace Interactions
