module;

export module potential_intra_pair;

import std;

import units;
import forcefield;
import potential_pair_derivatives;
import potential_pair_vdw;
import potential_pair_coulomb;
import potential_coulomb_real_space;

export namespace Potentials
{
/**
 * \brief Van der Waals interaction of a non-excluded intramolecular pair.
 *
 * The pair interacts with the regular force-field pair potential (the same truncated / shifted potential and the
 * same cutoff 'cutOffMoleculeVDW' as a pair of two different molecules) times the pair scaling (one for ordinary
 * pairs, the 1-4 scaling for the 1-4 pairs). The intramolecular interactions of a fractional molecule are not
 * scaled with lambda (its internal structure is independent of the coupling), so dUdlambda is zero.
 *
 * \param forceField The force field.
 * \param pairScaling The scaling of the pair (zero disables the pair).
 * \param rr The squared distance between the atoms.
 * \param typeA The pseudo-atom type of atom A.
 * \param typeB The pseudo-atom type of atom B.
 */
template <std::size_t Order>
[[clang::always_inline]] inline PairDerivatives<Order> intraMolecularVDW(const ForceField &forceField,
                                                                         double pairScaling, double rr,
                                                                         std::size_t typeA, std::size_t typeB)
{
  static_assert(Order <= 2, "intraMolecularVDW supports derivative orders 0, 1, and 2");

  if (pairScaling == 0.0 || rr >= forceField.cutOffMoleculeVDW * forceField.cutOffMoleculeVDW) return {};

  PairDerivatives<Order> derivatives = potentialVDW<Order>(forceField, 1.0, 1.0, rr, typeA, typeB);
  derivatives.energy *= pairScaling;
  derivatives.dUdlambda = 0.0;
  if constexpr (Order >= 1) derivatives.firstDerivativeFactor *= pairScaling;
  if constexpr (Order >= 2) derivatives.secondDerivativeFactor *= pairScaling;
  return derivatives;
}

/**
 * \brief Coulomb interaction of a non-excluded intramolecular pair.
 *
 * Inside the Coulomb cutoff the pair energy is
 *
 *     U = f C q_A q_B / r + lambda_A lambda_B C q_A q_B K(r)
 *
 * with f the pair scaling (one, or the 1-4 Coulomb scaling) and K the exclusion kernel of the charge method: for
 * Ewald K(r) = -erf(alpha s)/s (the Fourier sum counted this pair with the scaled charges, so its long-range part
 * is removed and the pair is left with the unscaled bare Coulomb interaction f/r); for the finite-cutoff shifted
 * methods (Wolf, damped-shifted-force, ...) K(r) = V(r) - 1/r, the completion of the shifted pair sum. For a
 * fully coupled ordinary pair (f = 1, lambda = 1) this is exactly the regular real-space pair potential, so the
 * same pair can be evaluated by a pair loop that treats it like any other pair. For the charge methods without
 * corrections (plain truncated Coulomb, Ewald without the Fourier part) U = f times the regular pair potential.
 *
 * Only the kernel term depends on lambda; 'dUdlambda' is the symmetric factor X with dU/d(lambda_A) = lambda_B X
 * (the convention of RunningEnergy::addDudlambdaEwald). Beyond the cutoff the pair does not interact (the Fourier
 * sum then supplies lambda_A lambda_B C q_A q_B erf(alpha r)/r, the long-range Ewald estimate of the pair).
 *
 * \param forceField The force field.
 * \param pairScaling The pair scaling f.
 * \param scalingA The Coulomb scaling of atom A.
 * \param scalingB The Coulomb scaling of atom B.
 * \param r The distance between the atoms.
 * \param chargeA The charge of atom A.
 * \param chargeB The charge of atom B.
 */
template <std::size_t Order>
[[clang::always_inline]] inline PairDerivatives<Order> intraMolecularCoulomb(const ForceField &forceField,
                                                                             double pairScaling, double scalingA,
                                                                             double scalingB, double r,
                                                                             double chargeA, double chargeB)
{
  static_assert(Order <= 2, "intraMolecularCoulomb supports derivative orders 0, 1, and 2");

  if (r >= forceField.cutOffCoulomb) return {};
  const double prefactor = Units::CoulombicConversionFactor * chargeA * chargeB;
  if (prefactor == 0.0) return {};

  const double scalingTotal = scalingA * scalingB;
  const double rr = r * r;
  const double inverseR = 1.0 / r;
  const double inverseR3 = inverseR / rr;
  const double inverseR5 = inverseR3 / rr;

  PairDerivatives<Order> result{};
  if (forceField.usesEwaldFourier())
  {
    const EwaldExclusionFactors exclusion = ewaldExclusionFactors(forceField.EwaldAlpha, scalingTotal, r);
    result.energy = prefactor * (pairScaling * inverseR - scalingTotal * exclusion.potential);
    result.dUdlambda = -prefactor * exclusion.dUdlambda;
    if constexpr (Order >= 1)
    {
      result.firstDerivativeFactor = -prefactor * (pairScaling * inverseR3 + scalingTotal * exclusion.firstDerivativeFactor);
    }
    if constexpr (Order >= 2)
    {
      result.secondDerivativeFactor =
          prefactor * (3.0 * pairScaling * inverseR5 - scalingTotal * exclusion.secondDerivativeFactor);
    }
    return result;
  }

  if (forceField.usesRealSpaceChargeCorrections())
  {
    const CoulombRealSpaceFactors factors = coulombRealSpaceFactors(forceField, r);
    const double completion = factors.potential - inverseR;
    result.energy = prefactor * (pairScaling * inverseR + scalingTotal * completion);
    result.dUdlambda = prefactor * completion;
    if constexpr (Order >= 1)
    {
      result.firstDerivativeFactor =
          prefactor * (-pairScaling * inverseR3 + scalingTotal * (factors.firstDerivativeFactor + inverseR3));
    }
    if constexpr (Order >= 2)
    {
      result.secondDerivativeFactor =
          prefactor * (3.0 * pairScaling * inverseR5 + scalingTotal * (factors.secondDerivativeFactor - 3.0 * inverseR5));
    }
    return result;
  }

  if (pairScaling == 0.0) return {};
  result = potentialCoulomb<Order>(forceField, 1.0, 1.0, r, chargeA, chargeB);
  result.energy *= pairScaling;
  result.dUdlambda = 0.0;
  if constexpr (Order >= 1) result.firstDerivativeFactor *= pairScaling;
  if constexpr (Order >= 2) result.secondDerivativeFactor *= pairScaling;
  return result;
}
}  // namespace Potentials
