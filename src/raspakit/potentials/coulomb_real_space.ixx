module;

export module potential_coulomb_real_space;

import std;

import forcefield;

export namespace Potentials
{
struct CoulombRealSpaceFactors
{
  double potential;
  double firstDerivativeFactor;
  double secondDerivativeFactor;
  double thirdDerivativeFactor;
};

// Distance offset used in the scaled (CFCMC fractional) Ewald real-space and exclusion terms:
// the pair separation r is replaced by r + Q with Q = delta * (1 - scalingA * scalingB) and
// delta = 0.01 Angstrom. For full interactions (scaling 1) the offset vanishes and plain Ewald is
// recovered; for scaled interactions the potential remains finite even at r = 0, removing the
// erfc(alpha r)/r divergence that otherwise produces 0 * inf = NaN for zero Coulomb scaling.
// See Hens et al., J. Chem. Inf. Model. 2020 (Brick-CFCMC), SI Eqs. (S13)-(S15) and (S19)-(S21).
inline constexpr double EwaldChargeOffsetDelta = 0.01;

/**
 * \brief Offset-form factors of the intramolecular Ewald exclusion term for scaled interactions.
 *
 * The exclusion energy of an intramolecular pair is U = -lambda_t * C q_A q_B * erf(alpha s)/s
 * with s = r + Q and Q = EwaldChargeOffsetDelta * (1 - lambda_t), lambda_t = scalingA * scalingB
 * (Eq. (S15) with erf = 1 - erfc). All fields are per Coulomb prefactor C q_A q_B:
 *   - potential: erf(alpha s)/s, so the exclusion energy is -lambda_t * C q_A q_B * potential.
 *   - dUdlambda: symmetric lambda-derivative factor, d(-U)/d(scalingA) = -scalingB * C q_A q_B * dUdlambda,
 *     including the chain-rule contribution of Q (Eq. (S21)).
 *   - firstDerivativeFactor: (1/r) d/dr [erf(alpha s)/s], so the gradient of U on atom A is
 *     -lambda_t * C q_A q_B * firstDerivativeFactor * dr with dr = posA - posB.
 *   - secondDerivativeFactor: (phi'' - phi'/r)/r^2 with phi(s) = erf(alpha s)/s, the RASPA Hessian
 *     factor of the (positive-erf) exclusion kernel; the Hessian factors of U carry a minus sign.
 */
struct EwaldExclusionFactors
{
  double potential;
  double dUdlambda;
  double firstDerivativeFactor;
  double secondDerivativeFactor;
};

[[clang::always_inline]] inline EwaldExclusionFactors ewaldExclusionFactors(double alpha, double scalingTotal,
                                                                            double r)
{
  const double offset = EwaldChargeOffsetDelta * (1.0 - scalingTotal);
  const double s = r + offset;
  const double inverseS = 1.0 / s;
  const double inverseR = 1.0 / r;
  const double potential = std::erf(alpha * s) * inverseS;
  const double gaussianTerm =
      2.0 * alpha * std::numbers::inv_sqrtpi_v<double> * std::exp(-alpha * alpha * s * s);
  // phi'(s) and phi''(s) of phi(s) = erf(alpha s)/s
  const double firstDerivative = (gaussianTerm - potential) * inverseS;
  const double secondDerivative =
      -2.0 * alpha * alpha * gaussianTerm - 2.0 * (gaussianTerm - potential) * inverseS * inverseS;
  return {potential, potential - scalingTotal * EwaldChargeOffsetDelta * firstDerivative,
          firstDerivative * inverseR, (secondDerivative - firstDerivative * inverseR) * inverseR * inverseR};
}

inline double coulombSelfEnergyPrefactor(const ForceField& forceField)
{
  const double alpha = forceField.EwaldAlpha;
  const double cutoff = forceField.cutOffCoulomb;
  const double alphaCutoff = alpha * cutoff;
  const double exponential = std::exp(-alphaCutoff * alphaCutoff);
  const double erfc = std::erfc(alphaCutoff);
  const double inverseSqrtPi = std::numbers::inv_sqrtpi_v<double>;

  switch (forceField.chargeMethod)
  {
    case ForceField::ChargeMethod::Wolf:
      return -alpha * inverseSqrtPi - 0.5 * erfc / cutoff;
    case ForceField::ChargeMethod::DampedShiftedForce:
      return -erfc / cutoff - alpha * (1.0 + exponential) * inverseSqrtPi;
    case ForceField::ChargeMethod::ModifiedShiftedForce:
    {
      const double beta = forceField.modifiedShiftedForceBeta;
      const double shortRangePrime =
          -erfc / (cutoff * cutoff) - 2.0 * alpha * exponential * inverseSqrtPi / cutoff;
      const double longRangePrime =
          2.0 * alpha * exponential * inverseSqrtPi / cutoff - std::erf(alphaCutoff) / (cutoff * cutoff);
      return -0.5 * (erfc / cutoff - shortRangePrime / beta + 2.0 * alpha * inverseSqrtPi -
                     longRangePrime * std::exp(-beta * cutoff) / beta);
    }
    case ForceField::ChargeMethod::ZeroDipole:
      return (-alphaCutoff * (1.0 + 0.5 * exponential) * inverseSqrtPi - 0.75 * erfc) / cutoff;
    case ForceField::ChargeMethod::Ewald:
    case ForceField::ChargeMethod::Coulomb:
      return 0.0;
  }
  std::unreachable();
}

[[clang::always_inline]] inline CoulombRealSpaceFactors coulombRealSpaceFactors(const ForceField& forceField,
                                                                                double r)
{
  const double alpha = forceField.EwaldAlpha;
  const double alphaSquared = alpha * alpha;
  const double rr = r * r;
  const double inverseR = 1.0 / r;
  const double inverseRR = inverseR * inverseR;
  const double erfc = std::erfc(alpha * r);
  const double exponential = std::exp(-alphaSquared * rr);
  const double twoAlphaInvSqrtPi = 2.0 * alpha * std::numbers::inv_sqrtpi_v<double>;

  const double base = erfc * inverseR;
  const double basePrime = -erfc * inverseRR - twoAlphaInvSqrtPi * exponential * inverseR;
  const double baseSecond =
      2.0 * erfc * inverseRR * inverseR +
      2.0 * twoAlphaInvSqrtPi * exponential * (inverseRR + alphaSquared);
  const double baseThird =
      -6.0 * erfc * inverseRR * inverseRR -
      twoAlphaInvSqrtPi * exponential *
          (4.0 * alphaSquared * alphaSquared * r + 4.0 * alphaSquared * inverseR + 6.0 * inverseR * inverseRR);

  // Convert raw radial derivatives (B', B'', B''') into the RASPA factor convention:
  //   f1 = B'/r, f2 = (B'' - B'/r)/r^2, f3 = (B''' - 3 B''/r + 3 B'/r^2)/r^3.
  const auto factorsFrom = [&](double value, double firstDerivative, double secondDerivative, double thirdDerivative)
  {
    return CoulombRealSpaceFactors{
        value, firstDerivative * inverseR, (secondDerivative - firstDerivative * inverseR) * inverseRR,
        (thirdDerivative - 3.0 * secondDerivative * inverseR + 3.0 * firstDerivative * inverseRR) * inverseR *
            inverseRR};
  };

  if (forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
  {
    return {base, basePrime * inverseR, (baseSecond - basePrime * inverseR) * inverseRR,
            (baseThird - 3.0 * baseSecond * inverseR + 3.0 * basePrime * inverseRR) * inverseR * inverseRR};
  }

  const double cutoff = forceField.cutOffCoulomb;
  const double inverseCutoff = 1.0 / cutoff;
  const double inverseCutoffSquared = inverseCutoff * inverseCutoff;
  const double alphaCutoff = alpha * cutoff;
  const double cutoffExponential = std::exp(-alphaCutoff * alphaCutoff);
  const double cutoffErfc = std::erfc(alphaCutoff);
  const double cutoffBase = cutoffErfc * inverseCutoff;
  const double cutoffPrime =
      -cutoffErfc * inverseCutoffSquared - twoAlphaInvSqrtPi * cutoffExponential * inverseCutoff;

  double potential = base;
  double firstDerivative = basePrime;
  double secondDerivative = baseSecond;
  double thirdDerivative = baseThird;

  switch (forceField.chargeMethod)
  {
    case ForceField::ChargeMethod::Wolf:
      potential -= cutoffBase;
      break;
    case ForceField::ChargeMethod::DampedShiftedForce:
      potential -= cutoffBase + cutoffPrime * (r - cutoff);
      firstDerivative -= cutoffPrime;
      break;
    case ForceField::ChargeMethod::ModifiedShiftedForce:
    {
      const double beta = forceField.modifiedShiftedForceBeta;
      const double switching = std::exp(beta * (r - cutoff));
      potential -= cutoffBase + cutoffPrime * (switching - 1.0) / beta;
      firstDerivative -= cutoffPrime * switching;
      secondDerivative -= cutoffPrime * beta * switching;
      thirdDerivative -= cutoffPrime * beta * beta * switching;
      break;
    }
    case ForceField::ChargeMethod::ZeroDipole:
    {
      const double correction =
          cutoffErfc + 2.0 * alphaCutoff * cutoffExponential * std::numbers::inv_sqrtpi_v<double>;
      const double inverseCutoffCubed = inverseCutoffSquared * inverseCutoff;
      potential += 0.5 * correction * inverseCutoff * (rr * inverseCutoffSquared - 1.0) - cutoffBase;
      firstDerivative += correction * r * inverseCutoffCubed;
      secondDerivative += correction * inverseCutoffCubed;
      break;
    }
    case ForceField::ChargeMethod::Coulomb:
      potential = inverseR;
      firstDerivative = -inverseRR;
      secondDerivative = 2.0 * inverseRR * inverseR;
      thirdDerivative = -6.0 * inverseRR * inverseRR;
      break;
    case ForceField::ChargeMethod::Ewald:
      std::unreachable();
  }

  return factorsFrom(potential, firstDerivative, secondDerivative, thirdDerivative);
}
}  // namespace Potentials
