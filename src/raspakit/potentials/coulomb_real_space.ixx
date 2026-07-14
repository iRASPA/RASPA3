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
};

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

  if (forceField.chargeMethod == ForceField::ChargeMethod::Ewald)
  {
    return {base, basePrime * inverseR, (baseSecond - basePrime * inverseR) * inverseRR};
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
      break;
    case ForceField::ChargeMethod::Ewald:
      std::unreachable();
  }

  return {potential, firstDerivative * inverseR,
          (secondDerivative - firstDerivative * inverseR) * inverseRR};
}
}  // namespace Potentials
