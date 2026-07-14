module;

export module potential_gradient_vdw;

import std;

import double4;

import vdwparameters;
import forcefield;
import gradient_factor;
import potential_vdw_rare_derivatives;

export namespace Potentials
{
/**
 * \brief Computes the gradient of the van der Waals (VDW) potential.
 *
 * This function calculates the gradient of the VDW potential between two atom types.
 * It returns D[U[r], r] / r to avoid computing the square root for Lennard-Jones (LJ) potential,
 * as only the squared distance (rr) is required.
 *
 * The returned GradientFactor.dUdlambda holds the symmetric derivative factor X such that
 *   dU/d(scalingA) = scalingB * X   and   dU/d(scalingB) = scalingA * X.
 *
 * \param forcefield The force field parameters defining the interaction.
 * \param scalingA Scaling factor for atom A.
 * \param scalingB Scaling factor for atom B.
 * \param rr The squared distance between the two atoms.
 * \param typeA The type identifier for atom A.
 * \param typeB The type identifier for atom B.
 *
 * \return A ForceFactor object containing the computed forces.
 */
[[clang::always_inline]] inline GradientFactor potentialVDWGradient(const ForceField& forcefield,
                                                                    const double& scalingA, const double& scalingB,
                                                                    const double& rr, const std::size_t& typeA,
                                                                    const std::size_t& typeB)
{
  VDWParameters::Type potentialType = forcefield(typeA, typeB).type;

  double scaling = scalingA * scalingB;
  switch (potentialType)
  {
    [[likely]] case VDWParameters::Type::LennardJones:
    {
      double arg1 = 4.0 * forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
      double arg3 = forcefield(typeA, typeB).shift;
      double temp = (rr / arg2);
      double temp3 = temp * temp * temp;
      double inv_scaling = 1.0 - scaling;
      double rri3 = 1.0 / (temp3 + 0.5 * inv_scaling * inv_scaling);
      double rri6 = rri3 * rri3;
      double term = arg1 * (rri3 * (rri3 - 1.0)) - arg3;
      double dlambda_term = arg1 * scaling * inv_scaling * (2.0 * rri6 * rri3 - rri6);
      return GradientFactor(scaling * term, term + dlambda_term,
                            12.0 * scaling * arg1 * (rri6 * temp3 * (0.5 - rri3)) / rr);
    }
    case VDWParameters::Type::LennardJonesShiftedForce:
    {
      const VDWParameters& p = forcefield(typeA, typeB);
      double eps4 = 4.0 * p.parameters.x;
      double sigma2 = p.parameters.y * p.parameters.y;
      double sigma6 = sigma2 * sigma2 * sigma2;
      double c6 = p.parameters2.x;
      double rc = p.parameters2.y;
      double linearCoefficient = 12.0 * c6 * c6 - 6.0 * c6;
      double invScaling = 1.0 - scaling;
      double x6 = rr * rr * rr + 0.5 * invScaling * invScaling * p.parameters2.w;
      double rs = std::sqrt(std::cbrt(x6));
      double u6 = sigma6 / x6;
      double u12 = u6 * u6;
      double term =
          eps4 * (u12 - u6 - c6 * (c6 - 1.0) + linearCoefficient * (rs - rc) / rc);
      double deriv = eps4 * (12.0 * u12 - 6.0 * u6 - linearCoefficient * rs / rc);
      double dlambdaTerm = scaling * invScaling * p.parameters2.w * deriv / (6.0 * x6);
      double gradientFactor = -scaling * deriv * rr * rr / x6;
      return GradientFactor(scaling * term, term + dlambdaTerm, gradientFactor);
    }
    case VDWParameters::Type::LennardJonesSecondOrderTaylorShifted:
    {
      const VDWParameters& p = forcefield(typeA, typeB);
      double eps4 = 4.0 * p.parameters.x;
      double sigma2 = p.parameters.y * p.parameters.y;
      double sigma6 = sigma2 * sigma2 * sigma2;
      double c6 = p.parameters2.x;
      double rc = p.parameters2.y;
      double linearCoefficient = 12.0 * c6 * c6 - 6.0 * c6;
      double quadraticCoefficient = 156.0 * c6 * c6 - 42.0 * c6;
      double invScaling = 1.0 - scaling;
      double x6 = rr * rr * rr + 0.5 * invScaling * invScaling * p.parameters2.w;
      double rs = std::sqrt(std::cbrt(x6));
      double displacement = (rs - rc) / rc;
      double u6 = sigma6 / x6;
      double u12 = u6 * u6;
      double term =
          eps4 * (u12 - u6 - c6 * (c6 - 1.0) + linearCoefficient * displacement -
                  0.5 * quadraticCoefficient * displacement * displacement);
      double deriv =
          eps4 * (12.0 * u12 - 6.0 * u6 - linearCoefficient * rs / rc +
                  quadraticCoefficient * rs * (rs - rc) / (rc * rc));
      double dlambdaTerm = scaling * invScaling * p.parameters2.w * deriv / (6.0 * x6);
      double gradientFactor = -scaling * deriv * rr * rr / x6;
      return GradientFactor(scaling * term, term + dlambdaTerm, gradientFactor);
    }
    case VDWParameters::Type::RepulsiveHarmonic:
    {
      Detail::RareVDWDerivatives derivatives =
          Detail::evaluateRareVDWDerivatives(forcefield, scalingA, scalingB, rr, typeA, typeB, potentialType);
      double gradientFactor =
          rr > 0.0 ? derivatives.radialFirstDerivative / std::sqrt(rr) : 0.0;
      return GradientFactor(derivatives.energy, derivatives.dUdlambda, gradientFactor);
    }
    default:
    {
      Detail::RareVDWDerivatives derivatives =
          Detail::evaluateRareVDWDerivatives(forcefield, scalingA, scalingB, rr, typeA, typeB, potentialType);
      double gradientFactor =
          rr > 0.0 ? derivatives.radialFirstDerivative / std::sqrt(rr) : 0.0;
      return GradientFactor(derivatives.energy, derivatives.dUdlambda, gradientFactor);
    }
  }
};
}  // namespace Potentials
