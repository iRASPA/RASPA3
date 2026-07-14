module;

export module potential_correction_vdw;

import std;

import vdwparameters;
import double4;

export namespace Potentials
{
/**
 * \brief Calculates the energy tail-correction integral for van der Waals (VDW) interactions.
 *
 * Returns Integrate[U[r]*r^2, {r, rc, Infinity}] for the fully-coupled, unshifted potential.
 * The result is multiplied elsewhere by 2 pi and the number densities of the two atom types.
 * Derived constants (e.g. the Feynmann-Hibbs pre-factor) must have been computed via
 * VDWParameters::computeDerivedParameters() before calling this function.
 *
 * \param parameters The van der Waals parameters of the pair interaction.
 * \param cutOffVDW The cutoff distance.
 *
 * \return The calculated VDW energy tail-correction integral.
 */
inline double potentialCorrectionVDW(const VDWParameters &parameters, double cutOffVDW)
{
  double rc = cutOffVDW;
  switch (parameters.type)
  {
    case VDWParameters::Type::LennardJones:
    {
      double arg1 = parameters.parameters.x;
      double arg2 = parameters.parameters.y;
      double term3 = (arg2 / rc) * (arg2 / rc) * (arg2 / rc);
      double term6 = term3 * term3;
      return (4.0 / 3.0) * arg1 * arg2 * arg2 * arg2 * ((1.0 / 3.0) * term6 * term3 - term3);
    }
    case VDWParameters::Type::FeynmannHibbs:
    {
      double arg1 = parameters.parameters.x;
      double arg2 = parameters.parameters.y;
      double fh = parameters.parameters2.x;
      double term3 = (arg2 / rc) * (arg2 / rc) * (arg2 / rc);
      double term9 = term3 * term3 * term3;
      double term6 = term3 * term3;
      double term12 = term6 * term6;
      return (4.0 / 3.0) * arg1 * arg2 * arg2 * arg2 * (term9 / 3.0 - term3) +
             fh * 24.0 * arg1 * rc * (2.0 * term12 - term6);
    }
    case VDWParameters::Type::BuckingHam:
    {
      double A = parameters.parameters.x;
      double b = parameters.parameters.y;
      double C6 = parameters.parameters.z;
      double br = b * rc;
      return A * std::exp(-br) * (2.0 + br * (2.0 + br)) / (b * b * b) - C6 / (3.0 * rc * rc * rc);
    }
    case VDWParameters::Type::Morse:
    {
      double D = parameters.parameters.x;
      double a = parameters.parameters.y;
      double r0 = parameters.parameters.z;
      double ar = a * rc;
      return (D * std::exp(a * (r0 - 2.0 * rc)) *
              (std::exp(a * r0) * (1.0 + 2.0 * ar * (1.0 + ar)) -
               8.0 * std::exp(ar) * (2.0 + ar * (2.0 + ar)))) /
             (4.0 * a * a * a);
    }
    case VDWParameters::Type::MM3:
    {
      double epsilon = parameters.parameters.x;
      double sigma = parameters.parameters.y;
      double sigma2 = sigma * sigma;
      double sigma6 = sigma2 * sigma2 * sigma2;
      return (1.0 / 864.0) *
                 (epsilon * sigma * 1.84e5 * std::exp(-12.0 * rc / sigma) *
                  (sigma2 + 12.0 * sigma * rc + 72.0 * rc * rc)) -
             epsilon * 2.25 * sigma6 / (3.0 * rc * rc * rc);
    }
    case VDWParameters::Type::BornHugginsMeyer:
    {
      double A = parameters.parameters.x;
      double b = parameters.parameters.y;
      double sigma = parameters.parameters.z;
      double C6 = parameters.parameters.w;
      double C8 = parameters.parameters2.x;
      double br = b * rc;
      double rc2 = rc * rc;
      double rc5 = rc2 * rc2 * rc;
      return (A * std::exp(b * (sigma - rc)) * (2.0 + br * (2.0 + br))) / (b * b * b) -
             (3.0 * C8 + 5.0 * C6 * rc2) / (15.0 * rc5);
    }
    case VDWParameters::Type::Potential12_6:
    {
      double p0 = parameters.parameters.x;
      double p1 = parameters.parameters.y;
      double rc3 = rc * rc * rc;
      double rc9 = rc3 * rc3 * rc3;
      return (p0 / rc9 - 3.0 * p1 / rc3) / 9.0;
    }
    case VDWParameters::Type::Potential12_6_2_0:
    {
      throw std::runtime_error("The 12-6-2-0 potential has no converging tail-correction.\n");
    }
    case VDWParameters::Type::CFF9_6:
    {
      double p0 = parameters.parameters.x;
      double p1 = parameters.parameters.y;
      double rc3 = rc * rc * rc;
      double rc6 = rc3 * rc3;
      return (p0 - 2.0 * p1 * rc3) / (6.0 * rc6);
    }
    case VDWParameters::Type::CFFEpsilonSigma:
    {
      double epsilon = parameters.parameters.x;
      double sigma = parameters.parameters.y;
      double sigma3 = sigma * sigma * sigma;
      double sigma6 = sigma3 * sigma3;
      double rc3 = rc * rc * rc;
      double rc6 = rc3 * rc3;
      return epsilon * sigma6 * (sigma3 - 3.0 * rc3) / (3.0 * rc6);
    }
    case VDWParameters::Type::MatsuokaClementiYoshimine:
    {
      double A = parameters.parameters.x;
      double a = parameters.parameters.y;
      double B = parameters.parameters.z;
      double b = parameters.parameters.w;
      double ar = a * rc;
      double br = b * rc;
      return A * std::exp(-ar) * (2.0 + ar * (2.0 + ar)) / (a * a * a) +
             B * std::exp(-br) * (2.0 + br * (2.0 + br)) / (b * b * b);
    }
    case VDWParameters::Type::Generic:
    {
      double A = parameters.parameters.x;
      double b = parameters.parameters.y;
      double C4 = parameters.parameters.z;
      double C6 = parameters.parameters.w;
      double C8 = parameters.parameters2.x;
      double C10 = parameters.parameters2.y;
      double br = b * rc;
      double rc3 = rc * rc * rc;
      double rc5 = rc3 * rc * rc;
      double rc7 = rc5 * rc * rc;
      return A * std::exp(-br) * (2.0 + br * (2.0 + br)) / (b * b * b) - C4 / rc - C6 / (3.0 * rc3) -
             C8 / (5.0 * rc5) - C10 / (7.0 * rc7);
    }
    case VDWParameters::Type::PellenqNicholson:
    {
      double A = parameters.parameters.x;
      double b = parameters.parameters.y;
      double C6 = parameters.parameters.z;
      double C8 = parameters.parameters.w;
      double C10 = parameters.parameters2.x;
      double br = b * rc;
      double rc3 = rc * rc * rc;
      double rc5 = rc3 * rc * rc;
      double rc7 = rc5 * rc * rc;
      double expTerm = std::exp(-br);
      double term1 = A * expTerm * (2.0 + br * (2.0 + br)) / (b * b * b);
      double term2 = (C6 * (-240.0 + 240.0 * std::exp(br) -
                            br * (240.0 + br * (120.0 + br * (38.0 + br * (8.0 + br)))))) /
                     (720.0 * rc3);
      double term3 =
          (C8 * (-8064.0 + 8064.0 * std::exp(br) -
                 br * (8064.0 + br * (4032.0 + br * (1344.0 + br * (336.0 + br * (66.0 + br * (10.0 + br)))))))) /
          (40320.0 * rc5);
      double term4 =
          (C10 *
           (-518400.0 + 518400.0 * std::exp(br) -
            br * (518400.0 +
                  br * (259200.0 +
                        br * (86400.0 +
                              br * (21600.0 + br * (4320.0 + br * (720.0 + br * (102.0 + br * (12.0 + br)))))))))) /
          (3.6288e6 * rc7);
      return term1 - expTerm * (term2 + term3 + term4);
    }
    case VDWParameters::Type::HydratedIonWater:
    {
      double A = parameters.parameters.x;
      double b = parameters.parameters.y;
      double C4 = parameters.parameters.z;
      double C6 = parameters.parameters.w;
      double C12 = parameters.parameters2.x;
      double br = b * rc;
      double rc3 = rc * rc * rc;
      double rc9 = rc3 * rc3 * rc3;
      return A * std::exp(-br) * (2.0 + br * (2.0 + br)) / (b * b * b) - C4 / rc - C6 / (3.0 * rc3) -
             C12 / (9.0 * rc9);
    }
    case VDWParameters::Type::Mie:
    {
      double p0 = parameters.parameters.x;
      double n = parameters.parameters.y;
      double p2 = parameters.parameters.z;
      double m = parameters.parameters.w;
      return p0 * std::pow(rc, 3.0 - n) / (n - 3.0) - p2 * std::pow(rc, 3.0 - m) / (m - 3.0);
    }
    case VDWParameters::Type::LennardJonesShiftedForce:
    case VDWParameters::Type::LennardJonesSecondOrderTaylorShifted:
    case VDWParameters::Type::WeeksChandlerAndersen:
      // zero at (and beyond) the cutoff by construction
      return 0.0;
    default:
      return 0.0;
  }
};
}  // namespace Potentials
