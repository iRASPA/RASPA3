module;

export module potential_correction_pressure;

import std;

import vdwparameters;
import double4;

export namespace Potentials
{
/**
 * \brief Calculates the pressure tail-correction integral for van der Waals (VDW) interactions.
 *
 * Returns Integrate[D[U[r],r]*r^3, {r, rc, Infinity}] for the fully-coupled, unshifted potential.
 * Derived constants (e.g. the Feynmann-Hibbs pre-factor) must have been computed via
 * VDWParameters::computeDerivedParameters() before calling this function.
 *
 * \param parameters The van der Waals parameters of the pair interaction.
 * \param cutOffVDW The cutoff distance.
 *
 * \return The calculated VDW pressure tail-correction integral.
 */
inline double potentialCorrectionPressure(const VDWParameters &parameters, double cutOffVDW)
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
      return 8.0 * arg1 * arg2 * arg2 * arg2 * (term3 - (2.0 / 3.0) * term6 * term3);
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
      return 8.0 * arg1 * arg2 * arg2 * arg2 * (term3 - (2.0 / 3.0) * term9) +
             fh * 96.0 * arg1 * rc * (-7.0 * term12 + 2.0 * term6);
    }
    case VDWParameters::Type::BuckingHam:
    {
      double A = parameters.parameters.x;
      double b = parameters.parameters.y;
      double C6 = parameters.parameters.z;
      double br = b * rc;
      return (2.0 * C6) / (rc * rc * rc) -
             (A * (6.0 + br * (6.0 + br * (3.0 + br)))) / (b * b * b * std::exp(br));
    }
    case VDWParameters::Type::Morse:
    {
      double D = parameters.parameters.x;
      double a = parameters.parameters.y;
      double r0 = parameters.parameters.z;
      double ar = a * rc;
      return (D * std::exp(a * (r0 - 2.0 * rc)) *
              (8.0 * std::exp(ar) * (6.0 + ar * (6.0 + ar * (3.0 + ar))) -
               std::exp(a * r0) * (3.0 + 2.0 * ar * (3.0 + ar * (3.0 + 2.0 * ar))))) /
             (4.0 * a * a * a);
    }
    case VDWParameters::Type::MM3:
    {
      double epsilon = parameters.parameters.x;
      double sigma = parameters.parameters.y;
      double sigma2 = sigma * sigma;
      double sigma3 = sigma2 * sigma;
      double sigma6 = sigma3 * sigma3;
      double rc3 = rc * rc * rc;
      return (2.0 * epsilon * sigma6 * 2.25) / rc3 -
             (1.84e5 * epsilon * (sigma3 + 12.0 * sigma2 * rc + 72.0 * sigma * rc * rc + 288.0 * rc3)) /
                 (288.0 * std::exp((12.0 * rc) / sigma));
    }
    case VDWParameters::Type::BornHugginsMeyer:
    {
      double A = parameters.parameters.x;
      double b = parameters.parameters.y;
      double sigma = parameters.parameters.z;
      double C6 = parameters.parameters.w;
      double C8 = parameters.parameters2.x;
      double br = b * rc;
      double rc3 = rc * rc * rc;
      double rc5 = rc3 * rc * rc;
      return (8.0 * C8) / (5.0 * rc5) + (2.0 * C6) / rc3 -
             (A * std::exp(b * (sigma - rc)) * (6.0 + br * (6.0 + br * (3.0 + br)))) / (b * b * b);
    }
    case VDWParameters::Type::Potential12_6:
    {
      double p0 = parameters.parameters.x;
      double p1 = parameters.parameters.y;
      double rc3 = rc * rc * rc;
      double rc9 = rc3 * rc3 * rc3;
      return (-4.0 * p0) / (3.0 * rc9) + (2.0 * p1) / rc3;
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
      return 2.0 * p1 / rc3 - 1.5 * p0 / rc6;
    }
    case VDWParameters::Type::CFFEpsilonSigma:
    {
      double epsilon = parameters.parameters.x;
      double sigma = parameters.parameters.y;
      double sigma3 = sigma * sigma * sigma;
      double term2 = 3.0 * epsilon * sigma3 * sigma3;
      double rc3 = rc * rc * rc;
      return 2.0 * term2 / rc3 - term2 * sigma3 / (rc3 * rc3);
    }
    case VDWParameters::Type::MatsuokaClementiYoshimine:
    {
      double A = parameters.parameters.x;
      double a = parameters.parameters.y;
      double B = parameters.parameters.z;
      double b = parameters.parameters.w;
      double ar = a * rc;
      double br = b * rc;
      double a3 = a * a * a;
      double b3 = b * b * b;
      return (std::exp((-a - b) * rc) *
              (-(A * b3 * std::exp(br) * (6.0 + ar * (6.0 + ar * (3.0 + ar)))) -
               a3 * B * std::exp(ar) * (6.0 + br * (6.0 + br * (3.0 + br))))) /
             (a3 * b3);
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
      return -(A * (6.0 + br * (6.0 + br * (3.0 + br))) * std::exp(-br)) / (b * b * b) + (4.0 * C4) / rc +
             (2.0 * C6) / rc3 + (8.0 * C8) / (5.0 * rc5) + (10.0 * C10) / (7.0 * rc7);
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
      double term1 = -A * (6.0 + br * (6.0 + br * (3.0 + br))) * expTerm / (b * b * b);
      double term2 = (C6 * (1440.0 - 1440.0 * std::exp(br) +
                            br * (1440.0 + br * (720.0 + br * (234.0 + br * (54.0 + br * (9.0 + br))))))) /
                     (720.0 * rc3);
      double term3 =
          (C8 *
           (64512.0 - 64512.0 * std::exp(br) +
            br * (64512.0 +
                  br * (32256.0 +
                        br * (10752.0 + br * (2688.0 + br * (534.0 + br * (86.0 + br * (11.0 + br))))))))) /
          (40320.0 * rc5);
      double term4 =
          (C10 *
           (5184000.0 - 5184000.0 * std::exp(br) +
            br * (5184000.0 +
                  br * (2592000.0 +
                        br * (864000.0 +
                              br * (216000.0 +
                                    br * (43200.0 +
                                          br * (7200.0 + br * (1026.0 + br * (126.0 + br * (13.0 + br))))))))))) /
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
      return (4.0 * C12) / (3.0 * rc9) + (2.0 * C6) / rc3 + (4.0 * C4) / rc -
             std::exp(-br) * (A * (6.0 + br * (6.0 + br * (3.0 + br)))) / (b * b * b);
    }
    case VDWParameters::Type::Mie:
    {
      double p0 = parameters.parameters.x;
      double n = parameters.parameters.y;
      double p2 = parameters.parameters.z;
      double m = parameters.parameters.w;
      return (p2 * m * std::pow(rc, 3.0 - m)) / (m - 3.0) - (p0 * n * std::pow(rc, 3.0 - n)) / (n - 3.0);
    }
    case VDWParameters::Type::LennardJonesShiftedForce:
    case VDWParameters::Type::LennardJonesSecondOrderTaylorShifted:
    case VDWParameters::Type::WeeksChandlerAndersen:
      return 0.0;
    default:
      return 0.0;
  }
};
}  // namespace Potentials
