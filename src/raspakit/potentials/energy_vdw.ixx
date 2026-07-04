module;

export module potential_energy_vdw;

import std;

import vdwparameters;
import forcefield;
import energy_factor;

import double4;

export namespace Potentials
{
/**
 * \brief Calculates the van der Waals potential energy between two atoms.
 *
 * This function computes the van der Waals (VDW) potential energy based on the specified
 * force field parameters, group identifiers, scaling factors, distance between atoms,
 * and their respective types. It supports the potential types ported from RASPA2:
 * Lennard-Jones, Buckingham, Morse, Feynman-Hibbs Lennard-Jones, MM3, Born-Huggins-Meyer,
 * shifted-force Lennard-Jones, 12-6, 12-6-2-0, CFF 9-6, CFF eps-sigma,
 * Matsuoka-Clementi-Yoshimine, generic, Pellenq-Nicholson, hydrated-ion-water, Mie,
 * and Weeks-Chandler-Andersen.
 *
 * The scaling is linear: it first activates van der Waals interactions from 0 to 0.5,
 * then activates electrostatic interactions from 0.5 to 1.0.
 *
 * All potentials use the same soft-core lambda-scaling as the Lennard-Jones potential:
 * the potential is evaluated at the softened distance
 *   r_soft^6 = r^6 + 0.5 (1 - lambda)^2 w^6
 * and multiplied by lambda, where w is the soft-core reference diameter of the pair
 * potential (precomputed in VDWParameters::computeDerivedParameters(), stored in
 * parameters2.w). At lambda = 1 the original potential is recovered exactly; for
 * lambda < 1 the energy remains finite in the limit r -> 0. Potentials that are
 * already finite at r = 0 (Morse, Matsuoka-Clementi-Yoshimine) use w = 0, i.e. plain
 * linear scaling.
 *
 * The lambda-derivative of a soft-core potential U(lambda) = lambda * u(r_soft(lambda)) is
 *   dU/dlambda = u(r_soft) + lambda * (-u'(r_soft)) * r_soft * (1 - lambda) * w^6 / (6 r_soft^6)
 * where -u'(r_soft) * r_soft is computed analytically per potential ('deriv' below).
 *
 * \param forcefield The force field parameters used for the calculation.
 * \param groupIdA Boolean indicating if the first atom is part of a specific group.
 * \param groupIdB Boolean indicating if the second atom is part of a specific group.
 * \param scalingA Scaling factor for the first atom's interactions.
 * \param scalingB Scaling factor for the second atom's interactions.
 * \param rr The squared distance between the two atoms.
 * \param typeA The type identifier for the first atom.
 * \param typeB The type identifier for the second atom.
 * \return An EnergyFactor object containing the calculated potential energy and lambda derivative.
 */

// Slow path with the potential types ported from RASPA2. Deliberately not inlined: keeping
// these out of the always-inlined hot function avoids code bloat at every call site and
// preserves the performance of the Lennard-Jones code path. The 'preserve_most' calling
// convention makes this rarely-taken call save the registers itself, so the caller's hot
// loop does not have to spill around it.
[[clang::noinline]] [[clang::preserve_most]] inline EnergyFactor potentialVDWEnergyRare(
    const ForceField& forcefield, const bool groupIdA, const bool groupIdB, const double scalingA,
    const double scalingB, const double rr, const std::size_t typeA, const std::size_t typeB,
    VDWParameters::Type potentialType)
{
  double scaling = scalingA * scalingB;
  switch (potentialType)
  {
    case VDWParameters::Type::BuckingHam:
    {
      // p_0*exp(-p_1*r) - p_2/r^6
      const VDWParameters& p = forcefield(typeA, typeB);
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double rs = std::sqrt(std::cbrt(x6));
      double repulsion = p.parameters.x * std::exp(-p.parameters.y * rs);
      double dispersion = p.parameters.z / x6;
      double term = repulsion - dispersion - p.shift;
      double deriv = p.parameters.y * rs * repulsion - 6.0 * dispersion;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::Morse:
    {
      // p_0*{(1.0-exp[-p_1*(r-p_2)])^2 - 1.0}; finite at r=0, plain linear lambda-scaling
      const VDWParameters& p = forcefield(typeA, typeB);
      double r = std::sqrt(rr);
      double expTerm = std::exp(-p.parameters.y * (r - p.parameters.z));
      double term = p.parameters.x * ((1.0 - expTerm) * (1.0 - expTerm) - 1.0) - p.shift;
      return EnergyFactor(scaling * term,
                          (groupIdA ? scalingB * term : 0.0) + (groupIdB ? scalingA * term : 0.0));
    }
    case VDWParameters::Type::FeynmannHibbs:
    {
      // 4*p_0*((p_1/r)^12-(p_1/r)^6) + (hbar^2/(24 mu kB T))*4*p_0*(132*(p_1/r)^12-30*(p_1/r)^6)/r^2
      const VDWParameters& p = forcefield(typeA, typeB);
      double eps4 = 4.0 * p.parameters.x;
      double sigma2 = p.parameters.y * p.parameters.y;
      double sigma6 = sigma2 * sigma2 * sigma2;
      double fh = p.parameters2.x;
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double rs2 = std::cbrt(x6);
      double u6 = sigma6 / x6;
      double u12 = u6 * u6;
      double quantumCorrection = eps4 * fh * (132.0 * u12 - 30.0 * u6) / rs2;
      double term = eps4 * (u12 - u6) + quantumCorrection - p.shift;
      double deriv =
          eps4 * (12.0 * u12 - 6.0 * u6) + eps4 * fh * (14.0 * 132.0 * u12 - 8.0 * 30.0 * u6) / rs2;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::MM3:
    {
      // eps*[1.84e5*exp(-12/P) - 2.25*P^6] with P = sigma/r; eps*192.27*P^2 for P > 3.02
      const VDWParameters& p = forcefield(typeA, typeB);
      double epsilon = p.parameters.x;
      double sigma = p.parameters.y;
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double rs = std::sqrt(std::cbrt(x6));
      double P = sigma / rs;
      double term;
      double deriv;
      if (P > 3.02)
      {
        double value = epsilon * 192.270 * P * P;
        term = value - p.shift;
        deriv = 2.0 * value;
      }
      else
      {
        double P2 = P * P;
        double P6 = P2 * P2 * P2;
        double repulsion = epsilon * 1.84e5 * std::exp(-12.0 * rs / sigma);
        double dispersion = epsilon * 2.25 * P6;
        term = repulsion - dispersion - p.shift;
        deriv = (12.0 * rs / sigma) * repulsion - 6.0 * dispersion;
      }
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::BornHugginsMeyer:
    {
      // p_0*exp(p_1*(p_2-r)) - p_3/r^6 - p_4/r^8
      const VDWParameters& p = forcefield(typeA, typeB);
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double rs2 = std::cbrt(x6);
      double rs = std::sqrt(rs2);
      double t6 = 1.0 / x6;
      double t8 = t6 / rs2;
      double repulsion = p.parameters.x * std::exp(p.parameters.y * (p.parameters.z - rs));
      double dispersion6 = p.parameters.w * t6;
      double dispersion8 = p.parameters2.x * t8;
      double term = repulsion - dispersion6 - dispersion8 - p.shift;
      double deriv = p.parameters.y * rs * repulsion - 6.0 * dispersion6 - 8.0 * dispersion8;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::LennardJonesShiftedForce:
    {
      // 4*p_0*{[(p_1/r)^12-(p_1/r)^6] - [(p_1/rc)^12-(p_1/rc)^6] + [12*(p_1/rc)^12-6*(p_1/rc)^6]*(r-rc)/rc}
      const VDWParameters& p = forcefield(typeA, typeB);
      double eps4 = 4.0 * p.parameters.x;
      double sigma2 = p.parameters.y * p.parameters.y;
      double sigma6 = sigma2 * sigma2 * sigma2;
      double c6 = p.parameters2.x;   // (sigma/rc)^6
      double rc = p.parameters2.y;
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double rs = std::sqrt(std::cbrt(x6));
      double u6 = sigma6 / x6;
      double u12 = u6 * u6;
      double linearCoefficient = 12.0 * c6 * c6 - 6.0 * c6;
      double term = eps4 * (u12 - u6 - c6 * (c6 - 1.0) + linearCoefficient * (rs - rc) / rc);
      double deriv = eps4 * (12.0 * u12 - 6.0 * u6) - eps4 * linearCoefficient * rs / rc;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::Potential12_6:
    {
      // p_0/r^12 - p_1/r^6
      const VDWParameters& p = forcefield(typeA, typeB);
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double t6 = 1.0 / x6;
      double repulsion = p.parameters.x * t6 * t6;
      double dispersion = p.parameters.y * t6;
      double term = repulsion - dispersion - p.shift;
      double deriv = 12.0 * repulsion - 6.0 * dispersion;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::Potential12_6_2_0:
    {
      // p_0/r^12 + p_1/r^6 + p_2/r^2 + p_3
      const VDWParameters& p = forcefield(typeA, typeB);
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double rs2 = std::cbrt(x6);
      double t6 = 1.0 / x6;
      double term12 = p.parameters.x * t6 * t6;
      double term6 = p.parameters.y * t6;
      double term2 = p.parameters.z / rs2;
      double term = term12 + term6 + term2 + p.parameters.w - p.shift;
      double deriv = 12.0 * term12 + 6.0 * term6 + 2.0 * term2;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::CFF9_6:
    {
      // p_0/r^9 - p_1/r^6
      const VDWParameters& p = forcefield(typeA, typeB);
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double t6 = 1.0 / x6;
      double t9 = t6 * std::sqrt(t6);
      double repulsion = p.parameters.x * t9;
      double dispersion = p.parameters.y * t6;
      double term = repulsion - dispersion - p.shift;
      double deriv = 9.0 * repulsion - 6.0 * dispersion;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::CFFEpsilonSigma:
    {
      // p_0*[2*(p_1/r)^9 - 3*(p_1/r)^6]
      const VDWParameters& p = forcefield(typeA, typeB);
      double epsilon = p.parameters.x;
      double sigma2 = p.parameters.y * p.parameters.y;
      double sigma6 = sigma2 * sigma2 * sigma2;
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double u6 = sigma6 / x6;
      double u9 = u6 * std::sqrt(u6);
      double repulsion = 2.0 * epsilon * u9;
      double dispersion = 3.0 * epsilon * u6;
      double term = repulsion - dispersion - p.shift;
      double deriv = 9.0 * repulsion - 6.0 * dispersion;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::MatsuokaClementiYoshimine:
    {
      // p_0*exp(-p_1*r) + p_2*exp(-p_3*r); finite at r=0, plain linear lambda-scaling
      const VDWParameters& p = forcefield(typeA, typeB);
      double r = std::sqrt(rr);
      double term = p.parameters.x * std::exp(-p.parameters.y * r) +
                    p.parameters.z * std::exp(-p.parameters.w * r) - p.shift;
      return EnergyFactor(scaling * term,
                          (groupIdA ? scalingB * term : 0.0) + (groupIdB ? scalingA * term : 0.0));
    }
    case VDWParameters::Type::Generic:
    {
      // p_0*exp(-p_1*r) - p_2/r^4 - p_3/r^6 - p_4/r^8 - p_5/r^10
      const VDWParameters& p = forcefield(typeA, typeB);
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double rs2 = std::cbrt(x6);
      double rs = std::sqrt(rs2);
      double t2 = 1.0 / rs2;
      double t4 = t2 * t2;
      double t6 = t4 * t2;
      double t8 = t6 * t2;
      double t10 = t8 * t2;
      double repulsion = p.parameters.x * std::exp(-p.parameters.y * rs);
      double dispersion4 = p.parameters.z * t4;
      double dispersion6 = p.parameters.w * t6;
      double dispersion8 = p.parameters2.x * t8;
      double dispersion10 = p.parameters2.y * t10;
      double term = repulsion - dispersion4 - dispersion6 - dispersion8 - dispersion10 - p.shift;
      double deriv = p.parameters.y * rs * repulsion - 4.0 * dispersion4 - 6.0 * dispersion6 -
                     8.0 * dispersion8 - 10.0 * dispersion10;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::PellenqNicholson:
    {
      // p_0*exp(-p_1*r) - f_6*p_2/r^6 - f_8*p_3/r^8 - f_10*p_4/r^10 with Tang-Toennies damping f_n
      const VDWParameters& p = forcefield(typeA, typeB);
      double b = p.parameters.y;
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double rs2 = std::cbrt(x6);
      double rs = std::sqrt(rs2);
      double f6, f8, f10;
      computeTangToenniesDampingCoefficients(rs, b, f6, f8, f10);
      double t2 = 1.0 / rs2;
      double t6 = t2 * t2 * t2;
      double t8 = t6 * t2;
      double t10 = t8 * t2;
      double expTerm = std::exp(-b * rs);
      double repulsion = p.parameters.x * expTerm;
      double dispersion6 = f6 * p.parameters.z * t6;
      double dispersion8 = f8 * p.parameters.w * t8;
      double dispersion10 = f10 * p.parameters2.x * t10;
      double term = repulsion - dispersion6 - dispersion8 - dispersion10 - p.shift;
      // damping derivatives: f_n'(r) = b (b r)^n exp(-b r)/n!
      double br = b * rs;
      double br2 = br * br;
      double br6 = br2 * br2 * br2;
      double br8 = br6 * br2;
      double br10 = br8 * br2;
      double dampingDeriv = b * rs * expTerm *
                            (p.parameters.z * t6 * br6 / 720.0 + p.parameters.w * t8 * br8 / 40320.0 +
                             p.parameters2.x * t10 * br10 / 3628800.0);
      double deriv = b * rs * repulsion - 6.0 * dispersion6 - 8.0 * dispersion8 - 10.0 * dispersion10 +
                     dampingDeriv;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::HydratedIonWater:
    {
      // p_0*exp(-p_1*r) - p_2/r^4 - p_3/r^6 - p_4/r^12
      const VDWParameters& p = forcefield(typeA, typeB);
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double rs2 = std::cbrt(x6);
      double rs = std::sqrt(rs2);
      double t2 = 1.0 / rs2;
      double t4 = t2 * t2;
      double t6 = t4 * t2;
      double t12 = t6 * t6;
      double repulsion = p.parameters.x * std::exp(-p.parameters.y * rs);
      double dispersion4 = p.parameters.z * t4;
      double dispersion6 = p.parameters.w * t6;
      double repulsion12 = p.parameters2.x * t12;
      double term = repulsion - dispersion4 - dispersion6 - repulsion12 - p.shift;
      double deriv =
          p.parameters.y * rs * repulsion - 4.0 * dispersion4 - 6.0 * dispersion6 - 12.0 * repulsion12;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::Mie:
    {
      // p_0/r^p_1 - p_2/r^p_3
      const VDWParameters& p = forcefield(typeA, typeB);
      double inv = 1.0 - scaling;
      double w6 = p.parameters2.w;
      double x6 = rr * rr * rr + 0.5 * inv * inv * w6;
      double rs = std::sqrt(std::cbrt(x6));
      double repulsion = p.parameters.x / std::pow(rs, p.parameters.y);
      double attraction = p.parameters.z / std::pow(rs, p.parameters.w);
      double term = repulsion - attraction - p.shift;
      double deriv = p.parameters.y * repulsion - p.parameters.w * attraction;
      double dlambda_term = scaling * inv * w6 * deriv / (6.0 * x6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::WeeksChandlerAndersen:
    {
      // 4*p_0*((p_1/r)^12-(p_1/r)^6)+p_0 if r < 2^(1/6)*p_1, otherwise 0.
      // The cutoff constant is evaluated at the same soft-core lambda so that the
      // potential remains continuous at the cutoff for every value of lambda.
      const VDWParameters& p = forcefield(typeA, typeB);
      double sigma2 = p.parameters.y * p.parameters.y;
      if (rr > std::cbrt(2.0) * sigma2)
      {
        return EnergyFactor(0.0, 0.0);
      }
      double eps4 = 4.0 * p.parameters.x;
      double sigma6 = sigma2 * sigma2 * sigma2;
      double inv = 1.0 - scaling;
      double x6 = rr * rr * rr + 0.5 * inv * inv * sigma6;
      double u6 = sigma6 / x6;
      double u12 = u6 * u6;
      double gcut = 1.0 / (2.0 + 0.5 * inv * inv);
      double term = eps4 * (u12 - u6) - eps4 * (gcut * gcut - gcut);
      double dlambda_term = scaling * (inv * sigma6 * eps4 * (12.0 * u12 - 6.0 * u6) / (6.0 * x6) -
                                       eps4 * inv * (2.0 * gcut - 1.0) * gcut * gcut);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::RepulsiveHarmonic:
    {
      double r = std::sqrt(rr);
      double arg1 = forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y;
      double temp = (1.0 - r / arg2);
      return EnergyFactor(r >= arg2 ? 0.0 : 0.5 * arg1 * temp * temp, 0.0);
    }
    default:
      return EnergyFactor(0.0, 0.0);
  }

  return EnergyFactor(0.0, 0.0);
}

[[clang::always_inline]] inline EnergyFactor potentialVDWEnergy(const ForceField& forcefield, const bool& groupIdA,
                                                                const bool& groupIdB, const double& scalingA,
                                                                const double& scalingB, const double& rr,
                                                                const std::size_t& typeA, const std::size_t& typeB)
{
  VDWParameters::Type potentialType = forcefield(typeA, typeB).type;

  double scaling = scalingA * scalingB;

  // The Lennard-Jones potential is dispatched with a single compare-and-branch; all other
  // potential types are handled in the non-inlined function above so that the hot code path
  // stays identical to the one before the RASPA2 potentials were added.
  if (potentialType == VDWParameters::Type::LennardJones) [[likely]]
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
    return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                            (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
  }

  return potentialVDWEnergyRare(forcefield, groupIdA, groupIdB, scalingA, scalingB, rr, typeA, typeB, potentialType);
};
}  // namespace Potentials
