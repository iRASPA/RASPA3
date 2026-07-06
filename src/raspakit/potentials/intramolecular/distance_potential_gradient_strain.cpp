module;

module distance_potential_gradient_strain;

import std;

import double3;
import double3x3;
import bond_potential;

std::tuple<double, std::array<double3, 2>, double3x3> Potentials::Internal::distancePotentialEnergyGradientStrain(
    BondType type, const std::array<double, maximumNumberOfBondParameters> &parameters, const double3 &posA,
    const double3 &posB, bool zero_gradient)
{
  double temp, temp2, exp_term, r1;
  double U{}, DF{};
  double3 du_dr, du_da, du_db;
  double3x3 strain_derivative{};

  double3 dr = posA - posB;
  double rr = double3::dot(dr, dr);
  double r = std::sqrt(rr);

  if (!zero_gradient && r > 0.0)
  {
    switch (type)
    {
      case BondType::None:
      case BondType::Fixed:
        U = 0.0;
        DF = 0.0;
        break;
      case BondType::Harmonic:
        U = 0.5 * parameters[0] * (r - parameters[1]) * (r - parameters[1]);
        DF = parameters[0] * (r - parameters[1]) / r;
        break;
      case BondType::CoreShellSpring:
        U = 0.5 * parameters[0] * r * r;
        DF = parameters[0];
        break;
      case BondType::Morse:
        temp = std::exp(parameters[1] * (parameters[2] - r));
        U = parameters[0] * ((1.0 - temp) * (1.0 - temp) - 1.0);
        DF = 2.0 * parameters[0] * parameters[1] * (1.0 - temp) * temp / r;
        break;
      case BondType::LJ_12_6:
        temp = 1.0 / (rr * rr * rr);
        U = parameters[0] * temp * temp - parameters[1] * temp;
        DF = 6.0 * (parameters[1] * temp - 2.0 * parameters[0] * temp * temp) / rr;
        break;
      case BondType::LennardJones:
        temp = (parameters[1] / rr) * (parameters[1] / rr) * (parameters[1] / rr);
        U = 4.0 * parameters[0] * (temp * (temp - 1.0));
        DF = 24.0 * parameters[0] * (temp * (1.0 - 2.0 * temp)) / rr;
        break;
      case BondType::Buckingham:
        temp = parameters[2] / (rr * rr * rr);
        exp_term = parameters[0] * std::exp(-parameters[1] * r);
        U = -temp + exp_term;
        DF = (6.0 / rr) * temp - parameters[1] * exp_term / r;
        break;
      case BondType::RestrainedHarmonic:
        r1 = r - parameters[1];
        U = 0.5 * parameters[0] * std::pow(std::min(std::fabs(r1), parameters[2]), 2) +
            parameters[0] * parameters[2] * std::max(std::fabs(r1) - parameters[2], 0.0);
        DF = -parameters[0] * (std::copysign(std::min(std::fabs(r1), parameters[2]), r1)) / r;
        break;
      case BondType::Quartic:
        temp = r - parameters[1];
        temp2 = temp * temp;
        U = 0.5 * parameters[0] * temp2 + (1.0 / 3.0) * parameters[2] * temp * temp2 +
            0.25 * parameters[3] * temp2 * temp2;
        DF = temp * (parameters[0] + parameters[2] * temp + parameters[3] * temp2) / r;
        break;
      case BondType::CFF_Quartic:
        temp = r - parameters[1];
        temp2 = temp * temp;
        U = parameters[0] * temp2 + parameters[2] * temp * temp2 + parameters[3] * temp2 * temp2;
        DF = temp * (2.0 * parameters[0] + 3.0 * parameters[2] * temp + 4.0 * parameters[3] * temp2) / r;
        break;
      case BondType::MM3:
        temp = r - parameters[1];
        temp2 = temp * temp;
        U = parameters[0] * temp2 * (1.0 - 2.55 * temp + (7.0 / 12.0) * 2.55 * 2.55 * temp2);
        DF = parameters[0] * (2.0 + 2.55 * (4.0 * 2.55 * (7.0 / 12.0) * temp - 3.0) * temp) * temp / r;
        break;
      default:
        std::unreachable();
    }
  }
  else
  {
    switch (type)
    {
      case BondType::None:
      case BondType::Fixed:
        U = 0.0;
        break;
      case BondType::Harmonic:
        U = 0.5 * parameters[0] * (r - parameters[1]) * (r - parameters[1]);
        break;
      case BondType::CoreShellSpring:
        U = 0.5 * parameters[0] * r * r;
        break;
      case BondType::Morse:
        temp = std::exp(parameters[1] * (parameters[2] - r));
        U = parameters[0] * ((1.0 - temp) * (1.0 - temp) - 1.0);
        break;
      case BondType::LJ_12_6:
        temp = 1.0 / (rr * rr * rr);
        U = parameters[0] * temp * temp - parameters[1] * temp;
        break;
      case BondType::LennardJones:
        temp = (parameters[1] / rr) * (parameters[1] / rr) * (parameters[1] / rr);
        U = 4.0 * parameters[0] * (temp * (temp - 1.0));
        break;
      case BondType::Buckingham:
        temp = parameters[2] / (rr * rr * rr);
        U = parameters[0] * std::exp(-parameters[1] * r) - temp;
        break;
      case BondType::RestrainedHarmonic:
        r1 = r - parameters[1];
        U = 0.5 * parameters[0] * std::pow(std::min(std::fabs(r1), parameters[2]), 2) +
            parameters[0] * parameters[2] * std::max(std::fabs(r1) - parameters[2], 0.0);
        break;
      case BondType::Quartic:
        temp = r - parameters[1];
        temp2 = temp * temp;
        U = 0.5 * parameters[0] * temp2 + (1.0 / 3.0) * parameters[2] * temp * temp2 +
            0.25 * parameters[3] * temp2 * temp2;
        break;
      case BondType::CFF_Quartic:
        temp = r - parameters[1];
        temp2 = temp * temp;
        U = parameters[0] * temp2 + parameters[2] * temp * temp2 + parameters[3] * temp2 * temp2;
        break;
      case BondType::MM3:
        temp = r - parameters[1];
        temp2 = temp * temp;
        U = parameters[0] * temp2 * (1.0 - 2.55 * temp + (7.0 / 12.0) * 2.55 * 2.55 * temp2);
        break;
      default:
        std::unreachable();
    }
    DF = 0.0;
  }

  du_dr = DF * dr;
  du_da = du_dr;
  du_db = -du_dr;

  strain_derivative.ax = dr.x * du_dr.x;
  strain_derivative.bx = dr.y * du_dr.x;
  strain_derivative.cx = dr.z * du_dr.x;
  strain_derivative.ay = dr.x * du_dr.y;
  strain_derivative.by = dr.y * du_dr.y;
  strain_derivative.cy = dr.z * du_dr.y;
  strain_derivative.az = dr.x * du_dr.z;
  strain_derivative.bz = dr.y * du_dr.z;
  strain_derivative.cz = dr.z * du_dr.z;

  return {U, {du_da, du_db}, strain_derivative};
}
