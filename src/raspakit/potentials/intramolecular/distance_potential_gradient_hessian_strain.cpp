module;

module distance_potential_gradient_hessian_strain;

import std;

import double3;
import double3x3;
import bond_potential;

namespace
{
struct DistanceEnergyDerivatives
{
  double energy{};
  double firstDerivative{};   // dU/dr
  double secondDerivative{};  // d^2U/dr^2
};

DistanceEnergyDerivatives distanceEnergyDerivatives(BondType type,
                                                    const std::array<double, maximumNumberOfBondParameters> &parameters,
                                                    double r, double rr)
{
  DistanceEnergyDerivatives result{};
  double temp{};
  double temp2{};
  double exp_term{};

  switch (type)
  {
    case BondType::None:
    case BondType::Fixed:
      break;
    case BondType::Harmonic:
      result.energy = 0.5 * parameters[0] * (r - parameters[1]) * (r - parameters[1]);
      result.firstDerivative = parameters[0] * (r - parameters[1]);
      result.secondDerivative = parameters[0];
      break;
    case BondType::CoreShellSpring:
      result.energy = 0.5 * parameters[0] * r * r;
      result.firstDerivative = parameters[0] * r;
      result.secondDerivative = parameters[0];
      break;
    case BondType::Morse:
      temp = std::exp(parameters[1] * (parameters[2] - r));
      result.energy = parameters[0] * ((1.0 - temp) * (1.0 - temp) - 1.0);
      result.firstDerivative = 2.0 * parameters[0] * parameters[1] * (1.0 - temp) * temp;
      result.secondDerivative =
          2.0 * parameters[0] * parameters[1] * parameters[1] * temp * (2.0 * temp - 1.0);
      break;
    case BondType::LJ_12_6:
      temp = 1.0 / (rr * rr * rr);
      result.energy = parameters[0] * temp * temp - parameters[1] * temp;
      result.firstDerivative = 6.0 * (parameters[1] * temp - 2.0 * parameters[0] * temp * temp) / r;
      result.secondDerivative = 24.0 * (7.0 * parameters[0] * temp * temp - 2.0 * parameters[1] * temp) / rr;
      break;
    case BondType::LennardJones:
      temp = (parameters[1] / rr) * (parameters[1] / rr) * (parameters[1] / rr);
      result.energy = 4.0 * parameters[0] * (temp * (temp - 1.0));
      result.firstDerivative = 24.0 * parameters[0] * (temp * (1.0 - 2.0 * temp)) / r;
      result.secondDerivative = 96.0 * parameters[0] * temp * (7.0 * temp - 2.0) / rr;
      break;
    case BondType::Buckingham:
      temp = parameters[2] / (rr * rr * rr);
      exp_term = parameters[0] * std::exp(-parameters[1] * r);
      result.energy = -temp + exp_term;
      result.firstDerivative = 6.0 * temp / r + parameters[1] * exp_term;
      result.secondDerivative = (-48.0 * temp / rr + parameters[1] * parameters[1] * exp_term) / r;
      break;
    case BondType::RestrainedHarmonic:
    {
      const double r1 = r - parameters[1];
      const double abs_r1 = std::fabs(r1);
      if (abs_r1 <= parameters[2])
      {
        result.energy = 0.5 * parameters[0] * r1 * r1;
        result.firstDerivative = parameters[0] * r1;
        result.secondDerivative = parameters[0];
      }
      else
      {
        result.energy = 0.5 * parameters[0] * parameters[2] * parameters[2] +
                        parameters[0] * parameters[2] * (abs_r1 - parameters[2]);
        result.firstDerivative = parameters[0] * parameters[2] * std::copysign(1.0, r1);
        result.secondDerivative = 0.0;
      }
      break;
    }
    case BondType::Quartic:
      temp = r - parameters[1];
      temp2 = temp * temp;
      result.energy = 0.5 * parameters[0] * temp2 + (1.0 / 3.0) * parameters[2] * temp * temp2 +
                      0.25 * parameters[3] * temp2 * temp2;
      result.firstDerivative = temp * (parameters[0] + parameters[2] * temp + parameters[3] * temp2);
      result.secondDerivative = parameters[0] + 2.0 * parameters[2] * temp + 3.0 * parameters[3] * temp2;
      break;
    case BondType::CFF_Quartic:
      temp = r - parameters[1];
      temp2 = temp * temp;
      result.energy = parameters[0] * temp2 + parameters[2] * temp * temp2 + parameters[3] * temp2 * temp2;
      result.firstDerivative = temp * (2.0 * parameters[0] + 3.0 * parameters[2] * temp + 4.0 * parameters[3] * temp2);
      result.secondDerivative = 2.0 * parameters[0] + 6.0 * parameters[2] * temp + 12.0 * parameters[3] * temp2;
      break;
    case BondType::MM3:
      temp = r - parameters[1];
      temp2 = temp * temp;
      result.energy = parameters[0] * temp2 * (1.0 - 2.55 * temp + (7.0 / 12.0) * 2.55 * 2.55 * temp2);
      result.firstDerivative =
          parameters[0] * (2.0 + 2.55 * (4.0 * 2.55 * (7.0 / 12.0) * temp - 3.0) * temp) * temp;
      result.secondDerivative =
          parameters[0] *
          (2.0 + 2.55 * (3.0 * 4.0 * 2.55 * (7.0 / 12.0) * temp - 3.0) * temp + 2.55 * (4.0 * 2.55 * (7.0 / 12.0) * temp - 3.0) * temp);
      break;
    default:
      std::unreachable();
  }

  return result;
}
}  // namespace

std::tuple<double, std::array<double3, 2>, double3x3, double, double>
Potentials::Internal::distancePotentialEnergyGradientHessianStrain(
    BondType type, const std::array<double, maximumNumberOfBondParameters> &parameters, const double3 &posA,
    const double3 &posB, bool zero_gradient)
{
  const double3 dr = posA - posB;
  const double rr = double3::dot(dr, dr);
  const double r = std::sqrt(rr);

  DistanceEnergyDerivatives derivatives{};
  if (r > 0.0)
  {
    derivatives = distanceEnergyDerivatives(type, parameters, r, rr);
  }
  else
  {
    auto energy_only = distanceEnergyDerivatives(type, parameters, r, rr);
    derivatives.energy = energy_only.energy;
  }

  double f1{};
  double f2{};
  double3 du_dr{};
  if (!zero_gradient && r > 0.0)
  {
    f1 = derivatives.firstDerivative / r;
    f2 = derivatives.secondDerivative / rr - derivatives.firstDerivative / (r * rr);
    du_dr = f1 * dr;
  }

  const double3 du_da = du_dr;
  const double3 du_db = -du_dr;

  double3x3 strain_derivative{};
  strain_derivative.ax = dr.x * du_dr.x;
  strain_derivative.bx = dr.y * du_dr.x;
  strain_derivative.cx = dr.z * du_dr.x;
  strain_derivative.ay = dr.x * du_dr.y;
  strain_derivative.by = dr.y * du_dr.y;
  strain_derivative.cy = dr.z * du_dr.y;
  strain_derivative.az = dr.x * du_dr.z;
  strain_derivative.bz = dr.y * du_dr.z;
  strain_derivative.cz = dr.z * du_dr.z;

  return {derivatives.energy, {du_da, du_db}, strain_derivative, f1, f2};
}
