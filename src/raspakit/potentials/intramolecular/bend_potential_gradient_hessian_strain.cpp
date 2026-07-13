module;

module bend_potential_gradient_hessian_strain;

import std;

import units;

namespace
{
double cube(double x) { return x * x * x; }
}  // namespace

std::tuple<double, std::array<double3, 3>, double3x3, Minimization::BendHessianGeometry>
Potentials::Internal::bendPotentialEnergyGradientHessianStrain(
    BendType type, const std::array<double, maximumNumberOfBendParameters> &parameters, const double3 &posA,
    const double3 &posB, const double3 &posC)
{
  Minimization::BendHessianGeometry geometry{};

  double3 dr_ab = posA - posB;
  const double r_ab = std::sqrt(double3::dot(dr_ab, dr_ab));
  dr_ab /= r_ab;

  double3 dr_cb = posC - posB;
  const double r_cb = std::sqrt(double3::dot(dr_cb, dr_cb));
  dr_cb /= r_cb;

  double cos_theta = double3::dot(dr_ab, dr_cb);
  cos_theta = std::clamp(cos_theta, -1.0, 1.0);

  const double theta = std::acos(cos_theta);
  const double DTDX = -1.0 / std::sqrt(1.0 - cos_theta * cos_theta);

  double U{};
  double DF{};
  double DDF{};
  double temp{};
  double temp2{};

  switch (type)
  {
    case BendType::Fixed:
    case BendType::Rigid:
      U = 0.0;
      DF = 0.0;
      DDF = 0.0;
      break;
    case BendType::Harmonic:
    case BendType::CoreShell:
      temp = theta - parameters[1];
      U = 0.5 * parameters[0] * temp * temp;
      DF = parameters[0] * temp * DTDX;
      DDF = parameters[0] * DTDX * DTDX + parameters[0] * temp * cos_theta * cube(DTDX);
      break;
    case BendType::Quartic:
      temp = theta - parameters[1];
      temp2 = temp * temp;
      U = 0.5 * parameters[0] * temp2 + (1.0 / 3.0) * parameters[2] * temp * temp2 +
          0.25 * parameters[3] * temp2 * temp2;
      DF = (parameters[0] * temp + parameters[2] * temp2 + parameters[3] * temp * temp2) * DTDX;
      DDF = parameters[0] * DTDX * DTDX + parameters[0] * temp * cos_theta * cube(DTDX) +
            2.0 * parameters[2] * temp * DTDX * DTDX + parameters[2] * temp2 * cos_theta * cube(DTDX) +
            3.0 * parameters[3] * temp2 * DTDX * DTDX + parameters[3] * temp * temp2 * cos_theta * cube(DTDX);
      break;
    case BendType::CFF_Quartic:
      temp = theta - parameters[1];
      temp2 = temp * temp;
      U = parameters[0] * temp2 + parameters[2] * temp * temp2 + parameters[3] * temp2 * temp2;
      DF = (2.0 * parameters[0] * temp + 3.0 * parameters[2] * temp2 + 4.0 * parameters[3] * temp * temp2) * DTDX;
      DDF = 2.0 * parameters[0] * DTDX * DTDX + 2.0 * parameters[0] * temp * cos_theta * cube(DTDX) +
            6.0 * parameters[2] * temp * DTDX * DTDX + 3.0 * parameters[2] * temp2 * cos_theta * cube(DTDX) +
            12.0 * parameters[3] * temp2 * DTDX * DTDX + 4.0 * parameters[3] * temp * temp2 * cos_theta * cube(DTDX);
      break;
    case BendType::HarmonicCosine:
      temp = cos_theta - parameters[1];
      temp2 = temp * temp;
      U = 0.5 * parameters[0] * temp2;
      DF = parameters[0] * temp;
      DDF = parameters[0];
      break;
    case BendType::Cosine:
      temp = parameters[1] * theta - parameters[2];
      U = parameters[0] * (1.0 + std::cos(temp));
      DF = -(parameters[0] * parameters[1] * std::sin(temp)) * DTDX;
      DDF = -parameters[0] * parameters[1] *
            (parameters[1] * std::cos(temp) + std::sin(temp) * cos_theta * DTDX) * DTDX * DTDX;
      break;
    case BendType::Tafipolsky:
      U = 0.5 * parameters[0] * (1.0 + std::cos(theta)) * (1.0 + std::cos(2.0 * theta));
      DF = parameters[0] * cos_theta * (2.0 + 3.0 * cos_theta);
      DDF = 2.0 * parameters[0] * (1.0 + 3.0 * cos_theta);
      break;
    case BendType::MM3:
    case BendType::MM3_inplane:
      temp = (theta - parameters[1]) * Units::RadiansToDegrees;
      temp2 = temp * temp;
      U = parameters[0] * temp2 *
          (1.0 - 0.014 * temp + 5.6e-5 * temp2 - 7.0e-7 * temp * temp2 + 2.2e-8 * temp2 * temp2);
      DF = parameters[0] * Units::RadiansToDegrees *
           (2.0 - (3.0 * 0.014 - (4.0 * 5.6e-5 - (5.0 * 7.0e-7 - 6.0 * 2.2e-8 * temp) * temp) * temp) * temp) * temp *
           DTDX;
      DDF = 0.0;
      break;
    default:
      std::unreachable();
  }

  geometry.DF = DF;
  geometry.DDF = DDF;
  geometry.Rab = dr_ab;
  geometry.Rbc = dr_cb;
  geometry.cosTheta = cos_theta;
  geometry.rab = r_ab;
  geometry.rbc = r_cb;

  geometry.dtA.x = (dr_cb.x - cos_theta * dr_ab.x) / r_ab;
  geometry.dtA.y = (dr_cb.y - cos_theta * dr_ab.y) / r_ab;
  geometry.dtA.z = (dr_cb.z - cos_theta * dr_ab.z) / r_ab;

  geometry.dtC.x = (dr_ab.x - cos_theta * dr_cb.x) / r_cb;
  geometry.dtC.y = (dr_ab.y - cos_theta * dr_cb.y) / r_cb;
  geometry.dtC.z = (dr_ab.z - cos_theta * dr_cb.z) / r_cb;

  geometry.dtB.x = -(geometry.dtA.x + geometry.dtC.x);
  geometry.dtB.y = -(geometry.dtA.y + geometry.dtC.y);
  geometry.dtB.z = -(geometry.dtA.z + geometry.dtC.z);

  const double3 du_da = DF * geometry.dtA;
  const double3 du_db = DF * geometry.dtB;
  const double3 du_dc = DF * geometry.dtC;

  double3x3 strain_derivative{};
  strain_derivative.ax = r_ab * dr_ab.x * du_da.x + r_cb * dr_cb.x * du_dc.x;
  strain_derivative.bx = r_ab * dr_ab.y * du_da.x + r_cb * dr_cb.y * du_dc.x;
  strain_derivative.cx = r_ab * dr_ab.z * du_da.x + r_cb * dr_cb.z * du_dc.x;

  strain_derivative.ay = r_ab * dr_ab.x * du_da.y + r_cb * dr_cb.x * du_dc.y;
  strain_derivative.by = r_ab * dr_ab.y * du_da.y + r_cb * dr_cb.y * du_dc.y;
  strain_derivative.cy = r_ab * dr_ab.z * du_da.y + r_cb * dr_cb.z * du_dc.y;

  strain_derivative.az = r_ab * dr_ab.x * du_da.z + r_cb * dr_cb.x * du_dc.z;
  strain_derivative.bz = r_ab * dr_ab.y * du_da.z + r_cb * dr_cb.y * du_dc.z;
  strain_derivative.cz = r_ab * dr_ab.z * du_da.z + r_cb * dr_cb.z * du_dc.z;

  return {U, {du_da, du_db, du_dc}, strain_derivative, geometry};
}
