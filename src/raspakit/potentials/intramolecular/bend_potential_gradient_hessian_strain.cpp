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

  // The sine from the cross product and the angle from the two together: the arc cosine resolves an angle near
  // 0 or 180 degrees only to about 1.5e-8 rad, and the square root of 1 - cos^2 no better.
  const double sin_true = double3::cross(dr_ab, dr_cb).length();
  const double theta = std::atan2(sin_true, cos_theta);

  double U{};
  double DF{};
  double DDF{};
  double temp{};
  double temp2{};

  // The forms written in theta: U and its first and second derivative in theta at the angle given. False for a
  // form that is not one of them.
  const auto thetaForm = [&](double angle, double &energy, double &first, double &second) -> bool
  {
    double t{};
    double t2{};
    switch (type)
    {
      case BendType::Harmonic:
      case BendType::CoreShell:
        t = angle - parameters[1];
        energy = 0.5 * parameters[0] * t * t;
        first = parameters[0] * t;
        second = parameters[0];
        return true;
      case BendType::Quartic:
        t = angle - parameters[1];
        t2 = t * t;
        energy = 0.5 * parameters[0] * t2 + (1.0 / 3.0) * parameters[2] * t * t2 + 0.25 * parameters[3] * t2 * t2;
        first = parameters[0] * t + parameters[2] * t2 + parameters[3] * t * t2;
        second = parameters[0] + 2.0 * parameters[2] * t + 3.0 * parameters[3] * t2;
        return true;
      case BendType::CFF_Quartic:
        t = angle - parameters[1];
        t2 = t * t;
        energy = parameters[0] * t2 + parameters[2] * t * t2 + parameters[3] * t2 * t2;
        first = 2.0 * parameters[0] * t + 3.0 * parameters[2] * t2 + 4.0 * parameters[3] * t * t2;
        second = 2.0 * parameters[0] + 6.0 * parameters[2] * t + 12.0 * parameters[3] * t2;
        return true;
      case BendType::Cosine:
        t = parameters[1] * angle - parameters[2];
        energy = parameters[0] * (1.0 + std::cos(t));
        first = -parameters[0] * parameters[1] * std::sin(t);
        second = -parameters[0] * parameters[1] * parameters[1] * std::cos(t);
        return true;
      case BendType::MM3:
      case BendType::MM3_inplane:
      {
        // U(t) = K (t^2 - a t^3 + b t^4 - c t^5 + d t^6), t = (theta - theta_0) in degrees
        t = (angle - parameters[1]) * Units::RadiansToDegrees;
        t2 = t * t;
        energy = parameters[0] * t2 * (1.0 - 0.014 * t + 5.6e-5 * t2 - 7.0e-7 * t * t2 + 2.2e-8 * t2 * t2);
        const double dUdt =
            parameters[0] * (2.0 - (3.0 * 0.014 - (4.0 * 5.6e-5 - (5.0 * 7.0e-7 - 6.0 * 2.2e-8 * t) * t) * t) * t) * t;
        const double d2Udt2 = parameters[0] *
                              (2.0 - (6.0 * 0.014 - (12.0 * 5.6e-5 - (20.0 * 7.0e-7 - 30.0 * 2.2e-8 * t) * t) * t) * t);
        first = dUdt * Units::RadiansToDegrees;
        second = d2Udt2 * Units::RadiansToDegrees * Units::RadiansToDegrees;
        return true;
      }
      default:
        return false;
    }
  };

  // dU/dcos and d2U/dcos2 from the derivatives in theta, at an angle whose sine is given.
  const auto inCosine = [](double first, double second, double cosine, double sine, double &dF, double &ddF)
  {
    const double dtdx = -1.0 / sine;
    dF = first * dtdx;
    ddF = second * dtdx * dtdx + first * cosine * cube(dtdx);
  };

  double dUdTheta{};
  double d2UdTheta2{};
  if (thetaForm(theta, U, dUdTheta, d2UdTheta2))
  {
    // How near a linear geometry the two terms of d2U/dcos2, each as large as 1/sin^2, still leave their
    // difference: beyond it the form is read at this distance from the line instead, where it is smooth.
    constexpr double linearWindow = 1.0e-3;
    if (sin_true >= linearWindow)
    {
      inCosine(dUdTheta, d2UdTheta2, cos_theta, sin_true, DF, DDF);
    }
    else
    {
      // A form whose minimum (or any stationary point) is the linear geometry itself has dU/dtheta falling off
      // as the distance from it, and dU/dcos has a limit there: U''. The sine used to be held at 1e-8 and the
      // angle taken from the arc cosine, so on the line dU/dcos was nought times 1e8, which is nought, and the
      // Hessian, which there is dU/dcos times the second derivative of the cosine and nothing else, was nought
      // with it: a linear molecule had no bending curvature and bending modes of no frequency.
      const double delta = std::asin(sin_true);
      const double towards = cos_theta < 0.0 ? 1.0 : -1.0;  // dU/dtheta = -+ U'' delta about 180 or 0 degrees
      const bool stationary = std::fabs(dUdTheta + towards * d2UdTheta2 * delta) <=
                              1.0e-3 * std::fabs(d2UdTheta2) * delta + 1.0e-14 * std::fabs(d2UdTheta2);
      if (stationary)
      {
        double unusedEnergy{};
        double firstThere{};
        double secondThere{};
        const double there = cos_theta < 0.0 ? std::numbers::pi - linearWindow : linearWindow;
        thetaForm(there, unusedEnergy, firstThere, secondThere);
        double unusedFirst{};
        inCosine(firstThere, secondThere, std::cos(there), std::sin(there), unusedFirst, DDF);
        DF = sin_true > 1.0e-12 ? -dUdTheta / sin_true : towards * d2UdTheta2;
      }
      else
      {
        // not stationary on the line: dU/dcos has no limit there, and is bounded as before
        inCosine(dUdTheta, d2UdTheta2, cos_theta, std::max(1.0e-8, sin_true), DF, DDF);
      }
    }
  }
  else
  {
    switch (type)
    {
      case BendType::Fixed:
      case BendType::Rigid:
        U = 0.0;
        DF = 0.0;
        DDF = 0.0;
        break;
      case BendType::HarmonicCosine:
        temp = cos_theta - parameters[1];
        temp2 = temp * temp;
        U = 0.5 * parameters[0] * temp2;
        DF = parameters[0] * temp;
        DDF = parameters[0];
        break;
      case BendType::Tafipolsky:
        U = 0.5 * parameters[0] * (1.0 + std::cos(theta)) * (1.0 + std::cos(2.0 * theta));
        DF = parameters[0] * cos_theta * (2.0 + 3.0 * cos_theta);
        DDF = 2.0 * parameters[0] * (1.0 + 3.0 * cos_theta);
        break;
      default:
        std::unreachable();
    }
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
