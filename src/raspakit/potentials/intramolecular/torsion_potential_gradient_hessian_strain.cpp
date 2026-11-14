module;

module torsion_potential_gradient_hessian_strain;

import std;

import double3;
import double3x3;
import torsion_potential;
import minimization_torsion_hessian_geometry;
import dihedral_derivatives;

std::tuple<double, std::array<double3, 4>, double3x3, Minimization::TorsionHessianGeometry>
Potentials::Internal::torsionPotentialEnergyGradientHessianStrain(
    TorsionType type, const std::array<double, maximumNumberOfTorsionParameters> &parameters, const double3 &posA,
    const double3 &posB, const double3 &posC, const double3 &posD)
{
  const double3 Dab = posA - posB;
  const double3 Dcb_raw = posC - posB;
  const double rbc = std::sqrt(double3::dot(Dcb_raw, Dcb_raw));
  const double3 Dcb = Dcb_raw / rbc;
  const double3 Ddc = posD - posC;

  const double dot_ab = double3::dot(Dab, Dcb);
  const double dot_cd = double3::dot(Ddc, Dcb);

  double3 dr = Dab - dot_ab * Dcb;
  const double r = std::sqrt(double3::dot(dr, dr));
  dr /= r;

  double3 ds = Ddc - dot_cd * Dcb;
  const double s = std::sqrt(double3::dot(ds, ds));
  ds /= s;

  // Phi is defined in protein convention Phi(trans) = Pi.
  double cos_phi = std::clamp(double3::dot(dr, ds), -1.0, 1.0);
  const double cos_phi2 = cos_phi * cos_phi;


  double U{};
  double DF{};
  double DDF{};
  // a potential written in phi: dU/dphi and d2U/dphi2, taken with the exact derivatives of phi (see below)
  std::optional<std::pair<double, double>> phiDerivatives{};
  const auto exactPhi = [&]() { return Potentials::dihedralAngleDerivatives(posA, posB, posC, posD, true); };
  std::optional<Potentials::DihedralDerivatives> dihedral{};

  switch (type)
  {
    case TorsionType::Fixed:
    case TorsionType::CVFFBlocked:
      break;
    case TorsionType::Harmonic:
    {
      // (1/2)*p_0*(phi-p_1)^2
      dihedral = exactPhi();
      double phi = dihedral->phi - parameters[1];
      phi -= std::rint(phi / (2.0 * std::numbers::pi)) * 2.0 * std::numbers::pi;
      U = 0.5 * parameters[0] * phi * phi;
      phiDerivatives = std::pair{parameters[0] * phi, parameters[0]};
      break;
    }
    case TorsionType::HarmonicCosine:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      U = 0.5 * parameters[0] * (cos_phi - parameters[1]) * (cos_phi - parameters[1]);
      DF = parameters[0] * (cos_phi - parameters[1]);
      DDF = parameters[0];
      break;
    case TorsionType::ThreeCosine:
    case TorsionType::MM3:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      U = 0.5 * parameters[0] * (1.0 + cos_phi) + parameters[1] * (1.0 - cos_phi2) +
          0.5 * parameters[2] * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2);
      DF = 0.5 * parameters[0] - 2.0 * parameters[1] * cos_phi + 1.5 * parameters[2] * (4.0 * cos_phi2 - 1.0);
      DDF = -2.0 * (parameters[1] - 6.0 * parameters[2] * cos_phi);
      break;
    case TorsionType::RyckaertBellemans:
      // Sum_i=0^5 p_i*cos(phi)^i (with alternating signs, polymer convention)
      U = parameters[0] - parameters[1] * cos_phi + parameters[2] * cos_phi2 - parameters[3] * cos_phi * cos_phi2 +
          parameters[4] * cos_phi2 * cos_phi2 - parameters[5] * cos_phi2 * cos_phi2 * cos_phi;
      DF = -parameters[1] + 2.0 * parameters[2] * cos_phi - 3.0 * parameters[3] * cos_phi2 +
           4.0 * parameters[4] * cos_phi2 * cos_phi - 5.0 * parameters[5] * cos_phi2 * cos_phi2;
      DDF = 2.0 * parameters[2] - 6.0 * parameters[3] * cos_phi + 12.0 * parameters[4] * cos_phi2 -
            20.0 * parameters[5] * cos_phi2 * cos_phi;
      break;
    case TorsionType::TraPPE:
      // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
      U = parameters[0] + (1.0 + cos_phi) * (parameters[1] + parameters[3] -
                                             2.0 * (cos_phi - 1.0) * (parameters[2] - 2.0 * parameters[3] * cos_phi));
      DF = parameters[1] - 4.0 * parameters[2] * cos_phi + 3.0 * parameters[3] * (4.0 * cos_phi2 - 1.0);
      DDF = -4.0 * (parameters[2] - 6.0 * parameters[3] * cos_phi);
      break;
    case TorsionType::TraPPE_Extended:
      // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
      U = parameters[0] - parameters[2] + parameters[4] + (parameters[1] - 3.0 * parameters[3]) * cos_phi +
          (2.0 * parameters[2] - 8.0 * parameters[4]) * cos_phi2 + 4.0 * parameters[3] * cos_phi2 * cos_phi +
          8.0 * parameters[4] * cos_phi2 * cos_phi2;
      DF = parameters[1] - 3.0 * parameters[3] + 4.0 * (parameters[2] - 4.0 * parameters[4]) * cos_phi +
           12.0 * parameters[3] * cos_phi2 + 32.0 * parameters[4] * cos_phi2 * cos_phi;
      DDF = 4.0 * parameters[2] - 16.0 * parameters[4] + 24.0 * parameters[3] * cos_phi +
            96.0 * parameters[4] * cos_phi2;
      break;
    case TorsionType::ModifiedTraPPE:
    {
      // p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
      dihedral = exactPhi();
      double phi = dihedral->phi - parameters[4];
      phi -= std::rint(phi / (2.0 * std::numbers::pi)) * 2.0 * std::numbers::pi;
      const double shifted_cos_phi = std::cos(phi);
      const double shifted_cos_phi2 = shifted_cos_phi * shifted_cos_phi;
      U = parameters[0] + parameters[1] + 2.0 * parameters[2] + parameters[3] +
          (parameters[1] - 3.0 * parameters[3]) * shifted_cos_phi - 2.0 * parameters[2] * shifted_cos_phi2 +
          4.0 * parameters[3] * shifted_cos_phi * shifted_cos_phi2;
      phiDerivatives = std::pair{-parameters[1] * std::sin(phi) + 2.0 * parameters[2] * std::sin(2.0 * phi) -
                                     3.0 * parameters[3] * std::sin(3.0 * phi),
                                 -parameters[1] * std::cos(phi) + 4.0 * parameters[2] * std::cos(2.0 * phi) -
                                     9.0 * parameters[3] * std::cos(3.0 * phi)};
      break;
    }
    case TorsionType::CVFF:
    {
      // p_0*(1+cos(p_1*phi-p_2))
      dihedral = exactPhi();
      const double argument = parameters[1] * dihedral->phi - parameters[2];
      U = parameters[0] * (1.0 + std::cos(argument));
      phiDerivatives = std::pair{-parameters[0] * parameters[1] * std::sin(argument),
                                 -parameters[0] * parameters[1] * parameters[1] * std::cos(argument)};
      break;
    }
    case TorsionType::CFF:
      // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
      U = parameters[0] * (1.0 - cos_phi) + 2.0 * parameters[1] * (1.0 - cos_phi2) +
          parameters[2] * (1.0 + 3.0 * cos_phi - 4.0 * cos_phi * cos_phi2);
      DF = -parameters[0] - 4.0 * parameters[1] * cos_phi + 3.0 * parameters[2] * (1.0 - 4.0 * cos_phi2);
      DDF = -4.0 * (parameters[1] + 6.0 * parameters[2] * cos_phi);
      break;
    case TorsionType::CFF2:
      // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
      U = parameters[0] * (1.0 + cos_phi) + parameters[2] +
          cos_phi * (-3.0 * parameters[2] + 2.0 * cos_phi * (parameters[1] + 2.0 * parameters[2] * cos_phi));
      DF = parameters[0] - 3.0 * parameters[2] + 4.0 * cos_phi * (parameters[1] + 3.0 * parameters[2] * cos_phi);
      DDF = 4.0 * (parameters[1] + 6.0 * parameters[2] * cos_phi);
      break;
    case TorsionType::OPLS:
      // (1/2)p_0+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
      U = 0.5 * (parameters[0] + (1.0 + cos_phi) * (parameters[1] + parameters[3] -
                                                    2.0 * (cos_phi - 1.0) * (parameters[2] - 2.0 * parameters[3] * cos_phi)));
      DF = 0.5 * parameters[1] - 2.0 * parameters[2] * cos_phi + 1.5 * parameters[3] * (4.0 * cos_phi2 - 1.0);
      DDF = -2.0 * (parameters[2] - 6.0 * parameters[3] * cos_phi);
      break;
    case TorsionType::FourierSeries:
      U = 0.5 * (parameters[0] + 2.0 * parameters[1] + parameters[2] + parameters[4] + 2.0 * parameters[5] +
                 (parameters[0] - 3.0 * parameters[2] + 5.0 * parameters[4]) * cos_phi -
                 2.0 * (parameters[1] - 4.0 * parameters[3] + 9.0 * parameters[5]) * cos_phi2 +
                 4.0 * (parameters[2] - 5.0 * parameters[4]) * cos_phi2 * cos_phi -
                 8.0 * (parameters[3] - 6.0 * parameters[5]) * cos_phi2 * cos_phi2 +
                 16.0 * parameters[4] * cos_phi2 * cos_phi2 * cos_phi -
                 32.0 * parameters[5] * cos_phi2 * cos_phi2 * cos_phi2);
      DF = 0.5 * (parameters[0] - 3.0 * parameters[2] + 5.0 * parameters[4]) -
           2.0 * (parameters[1] - 4.0 * parameters[3] + 9.0 * parameters[5]) * cos_phi +
           6.0 * (parameters[2] - 5.0 * parameters[4]) * cos_phi2 -
           16.0 * (parameters[3] - 6.0 * parameters[5]) * cos_phi2 * cos_phi +
           40.0 * parameters[4] * cos_phi2 * cos_phi2 - 96.0 * parameters[5] * cos_phi2 * cos_phi2 * cos_phi;
      DDF = -2.0 * parameters[1] + 8.0 * parameters[3] - 18.0 * parameters[5] +
            12.0 * (parameters[2] - 5.0 * parameters[4]) * cos_phi -
            48.0 * (parameters[3] - 6.0 * parameters[5]) * cos_phi2 +
            160.0 * parameters[4] * cos_phi2 * cos_phi - 480.0 * parameters[5] * cos_phi2 * cos_phi2;
      break;
    case TorsionType::FourierSeries2:
      U = 0.5 * (parameters[2] + 2.0 * parameters[3] + parameters[4] - 3.0 * parameters[2] * cos_phi +
                 5.0 * parameters[4] * cos_phi + parameters[0] * (1.0 + cos_phi) +
                 2.0 * (parameters[1] - parameters[1] * cos_phi2 +
                        cos_phi2 * (parameters[5] * (3.0 - 4.0 * cos_phi2) * (3.0 - 4.0 * cos_phi2) +
                                    4.0 * parameters[3] * (cos_phi2 - 1.0) +
                                    2.0 * cos_phi * (parameters[2] + parameters[4] * (4.0 * cos_phi2 - 5.0)))));
      DF = 0.5 * parameters[0] + parameters[2] * (6.0 * cos_phi2 - 1.5) +
           parameters[4] * (2.5 - 30.0 * cos_phi2 + 40.0 * cos_phi2 * cos_phi2) +
           cos_phi * (-2.0 * parameters[1] + parameters[3] * (16.0 * cos_phi2 - 8.0) +
                      parameters[5] * (18.0 - 96.0 * cos_phi2 + 96.0 * cos_phi2 * cos_phi2));
      DDF = -2.0 * parameters[1] + 12.0 * parameters[2] * cos_phi + parameters[3] * (48.0 * cos_phi2 - 8.0) +
            parameters[4] * cos_phi * (160.0 * cos_phi2 - 60.0) +
            parameters[5] * (18.0 - 288.0 * cos_phi2 + 480.0 * cos_phi2 * cos_phi2);
      break;
    case TorsionType::Polynomial:
      // U = sum_i p_i c^i; DF = U'(c); DDF = U''(c), all by Horner's rule
      U = parameters[5];
      DF = 0.0;
      DDF = 0.0;
      for (std::size_t i = 5; i-- > 0;)
      {
        DDF = DDF * cos_phi + 2.0 * DF;
        DF = DF * cos_phi + U;
        U = U * cos_phi + parameters[i];
      }
      break;
    default:
      std::unreachable();
  }

  const double d = dot_ab / rbc;
  const double e = dot_cd / rbc;

  const double3 dtA = (ds - cos_phi * dr) / r;
  const double3 dtD = (dr - cos_phi * ds) / s;
  const double3 dtB = dtA * (d - 1.0) + e * dtD;
  const double3 dtC = -dtD * (e + 1.0) - d * dtA;

  double3 du_da = DF * dtA;
  double3 du_db = DF * dtB;
  double3 du_dc = DF * dtC;
  double3 du_dd = DF * dtD;
  std::optional<std::array<std::array<double, 12>, 12>> cartesianHessian{};
  if (phiDerivatives.has_value())
  {
    // gradient U' grad(phi), Hessian U'' grad(phi) grad(phi)^T + U' Hess(phi): finite at cis and trans
    const auto [firstDerivative, secondDerivative] = phiDerivatives.value();
    du_da = firstDerivative * dihedral->gradient[0];
    du_db = firstDerivative * dihedral->gradient[1];
    du_dc = firstDerivative * dihedral->gradient[2];
    du_dd = firstDerivative * dihedral->gradient[3];
    std::array<std::array<double, 12>, 12> cartesian{};
    for (std::size_t i = 0; i < 12; ++i)
    {
      const double gradientI = (&dihedral->gradient[i / 3].x)[i % 3];
      for (std::size_t j = 0; j < 12; ++j)
      {
        const double gradientJ = (&dihedral->gradient[j / 3].x)[j % 3];
        cartesian[i][j] = secondDerivative * gradientI * gradientJ + firstDerivative * dihedral->hessian[i][j];
      }
    }
    cartesianHessian = cartesian;
  }

  // Note: Dcb was normalized, so its contribution is scaled back by rbc.
  double3x3 strain_derivative{};
  strain_derivative.ax = Dab.x * du_da.x + Dcb.x * rbc * (du_dc.x + du_dd.x) + Ddc.x * du_dd.x;
  strain_derivative.bx = Dab.y * du_da.x + Dcb.y * rbc * (du_dc.x + du_dd.x) + Ddc.y * du_dd.x;
  strain_derivative.cx = Dab.z * du_da.x + Dcb.z * rbc * (du_dc.x + du_dd.x) + Ddc.z * du_dd.x;
  strain_derivative.ay = Dab.x * du_da.y + Dcb.x * rbc * (du_dc.y + du_dd.y) + Ddc.x * du_dd.y;
  strain_derivative.by = Dab.y * du_da.y + Dcb.y * rbc * (du_dc.y + du_dd.y) + Ddc.y * du_dd.y;
  strain_derivative.cy = Dab.z * du_da.y + Dcb.z * rbc * (du_dc.y + du_dd.y) + Ddc.z * du_dd.y;
  strain_derivative.az = Dab.x * du_da.z + Dcb.x * rbc * (du_dc.z + du_dd.z) + Ddc.x * du_dd.z;
  strain_derivative.bz = Dab.y * du_da.z + Dcb.y * rbc * (du_dc.z + du_dd.z) + Ddc.y * du_dd.z;
  strain_derivative.cz = Dab.z * du_da.z + Dcb.z * rbc * (du_dc.z + du_dd.z) + Ddc.z * du_dd.z;

  Minimization::TorsionHessianGeometry geometry{};
  geometry.DF = DF;
  geometry.DDF = DDF;
  geometry.dtA = dtA;
  geometry.dtB = dtB;
  geometry.dtC = dtC;
  geometry.dtD = dtD;
  geometry.Dab = Dab;
  geometry.Dcb = Dcb;
  geometry.Ddc = Ddc;
  geometry.rbc = rbc;
  geometry.dr = dr;
  geometry.ds = ds;
  geometry.r = r;
  geometry.s = s;
  geometry.d = d;
  geometry.e = e;
  geometry.cosPhi = cos_phi;
  geometry.cartesianHessian = cartesianHessian;

  return {U, {du_da, du_db, du_dc, du_dd}, strain_derivative, geometry};
}
