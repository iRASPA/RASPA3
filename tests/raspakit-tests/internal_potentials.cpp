#include <gtest/gtest.h>

import std;

import double3;
import units;
import bond_potential;
import urey_bradley_potential;
import bend_potential;
import inversion_bend_potential;
import torsion_potential;
import bond_bond_potential;
import bond_bend_potential;
import bond_torsion_potential;
import bend_bend_potential;
import bend_torsion_potential;

// Parity tests for the internal (intramolecular) potentials against the functional forms in
// RASPA2 'internal_energy.c'. The reference values below are computed with independent
// implementations of the RASPA2 formulas, written directly in this file.

namespace
{

constexpr double DEG2RAD = std::numbers::pi / 180.0;
constexpr double RAD2DEG = 180.0 / std::numbers::pi;

// Generic, non-degenerate test geometry.
const double3 posA{1.10, 0.20, -0.30};
const double3 posB{0.00, 0.00, 0.00};
const double3 posC{1.50, 1.40, 0.10};
const double3 posD{2.30, 1.20, 1.00};

double distance(const double3 &a, const double3 &b)
{
  double3 dr = a - b;
  return std::sqrt(double3::dot(dr, dr));
}

// Bend angle A-B-C as defined in RASPA2 CalculateBendEnergy.
double bendAngle(const double3 &a, const double3 &b, const double3 &c)
{
  double3 rab = (a - b).normalized();
  double3 rbc = (c - b).normalized();
  double cos_theta = std::clamp(double3::dot(rab, rbc), -1.0, 1.0);
  return std::acos(cos_theta);
}

// Dihedral angle A-B-C-D (protein convention), sign included, as in RASPA2 CalculateTorsionEnergy.
double dihedralAngle(const double3 &a, const double3 &b, const double3 &c, const double3 &d)
{
  double3 Dab = a - b;
  double3 Dbc = (c - b).normalized();
  double3 Dcd = d - c;

  double dot_ab = double3::dot(Dab, Dbc);
  double dot_cd = double3::dot(Dcd, Dbc);

  double3 dr = (Dab - dot_ab * Dbc).normalized();
  double3 ds = (Dcd - dot_cd * Dbc).normalized();

  double cos_phi = std::clamp(double3::dot(dr, ds), -1.0, 1.0);

  double3 Pb = double3::cross(Dab, Dbc);
  double3 Pc = double3::cross(Dbc, Dcd);
  // RASPA2 sign, written component-wise, equals Dbc . (Pb x Pc)
  double sign = double3::dot(Dbc, double3::cross(Pb, Pc));
  return std::copysign(std::acos(cos_phi), sign);
}

// Inversion angle chi as defined in RASPA2 CalculateInversionBendEnergy.
// planeBCD=true : reference plane through B-C-D; false: through A-C-D.
double inversionCosChi(const double3 &a, const double3 &b, const double3 &c, const double3 &d, bool planeBCD)
{
  double3 Rab = a - b;
  double3 Rbc = c - b;
  double3 Rbd = d - b;
  double3 Rcd = d - c;
  double3 Rad = d - a;

  double cval;
  if (planeBCD)
  {
    double dot = double3::dot(Rbc, Rbd);
    cval = double3::dot(Rbc, Rbc) * double3::dot(Rbd, Rbd) - dot * dot;
  }
  else
  {
    double dot = double3::dot(Rad, Rcd);
    cval = double3::dot(Rcd, Rcd) * double3::dot(Rad, Rad) - dot * dot;
  }

  double e = double3::dot(Rab, double3::cross(Rbd, Rbc));
  double rrab = std::sqrt(double3::dot(Rab, Rab));
  double cos_chi = std::sqrt(double3::dot(Rab, Rab) - e * e / cval) / rrab;
  return std::clamp(cos_chi, -1.0, 1.0);
}

// RASPA2 Smoothing() from utils.c
double smoothing(double theta)
{
  double on = 170.0 * DEG2RAD;
  double off = 180.0 * DEG2RAD;
  if (theta < on) return 1.0;
  return (off - theta) * (off - theta) * (off + 2.0 * theta - 3.0 * on) / ((off - on) * (off - on) * (off - on));
}

double K2E(double x) { return x * Units::KelvinToEnergy; }

}  // namespace

TEST(internal_potentials, bond_forms_match_RASPA2)
{
  double r = distance(posA, posB);
  double tol = 1e-8;

  // HARMONIC: 0.5*p0*(r-p1)^2
  {
    BondPotential bond({0, 1}, BondType::Harmonic, {96500.0, 1.10});
    double ref = 0.5 * K2E(96500.0) * (r - 1.10) * (r - 1.10);
    EXPECT_NEAR(bond.calculateEnergy(posA, posB), ref, tol * std::abs(ref));
  }

  // CORE_SHELL_SPRING: 0.5*p0*r^2
  {
    BondPotential bond({0, 1}, BondType::CoreShellSpring, {5000.0});
    double ref = 0.5 * K2E(5000.0) * r * r;
    EXPECT_NEAR(bond.calculateEnergy(posA, posB), ref, tol * std::abs(ref));
  }

  // MORSE: p0*[(1-exp(-p1*(r-p2)))^2-1]
  {
    BondPotential bond({0, 1}, BondType::Morse, {600.0, 2.0, 1.2});
    double temp = std::exp(2.0 * (1.2 - r));
    double ref = K2E(600.0) * ((1.0 - temp) * (1.0 - temp) - 1.0);
    EXPECT_NEAR(bond.calculateEnergy(posA, posB), ref, tol * std::abs(ref));
  }

  // LJ_12_6: p0/r^12 - p1/r^6
  {
    BondPotential bond({0, 1}, BondType::LJ_12_6, {1.0e7, 3.0e3});
    double rr = r * r;
    double t = 1.0 / (rr * rr * rr);
    double ref = K2E(1.0e7) * t * t - K2E(3.0e3) * t;
    EXPECT_NEAR(bond.calculateEnergy(posA, posB), ref, tol * std::abs(ref));
  }

  // LENNARD_JONES: 4*p0*((p1/r)^12-(p1/r)^6)
  {
    BondPotential bond({0, 1}, BondType::LennardJones, {120.0, 1.2});
    double t = std::pow((1.2 * 1.2) / (r * r), 3);
    double ref = 4.0 * K2E(120.0) * (t * (t - 1.0));
    EXPECT_NEAR(bond.calculateEnergy(posA, posB), ref, tol * std::abs(ref));
  }

  // BUCKINGHAM: p0*exp(-p1*r)-p2/r^6
  {
    BondPotential bond({0, 1}, BondType::Buckingham, {3.0e5, 3.5, 1.2e4});
    double rr = r * r;
    double ref = K2E(3.0e5) * std::exp(-3.5 * r) - K2E(1.2e4) / (rr * rr * rr);
    EXPECT_NEAR(bond.calculateEnergy(posA, posB), ref, tol * std::abs(ref));
  }

  // RESTRAINED_HARMONIC
  {
    BondPotential bond({0, 1}, BondType::RestrainedHarmonic, {96500.0, 1.0, 0.1});
    double r1 = r - 1.0;
    double ref = 0.5 * K2E(96500.0) * std::pow(std::min(std::fabs(r1), 0.1), 2) +
                 K2E(96500.0) * 0.1 * std::max(std::fabs(r1) - 0.1, 0.0);
    EXPECT_NEAR(bond.calculateEnergy(posA, posB), ref, tol * std::abs(ref));
  }

  // QUARTIC: (1/2)p0*dr^2+(1/3)p2*dr^3+(1/4)p3*dr^4
  {
    BondPotential bond({0, 1}, BondType::Quartic, {96500.0, 1.0, -3000.0, 4000.0});
    double dr = r - 1.0;
    double ref = 0.5 * K2E(96500.0) * dr * dr + (1.0 / 3.0) * K2E(-3000.0) * dr * dr * dr +
                 0.25 * K2E(4000.0) * dr * dr * dr * dr;
    EXPECT_NEAR(bond.calculateEnergy(posA, posB), ref, tol * std::abs(ref));
  }

  // CFF_QUARTIC: p0*dr^2+p2*dr^3+p3*dr^4
  {
    BondPotential bond({0, 1}, BondType::CFF_Quartic, {96500.0, 1.0, -3000.0, 4000.0});
    double dr = r - 1.0;
    double ref = K2E(96500.0) * dr * dr + K2E(-3000.0) * dr * dr * dr + K2E(4000.0) * dr * dr * dr * dr;
    EXPECT_NEAR(bond.calculateEnergy(posA, posB), ref, tol * std::abs(ref));
  }

  // MM3: p0*dr^2*(1-2.55*dr+(7/12)*2.55^2*dr^2), p0 in mdyne/A (x 71.94 kcal/mol)
  {
    BondPotential bond({0, 1}, BondType::MM3, {4.49, 1.0});
    double dr = r - 1.0;
    double p0 = 4.49 * 71.94 * Units::KCalPerMolToEnergy;
    double ref = p0 * dr * dr * (1.0 - 2.55 * dr + (7.0 / 12.0) * 2.55 * 2.55 * dr * dr);
    EXPECT_NEAR(bond.calculateEnergy(posA, posB), ref, tol * std::abs(ref));
  }
}

TEST(internal_potentials, urey_bradley_forms_match_RASPA2)
{
  // Urey-Bradley is a two-body distance potential between atoms A and C of a bend.
  double r = distance(posA, posC);
  double tol = 1e-8;

  {
    UreyBradleyPotential ub({0, 2}, UreyBradleyType::Harmonic, {40000.0, 1.5});
    double ref = 0.5 * K2E(40000.0) * (r - 1.5) * (r - 1.5);
    EXPECT_NEAR(ub.calculateEnergy(posA, posC), ref, tol * std::abs(ref));
  }

  {
    UreyBradleyPotential ub({0, 2}, UreyBradleyType::Morse, {600.0, 2.0, 1.2});
    double temp = std::exp(2.0 * (1.2 - r));
    double ref = K2E(600.0) * ((1.0 - temp) * (1.0 - temp) - 1.0);
    EXPECT_NEAR(ub.calculateEnergy(posA, posC), ref, tol * std::abs(ref));
  }

  {
    UreyBradleyPotential ub({0, 2}, UreyBradleyType::MM3, {4.49, 1.4});
    double dr = r - 1.4;
    double p0 = 4.49 * 71.94 * Units::KCalPerMolToEnergy;
    double ref = p0 * dr * dr * (1.0 - 2.55 * dr + (7.0 / 12.0) * 2.55 * 2.55 * dr * dr);
    EXPECT_NEAR(ub.calculateEnergy(posA, posC), ref, tol * std::abs(ref));
  }
}

TEST(internal_potentials, bend_forms_match_RASPA2)
{
  double theta = bendAngle(posA, posB, posC);
  double cos_theta = std::cos(theta);
  double tol = 1e-8;

  // HARMONIC: 0.5*p0*(theta-p1)^2, p1 in degrees
  {
    BendPotential bend({0, 1, 2}, BendType::Harmonic, {62500.0, 114.0});
    double ref = 0.5 * K2E(62500.0) * std::pow(theta - 114.0 * DEG2RAD, 2);
    EXPECT_NEAR(bend.calculateEnergy(posA, posB, posC, std::nullopt), ref, tol * std::abs(ref));
  }

  // CORE_SHELL
  {
    BendPotential bend({0, 1, 2}, BendType::CoreShell, {62500.0, 114.0});
    double ref = 0.5 * K2E(62500.0) * std::pow(theta - 114.0 * DEG2RAD, 2);
    EXPECT_NEAR(bend.calculateEnergy(posA, posB, posC, std::nullopt), ref, tol * std::abs(ref));
  }

  // QUARTIC
  {
    BendPotential bend({0, 1, 2}, BendType::Quartic, {62500.0, 114.0, -1000.0, 2000.0});
    double dt = theta - 114.0 * DEG2RAD;
    double ref = 0.5 * K2E(62500.0) * dt * dt + (1.0 / 3.0) * K2E(-1000.0) * dt * dt * dt +
                 0.25 * K2E(2000.0) * dt * dt * dt * dt;
    EXPECT_NEAR(bend.calculateEnergy(posA, posB, posC, std::nullopt), ref, tol * std::abs(ref));
  }

  // CFF_QUARTIC
  {
    BendPotential bend({0, 1, 2}, BendType::CFF_Quartic, {62500.0, 114.0, -1000.0, 2000.0});
    double dt = theta - 114.0 * DEG2RAD;
    double ref = K2E(62500.0) * dt * dt + K2E(-1000.0) * dt * dt * dt + K2E(2000.0) * dt * dt * dt * dt;
    EXPECT_NEAR(bend.calculateEnergy(posA, posB, posC, std::nullopt), ref, tol * std::abs(ref));
  }

  // HARMONIC_COSINE: 0.5*p0*(cos(theta)-cos(p1))^2
  {
    BendPotential bend({0, 1, 2}, BendType::HarmonicCosine, {1000.0, 114.0});
    double dt = cos_theta - std::cos(114.0 * DEG2RAD);
    double ref = 0.5 * K2E(1000.0) * dt * dt;
    EXPECT_NEAR(bend.calculateEnergy(posA, posB, posC, std::nullopt), ref, tol * std::abs(ref));
  }

  // COSINE: p0*(1+cos(p1*theta-p2))
  {
    BendPotential bend({0, 1, 2}, BendType::Cosine, {1000.0, 2.0, 60.0});
    double ref = K2E(1000.0) * (1.0 + std::cos(2.0 * theta - 60.0 * DEG2RAD));
    EXPECT_NEAR(bend.calculateEnergy(posA, posB, posC, std::nullopt), ref, tol * std::abs(ref));
  }

  // TAFIPOLSKY: 0.5*p0*(1+cos(theta))*(1+cos(2*theta))
  {
    BendPotential bend({0, 1, 2}, BendType::Tafipolsky, {1000.0});
    double ref = 0.5 * K2E(1000.0) * (1.0 + std::cos(theta)) * (1.0 + std::cos(2.0 * theta));
    EXPECT_NEAR(bend.calculateEnergy(posA, posB, posC, std::nullopt), ref, tol * std::abs(ref));
  }

  // MM3: p0*dt^2*(1-0.014*dt+5.6e-5*dt^2-7e-7*dt^3+2.2e-8*dt^4), dt in degrees, p0 in mdyne A/rad^2
  {
    BendPotential bend({0, 1, 2}, BendType::MM3, {0.59, 114.0});
    double dt = RAD2DEG * (theta - 114.0 * DEG2RAD);
    double p0 = 0.59 * 0.021914 * Units::KCalPerMolToEnergy;
    double ref =
        p0 * dt * dt * (1.0 - 0.014 * dt + 5.6e-5 * dt * dt - 7.0e-7 * dt * dt * dt + 2.2e-8 * dt * dt * dt * dt);
    EXPECT_NEAR(bend.calculateEnergy(posA, posB, posC, std::nullopt), ref, tol * std::abs(ref));
  }

  // MM3_IN_PLANE: the angle is computed after projecting B onto the A-C-D plane
  {
    BendPotential bend({0, 1, 2}, BendType::MM3_inplane, {0.59, 114.0});

    double3 Rad = posA - posD;
    double3 Rbd = posB - posD;
    double3 Rcd = posC - posD;
    double3 t = double3::cross(Rad, Rcd);
    double rt2 = double3::dot(t, t);
    double delta = -double3::dot(t, Rbd) / rt2;
    double3 ip = posB + delta * t;
    double3 ap = posA - ip;
    double3 cp = posC - ip;
    double cos_in_plane =
        std::clamp(double3::dot(ap, cp) / std::sqrt(double3::dot(ap, ap) * double3::dot(cp, cp)), -1.0, 1.0);
    double theta_in_plane = std::acos(cos_in_plane);

    double dt = RAD2DEG * (theta_in_plane - 114.0 * DEG2RAD);
    double p0 = 0.59 * 0.021914 * Units::KCalPerMolToEnergy;
    double ref =
        p0 * dt * dt * (1.0 - 0.014 * dt + 5.6e-5 * dt * dt - 7.0e-7 * dt * dt * dt + 2.2e-8 * dt * dt * dt * dt);
    EXPECT_NEAR(bend.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }
}

TEST(internal_potentials, inversion_bend_forms_match_RASPA2)
{
  double tol = 1e-8;

  // HARMONIC (plane B-C-D): 0.5*p0*(chi-p1)^2
  {
    InversionBendPotential inv({0, 1, 2, 3}, InversionBendType::Harmonic, {30000.0, 10.0});
    double chi = std::acos(inversionCosChi(posA, posB, posC, posD, true));
    double ref = 0.5 * K2E(30000.0) * std::pow(chi - 10.0 * DEG2RAD, 2);
    EXPECT_NEAR(inv.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // HARMONIC_COSINE (plane B-C-D): 0.5*p0*(cos(chi)-cos(p1))^2
  {
    InversionBendPotential inv({0, 1, 2, 3}, InversionBendType::HarmonicCosine, {30000.0, 10.0});
    double cos_chi = inversionCosChi(posA, posB, posC, posD, true);
    double dt = cos_chi - std::cos(10.0 * DEG2RAD);
    double ref = 0.5 * K2E(30000.0) * dt * dt;
    EXPECT_NEAR(inv.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // PLANAR (plane B-C-D): p0*(1-cos(chi))
  {
    InversionBendPotential inv({0, 1, 2, 3}, InversionBendType::Planar, {30000.0});
    double cos_chi = inversionCosChi(posA, posB, posC, posD, true);
    double ref = K2E(30000.0) * (1.0 - cos_chi);
    EXPECT_NEAR(inv.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // HARMONIC2/HARMONIC_COSINE2/PLANAR2 use the plane A-C-D
  {
    InversionBendPotential inv({0, 1, 2, 3}, InversionBendType::Harmonic2, {30000.0, 10.0});
    double chi = std::acos(inversionCosChi(posA, posB, posC, posD, false));
    double ref = 0.5 * K2E(30000.0) * std::pow(chi - 10.0 * DEG2RAD, 2);
    EXPECT_NEAR(inv.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  {
    InversionBendPotential inv({0, 1, 2, 3}, InversionBendType::Planar2, {30000.0});
    double cos_chi = inversionCosChi(posA, posB, posC, posD, false);
    double ref = K2E(30000.0) * (1.0 - cos_chi);
    EXPECT_NEAR(inv.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // MM3 (plane A-C-D): MM3 sextic expansion in degrees, factor 0.02191418 kcal/mol
  {
    InversionBendPotential inv({0, 1, 2, 3}, InversionBendType::MM3, {0.59, 10.0});
    double chi = std::acos(inversionCosChi(posA, posB, posC, posD, false));
    double dt = RAD2DEG * (chi - 10.0 * DEG2RAD);
    double p0 = 0.59 * 0.02191418 * Units::KCalPerMolToEnergy;
    double ref =
        p0 * dt * dt * (1.0 - 0.014 * dt + 5.6e-5 * dt * dt - 7.0e-7 * dt * dt * dt + 2.2e-8 * dt * dt * dt * dt);
    EXPECT_NEAR(inv.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }
}

TEST(internal_potentials, torsion_forms_match_RASPA2)
{
  double phi = dihedralAngle(posA, posB, posC, posD);
  double cos_phi = std::cos(phi);
  double cos_phi2 = cos_phi * cos_phi;
  double tol = 1e-8;

  // HARMONIC: 0.5*p0*(phi-p1)^2 (with periodic wrapping)
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::Harmonic, {5000.0, 35.0});
    double dphi = phi - 35.0 * DEG2RAD;
    dphi -= std::rint(dphi / (2.0 * std::numbers::pi)) * 2.0 * std::numbers::pi;
    double ref = 0.5 * K2E(5000.0) * dphi * dphi;
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // HARMONIC_COSINE: 0.5*p0*(cos(phi)-cos(p1))^2
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::HarmonicCosine, {1000.0, 35.0});
    double dt = cos_phi - std::cos(35.0 * DEG2RAD);
    double ref = 0.5 * K2E(1000.0) * dt * dt;
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // THREE_COSINE
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::ThreeCosine, {700.0, -200.0, 300.0});
    double ref = 0.5 * K2E(700.0) * (1.0 + cos_phi) + K2E(-200.0) * (1.0 - cos_phi2) +
                 0.5 * K2E(300.0) * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2);
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // SIX_COSINE (Ryckaert-Bellemans): sum_i p_i cos^i with alternating signs as in RASPA2
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::RyckaertBellemans,
                             {1204.654, 1947.740, -357.845, -1944.666, 715.690, -1565.572});
    double ref = K2E(1204.654) - K2E(1947.740) * cos_phi + K2E(-357.845) * cos_phi2 -
                 K2E(-1944.666) * cos_phi * cos_phi2 + K2E(715.690) * cos_phi2 * cos_phi2 -
                 K2E(-1565.572) * cos_phi2 * cos_phi2 * cos_phi;
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // TRAPPE: p0+p1*(1+cos(phi))+p2*(1-cos(2*phi))+p3*(1+cos(3*phi))
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32});
    double ref = K2E(0.0) + K2E(355.03) * (1.0 + std::cos(phi)) + K2E(-68.19) * (1.0 - std::cos(2.0 * phi)) +
                 K2E(791.32) * (1.0 + std::cos(3.0 * phi));
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, 1e-6 * std::abs(ref));
  }

  // TRAPPE_EXTENDED: p0+p1*cos(phi)+p2*cos(2*phi)+p3*cos(3*phi)+p4*cos(4*phi)
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::TraPPE_Extended, {100.0, 200.0, -300.0, 400.0, -50.0});
    double ref = K2E(100.0) + K2E(200.0) * std::cos(phi) + K2E(-300.0) * std::cos(2.0 * phi) +
                 K2E(400.0) * std::cos(3.0 * phi) + K2E(-50.0) * std::cos(4.0 * phi);
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, 1e-6 * std::abs(ref));
  }

  // MOD_TRAPPE: TraPPE with a phase shift p4
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::ModifiedTraPPE, {0.0, 355.03, -68.19, 791.32, 15.0});
    double shifted = phi - 15.0 * DEG2RAD;
    shifted -= std::rint(shifted / (2.0 * std::numbers::pi)) * 2.0 * std::numbers::pi;
    // Note: RASPA2's MOD_TRAPPE polynomial implements the 2nd harmonic as -p2*(1+cos(2*phi))
    // (i.e. -2*p2*cos^2(phi)), even though the comment in internal_energy.c writes p2*(1-cos(2*phi)).
    double ref = K2E(0.0) + K2E(355.03) * (1.0 + std::cos(shifted)) - K2E(-68.19) * (1.0 + std::cos(2.0 * shifted)) +
                 K2E(791.32) * (1.0 + std::cos(3.0 * shifted));
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, 1e-6 * std::abs(ref));
  }

  // CVFF: p0*(1+cos(p1*phi-p2))
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::CVFF, {1000.0, 3.0, 30.0});
    double ref = K2E(1000.0) * (1.0 + std::cos(3.0 * phi - 30.0 * DEG2RAD));
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // CFF: p0*(1-cos(phi))+p1*(1-cos(2*phi))+p2*(1-cos(3*phi))
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::CFF, {700.0, -200.0, 300.0});
    double ref = K2E(700.0) * (1.0 - std::cos(phi)) + K2E(-200.0) * (1.0 - std::cos(2.0 * phi)) +
                 K2E(300.0) * (1.0 - std::cos(3.0 * phi));
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, 1e-6 * std::abs(ref));
  }

  // CFF2: p0*(1+cos(phi))+p1*(1+cos(2*phi))+p2*(1+cos(3*phi))
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::CFF2, {700.0, -200.0, 300.0});
    double ref = K2E(700.0) * (1.0 + std::cos(phi)) + K2E(-200.0) * (1.0 + std::cos(2.0 * phi)) +
                 K2E(300.0) * (1.0 + std::cos(3.0 * phi));
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, 1e-6 * std::abs(ref));
  }

  // OPLS: (1/2)(p0+p1*(1+cos(phi))+p2*(1-cos(2*phi))+p3*(1+cos(3*phi)))
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::OPLS, {0.0, 710.06, -136.38, 1582.64});
    double ref = 0.5 * (K2E(0.0) + K2E(710.06) * (1.0 + std::cos(phi)) + K2E(-136.38) * (1.0 - std::cos(2.0 * phi)) +
                        K2E(1582.64) * (1.0 + std::cos(3.0 * phi)));
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, 1e-6 * std::abs(ref));
  }

  // MM3: three-cosine with parameters in kcal/mol
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::MM3, {0.2, 0.27, 0.093});
    double p0 = 0.2 * Units::KCalPerMolToEnergy;
    double p1 = 0.27 * Units::KCalPerMolToEnergy;
    double p2 = 0.093 * Units::KCalPerMolToEnergy;
    double ref = 0.5 * p0 * (1.0 + cos_phi) + p1 * (1.0 - cos_phi2) +
                 0.5 * p2 * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2);
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // FOURIER_SERIES
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::FourierSeries, {100.0, 200.0, -300.0, 400.0, -50.0, 60.0});
    // Note: RASPA2's FOURIER_SERIES polynomial implements the 6th term as (1-cos(6*phi)),
    // even though the comment in internal_energy.c writes (1+cos(6*phi)).
    double ref = 0.5 * (K2E(100.0) * (1.0 + std::cos(phi)) + K2E(200.0) * (1.0 - std::cos(2.0 * phi)) +
                        K2E(-300.0) * (1.0 + std::cos(3.0 * phi)) + K2E(400.0) * (1.0 - std::cos(4.0 * phi)) +
                        K2E(-50.0) * (1.0 + std::cos(5.0 * phi)) + K2E(60.0) * (1.0 - std::cos(6.0 * phi)));
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, 1e-6 * std::abs(ref));
  }

  // FOURIER_SERIES2
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::FourierSeries2, {100.0, 200.0, -300.0, 400.0, -50.0, 60.0});
    double ref = 0.5 * (K2E(100.0) * (1.0 + std::cos(phi)) + K2E(200.0) * (1.0 - std::cos(2.0 * phi)) +
                        K2E(-300.0) * (1.0 + std::cos(3.0 * phi)) + K2E(400.0) * (1.0 + std::cos(4.0 * phi)) +
                        K2E(-50.0) * (1.0 + std::cos(5.0 * phi)) + K2E(60.0) * (1.0 + std::cos(6.0 * phi)));
    EXPECT_NEAR(torsion.calculateEnergy(posA, posB, posC, posD), ref, 1e-6 * std::abs(ref));
  }

  // CVFF_BLOCKED: energy is zero by definition
  {
    TorsionPotential torsion({0, 1, 2, 3}, TorsionType::CVFFBlocked, {1.0, 1000.0, 2.0, 0.5, 1.5});
    EXPECT_EQ(torsion.calculateEnergy(posA, posB, posC, posD), 0.0);
  }
}

TEST(internal_potentials, bond_bond_forms_match_RASPA2)
{
  double r_ab = distance(posA, posB);
  double r_cb = distance(posC, posB);
  double tol = 1e-8;

  // CVFF/CFF: p0*(rab-p1)*(rbc-p2)
  {
    BondBondPotential bondBond({0, 1, 2}, BondBondType::CVFF, {10000.0, 1.1, 1.5});
    double ref = K2E(10000.0) * (r_ab - 1.1) * (r_cb - 1.5);
    EXPECT_NEAR(bondBond.calculateEnergy(posA, posB, posC), ref, tol * std::abs(ref));
  }
}

TEST(internal_potentials, bond_bend_forms_match_RASPA2)
{
  double r_ab = distance(posA, posB);
  double r_cb = distance(posC, posB);
  double theta = bendAngle(posA, posB, posC);
  double tol = 1e-8;

  // CVFF/CFF: (theta-p0)*(p1*(rab-p2)+p3*(rbc-p4))
  {
    BondBendPotential bondBend({0, 1, 2, 3}, BondBendType::CVFF, {110.0, 5000.0, 1.1, 6000.0, 1.5});
    double ref = (theta - 110.0 * DEG2RAD) * (K2E(5000.0) * (r_ab - 1.1) + K2E(6000.0) * (r_cb - 1.5));
    EXPECT_NEAR(bondBend.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // MM3: p0*((rab-p1)+(rbc-p2))*(theta-p3), with theta difference in degrees, p0 in mdyne/rad
  {
    BondBendPotential bondBend({0, 1, 2, 3}, BondBendType::MM3, {0.1, 1.1, 1.5, 110.0});
    double p0 = 0.1 * 2.51118 * Units::KCalPerMolToEnergy;
    double ref = p0 * ((r_ab - 1.1) + (r_cb - 1.5)) * RAD2DEG * (theta - 110.0 * DEG2RAD);
    EXPECT_NEAR(bondBend.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // TRUNCATED_HARMONIC
  {
    BondBendPotential bondBend({0, 1, 2, 3}, BondBendType::TruncatedHarmonic, {40000.0, 110.0, 2.0});
    double dt = theta - 110.0 * DEG2RAD;
    double ref = 0.5 * K2E(40000.0) * dt * dt *
                 std::exp(-(std::pow(r_ab, 8) + std::pow(r_cb, 8)) / std::pow(2.0, 8));
    EXPECT_NEAR(bondBend.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // SCREENED_HARMONIC
  {
    BondBendPotential bondBend({0, 1, 2, 3}, BondBendType::ScreenedHarmonic, {40000.0, 110.0, 2.0, 2.5});
    double dt = theta - 110.0 * DEG2RAD;
    double ref = 0.5 * K2E(40000.0) * dt * dt * std::exp(-(r_ab / 2.0 + r_cb / 2.5));
    EXPECT_NEAR(bondBend.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // SCREENED_VESSAL
  {
    BondBendPotential bondBend({0, 1, 2, 3}, BondBendType::ScreenedVessal, {40000.0, 110.0, 2.0, 2.5});
    double p1 = 110.0 * DEG2RAD;
    double t1 = (p1 - std::numbers::pi) * (p1 - std::numbers::pi) -
                (theta - std::numbers::pi) * (theta - std::numbers::pi);
    double ref = (K2E(40000.0) / (8.0 * (theta - std::numbers::pi) * (theta - std::numbers::pi))) * t1 * t1 *
                 std::exp(-(r_ab / 2.0 + r_cb / 2.5));
    EXPECT_NEAR(bondBend.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // TRUNCATED_VESSAL
  {
    BondBendPotential bondBend({0, 1, 2, 3}, BondBendType::TruncatedVessal, {40000.0, 110.0, 2.0, 2.5});
    double p1 = 110.0 * DEG2RAD;
    double p2 = 2.0;
    double ref = K2E(40000.0) *
                 (std::pow(theta, p2) * std::pow(theta - p1, 2) * std::pow(theta + p1 - 2.0 * std::numbers::pi, 2) -
                  0.5 * p2 * std::pow(std::numbers::pi, p2 - 1.0) * std::pow(theta - p1, 2) *
                      std::pow(std::numbers::pi - p1, 3.0)) *
                 std::exp(-(std::pow(r_ab, 8) + std::pow(r_cb, 8)) / std::pow(2.5, 8));
    EXPECT_NEAR(bondBend.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }
}

TEST(internal_potentials, bend_bend_forms_match_RASPA2)
{
  // Theta1 is the angle A-B-C, Theta2 is the angle A-B-D (atom B is the central atom)
  double theta1 = bendAngle(posA, posB, posC);
  double theta2 = bendAngle(posA, posB, posD);
  double tol = 1e-8;

  // CVFF/CFF: p0*(theta1-p1)*(theta2-p2)
  {
    BendBendPotential bendBend({0, 1, 2, 3}, BendBendType::CVFF, {10000.0, 110.0, 115.0});
    double ref = K2E(10000.0) * (theta1 - 110.0 * DEG2RAD) * (theta2 - 115.0 * DEG2RAD);
    EXPECT_NEAR(bendBend.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // MM3: -p0*(theta1-p1)*(theta2-p2), theta differences in degrees, p0 in mdyne A/rad^2
  {
    BendBendPotential bendBend({0, 1, 2, 3}, BendBendType::MM3, {0.24, 110.0, 115.0});
    double p0 = 0.24 * 0.02191418 * Units::KCalPerMolToEnergy;
    double ref = -p0 * (RAD2DEG * RAD2DEG) * (theta1 - 110.0 * DEG2RAD) * (theta2 - 115.0 * DEG2RAD);
    EXPECT_NEAR(bendBend.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }
}

TEST(internal_potentials, bond_torsion_forms_match_RASPA2)
{
  double r_bc = distance(posC, posB);
  double phi = dihedralAngle(posA, posB, posC, posD);
  double cos_phi = std::cos(phi);
  double cos_phi2 = cos_phi * cos_phi;
  double tol = 1e-8;

  // MM3: p0*(rbc-p3)*cos(phi)+p1*(rbc-p3)*cos(2phi)+p2*(rbc-p3)*cos(3phi), p_i in kcal/mol
  {
    BondTorsionPotential bondTorsion({0, 1, 2, 3}, BondTorsionType::MM3, {0.1, 0.2, 0.3, 1.5});
    double temp = r_bc - 1.5;
    double p0 = 0.1 * Units::KCalPerMolToEnergy;
    double p1 = 0.2 * Units::KCalPerMolToEnergy;
    double p2 = 0.3 * Units::KCalPerMolToEnergy;
    double ref = p0 * temp * cos_phi + p1 * temp * (2.0 * cos_phi2 - 1.0) +
                 p2 * temp * (4.0 * cos_phi2 * cos_phi - 3.0 * cos_phi);
    EXPECT_NEAR(bondTorsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }
}

TEST(internal_potentials, bend_torsion_forms_match_RASPA2)
{
  // Theta1 is the bend angle A-B-C, Theta2 the bend angle B-C-D, phi the dihedral.
  double3 Dab = posA - posB;
  double r_ab = std::sqrt(double3::dot(Dab, Dab));
  double3 Dbc = (posC - posB).normalized();
  double3 Dcd = posD - posC;
  double r_cd = std::sqrt(double3::dot(Dcd, Dcd));

  double theta1 = std::acos(std::clamp(double3::dot(Dab, Dbc) / r_ab, -1.0, 1.0));
  double theta2 = std::acos(std::clamp(-double3::dot(Dcd, Dbc) / r_cd, -1.0, 1.0));
  double phi = dihedralAngle(posA, posB, posC, posD);
  double cos_phi = std::cos(phi);
  double cos_phi2 = cos_phi * cos_phi;
  double tol = 1e-8;

  // CVFF/CFF cross: p0*(theta1-p1)*(theta2-p2)*cos(phi)
  {
    BendTorsionPotential bendTorsion({0, 1, 2, 3}, BendTorsionType::CVFF, {10000.0, 110.0, 115.0});
    double ref = K2E(10000.0) * (theta1 - 110.0 * DEG2RAD) * (theta2 - 115.0 * DEG2RAD) * cos_phi;
    EXPECT_NEAR(bendTorsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // SMOOTHED: S(theta1)*p0*(1+cos(p1*phi-p2))*S(theta2)
  {
    BendTorsionPotential bendTorsion({0, 1, 2, 3}, BendTorsionType::Smoothed, {1000.0, 3.0, 30.0});
    double ref = K2E(1000.0) * (1.0 + std::cos(3.0 * phi - 30.0 * DEG2RAD)) * smoothing(theta1) * smoothing(theta2);
    EXPECT_NEAR(bendTorsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // SMOOTHED_THREE_COSINE
  {
    BendTorsionPotential bendTorsion({0, 1, 2, 3}, BendTorsionType::SmoothedThreeCosine, {700.0, -200.0, 300.0});
    double ref = (0.5 * K2E(700.0) * (1.0 + cos_phi) + K2E(-200.0) * (1.0 - cos_phi2) +
                  0.5 * K2E(300.0) * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2)) *
                 smoothing(theta1) * smoothing(theta2);
    EXPECT_NEAR(bendTorsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // NICHOLAS: like SMOOTHED_THREE_COSINE, but smoothed on theta1 only
  {
    BendTorsionPotential bendTorsion({0, 1, 2, 3}, BendTorsionType::Nicholas, {700.0, -200.0, 300.0});
    double ref = (0.5 * K2E(700.0) * (1.0 + cos_phi) + K2E(-200.0) * (1.0 - cos_phi2) +
                  0.5 * K2E(300.0) * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2)) *
                 smoothing(theta1);
    EXPECT_NEAR(bendTorsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }

  // SMOOTHED_CFF
  {
    BendTorsionPotential bendTorsion({0, 1, 2, 3}, BendTorsionType::SmoothedCFF, {700.0, -200.0, 300.0});
    double ref = (K2E(700.0) * (1.0 - std::cos(phi)) + K2E(-200.0) * (1.0 - std::cos(2.0 * phi)) +
                  K2E(300.0) * (1.0 - std::cos(3.0 * phi))) *
                 smoothing(theta1) * smoothing(theta2);
    EXPECT_NEAR(bendTorsion.calculateEnergy(posA, posB, posC, posD), ref, 1e-6 * std::abs(ref));
  }

  // SMOOTHED_CFF2
  {
    BendTorsionPotential bendTorsion({0, 1, 2, 3}, BendTorsionType::SmoothedCFF2, {700.0, -200.0, 300.0});
    double ref = (K2E(700.0) * (1.0 + std::cos(phi)) + K2E(-200.0) * (1.0 + std::cos(2.0 * phi)) +
                  K2E(300.0) * (1.0 + std::cos(3.0 * phi))) *
                 smoothing(theta1) * smoothing(theta2);
    EXPECT_NEAR(bendTorsion.calculateEnergy(posA, posB, posC, posD), ref, 1e-6 * std::abs(ref));
  }

  // SMOOTHED_CFF_BEND_TORSION_CROSS
  {
    BendTorsionPotential bendTorsion({0, 1, 2, 3}, BendTorsionType::SmoothedCFF3, {10000.0, 110.0, 115.0});
    double ref = K2E(10000.0) * (theta1 - 110.0 * DEG2RAD) * (theta2 - 115.0 * DEG2RAD) * cos_phi *
                 smoothing(theta1) * smoothing(theta2);
    EXPECT_NEAR(bendTorsion.calculateEnergy(posA, posB, posC, posD), ref, tol * std::abs(ref));
  }
}
