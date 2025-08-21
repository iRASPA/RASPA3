module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <fstream>
#include <map>
#include <numbers>
#include <print>
#include <source_location>
#include <tuple>
#include <utility>
#include <vector>
#endif

module torsion_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import double3;

TorsionPotential::TorsionPotential(std::array<std::size_t, 4> identifiers, TorsionType type,
                                   std::vector<double> vector_parameters)
    : identifiers(identifiers), type(type)
{
  for (std::size_t i = 0; i < std::min(vector_parameters.size(), maximumNumberOfTorsionParameters); ++i)
  {
    parameters[i] = vector_parameters[i];
  }
  switch (type)
  {
    case TorsionType::Fixed:
      parameters[0] *= Units::DegreesToRadians;
      break;
    case TorsionType::Harmonic:
      // (1/2)*p_0*(phi-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    case TorsionType::HarmonicCosine:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] = std::cos(parameters[1] * Units::DegreesToRadians);
      break;
    case TorsionType::ThreeCosine:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      break;
    case TorsionType::RyckaertBellemans:
      // Prod_i=0^5 p_i*cos(phi)^i
      // =========================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      parameters[4] *= Units::KelvinToEnergy;
      parameters[5] *= Units::KelvinToEnergy;
      break;
    case TorsionType::TraPPE:
      // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
      // =============================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      break;
    case TorsionType::TraPPE_Extended:
      // p_0[0]+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
      // ================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      parameters[4] *= Units::KelvinToEnergy;
      break;
    case TorsionType::ModifiedTraPPE:
      // p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
      // ==========================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      parameters[4] *= Units::KelvinToEnergy;
      parameters[5] *= Units::DegreesToRadians;
      break;
    case TorsionType::CVFF:
      // p_0*(1+cos(p_1*phi-p_2))
      // ========================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::DegreesToRadians;
      break;
    case TorsionType::CFF:
      // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      break;
    case TorsionType::CFF2:
      // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      break;
    case TorsionType::OPLS:
      // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
      // =================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      break;
    case TorsionType::MM3:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0     [kcal/mol]
      // p_1     [kcal/mol]
      // p_2     [kcal/mol]
      parameters[0] *= Units::KCalPerMolToEnergy;
      parameters[1] *= Units::KCalPerMolToEnergy;
      parameters[2] *= Units::KCalPerMolToEnergy;
      break;
    case TorsionType::FourierSeries:
      // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
      // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
      // =======================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      parameters[4] *= Units::KelvinToEnergy;
      parameters[5] *= Units::KelvinToEnergy;
      break;
    case TorsionType::FourierSeries2:
      // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
      // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
      // =======================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      parameters[0] *= Units::KelvinToEnergy;
      parameters[1] *= Units::KelvinToEnergy;
      parameters[2] *= Units::KelvinToEnergy;
      parameters[3] *= Units::KelvinToEnergy;
      parameters[4] *= Units::KelvinToEnergy;
      parameters[5] *= Units::KelvinToEnergy;
      break;
  }
}

std::string TorsionPotential::print() const
{
  switch (type)
  {
    case TorsionType::Fixed:
      return std::format("{} - {} - {} - {} : FIXED\n", identifiers[0], identifiers[1], identifiers[2], identifiers[3]);
    case TorsionType::Harmonic:
      // (1/2)*p_0*(phi-p_1)^2
      // ===============================================
      // p_0/k_B [K/rad^2]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : HARMONIC p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n", identifiers[0],
                         identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::RadiansToDegrees);
    case TorsionType::HarmonicCosine:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      return std::format("{} - {} - {} - {} : HARMONIC_COSINE p_0/k_B={:g} [K/rad^2], p_1={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, parameters[1] * Units::RadiansToDegrees);
    case TorsionType::ThreeCosine:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      return std::format("{} - {} - {} - {} : THREE_COSINE p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, parameters[1] * Units::EnergyToKelvin,
                         parameters[2] * Units::EnergyToKelvin);
    case TorsionType::RyckaertBellemans:
      // Prod_i=0^5 p_i*cos(phi)^i
      // =========================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      return std::format(
          "{} - {} - {} - {} : RYCKAERT_BELLEMANS p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K], p_3/k_B={:g} "
          "[K], p_4/k_B={:g} [K], p_5/k_B={:g} [K]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::EnergyToKelvin, parameters[2] * Units::EnergyToKelvin,
          parameters[3] * Units::EnergyToKelvin, parameters[4] * Units::EnergyToKelvin,
          parameters[5] * Units::EnergyToKelvin);
    case TorsionType::TraPPE:
      // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
      // =============================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      return std::format(
          "{} - {} - {} - {} : TRAPPE p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K], p_3/k_B={:g} [K]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::EnergyToKelvin, parameters[2] * Units::EnergyToKelvin,
          parameters[3] * Units::EnergyToKelvin);
    case TorsionType::TraPPE_Extended:
      // p_0[0]+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
      // ================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      return std::format(
          "{} - {} - {} - {} : TRAPPE_EXTENDED p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K], p_3/k_B={:g} [K], "
          "p_4/k_B={:g} [K]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::EnergyToKelvin, parameters[2] * Units::EnergyToKelvin,
          parameters[3] * Units::EnergyToKelvin, parameters[4] * Units::EnergyToKelvin);
    case TorsionType::ModifiedTraPPE:
      // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
      // =============================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4     [degrees]
      return std::format(
          "{} - {} - {} - {} : TRAPPE_MODIFIED p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K], p_3/k_B={:g} [K], "
          "p_4={:g} [degrees]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::EnergyToKelvin, parameters[2] * Units::EnergyToKelvin,
          parameters[3] * Units::EnergyToKelvin, parameters[4] * Units::RadiansToDegrees);
    case TorsionType::CVFF:
      // p_0*(1+cos(p_1*phi-p_2))
      // ========================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      return std::format("{} - {} - {} - {} : CVFF p_0/k_B={:g} [K], p_1={:g} [-], p_2={:g} [degrees]\n",
                         identifiers[0], identifiers[1], identifiers[2], identifiers[3],
                         parameters[0] * Units::EnergyToKelvin, parameters[1], parameters[2] * Units::RadiansToDegrees);
    case TorsionType::CFF:
      // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      return std::format("{} - {} - {} - {} : CFF p_0/k_B={:g} [K], p_1={:g} [-], p_2/k_B={:g} [K]\n", identifiers[0],
                         identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::EnergyToKelvin, parameters[2] * Units::EnergyToKelvin);
    case TorsionType::CFF2:
      // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      return std::format("{} - {} - {} - {} : CFF2 p_0/k_B={:g} [K], p_1={:g} [-], p_2/k_B={:g} [K]\n", identifiers[0],
                         identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
                         parameters[1] * Units::EnergyToKelvin, parameters[2] * Units::EnergyToKelvin);
    case TorsionType::OPLS:
      // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
      // =================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      return std::format(
          "{} - {} - {} - {} : OPLS p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K], p_3/k_B={:g} [K]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::EnergyToKelvin, parameters[2] * Units::EnergyToKelvin,
          parameters[3] * Units::EnergyToKelvin);
    case TorsionType::MM3:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0     [kcal/mol]
      // p_1     [kcal/mol]
      // p_2     [kcal/mol]
      return std::format(
          "{} - {} - {} - {} : MM3 p_0/k_B={:g} [Kcal/mol], p_1/k_B={:g} [Kcal/mol], p_2/k_B={:g} [Kcal/mol]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKCalPerMol,
          parameters[1] * Units::EnergyToKCalPerMol, parameters[2] * Units::EnergyToKCalPerMol);
    case TorsionType::FourierSeries:
      // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
      // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
      // =======================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      return std::format(
          "{} - {} - {} - {} : FOURIER_SERIES p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K], p_3/k_B={:g} [K], "
          "p_4/k_B={:g} [K], p_5/k_B={:g} [K]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::EnergyToKelvin, parameters[2] * Units::EnergyToKelvin,
          parameters[3] * Units::EnergyToKelvin, parameters[4] * Units::EnergyToKelvin,
          parameters[5] * Units::EnergyToKelvin);
    case TorsionType::FourierSeries2:
      // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
      // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
      // =======================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      return std::format(
          "{} - {} - {} - {} : FOURIER_SERIES2 p_0/k_B={:g} [K], p_1/k_B={:g} [K], p_2/k_B={:g} [K], p_3/k_B={:g} [K], "
          "p_4/k_B={:g} [K], p_5/k_B={:g} [K]\n",
          identifiers[0], identifiers[1], identifiers[2], identifiers[3], parameters[0] * Units::EnergyToKelvin,
          parameters[1] * Units::EnergyToKelvin, parameters[2] * Units::EnergyToKelvin,
          parameters[3] * Units::EnergyToKelvin, parameters[4] * Units::EnergyToKelvin,
          parameters[5] * Units::EnergyToKelvin);
  }
}

double TorsionPotential::calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posC,
                                         const double3 &posD) const
{
  double cos_phi, cos_phi2, phi;
  double temp, temp2, sign;
  double shifted_cos_phi, shifted_cos_phi2;

  double3 Dab = posA - posB;

  double3 Dcb = (posC - posB).normalized();

  double3 Ddc = posD - posC;

  double dot_ab = double3::dot(Dab, Dcb);
  double dot_dc = double3::dot(Ddc, Dcb);

  double3 dr = (Dab - dot_ab * Dcb).normalized();
  double3 ds = (Ddc - dot_dc * Dcb).normalized();

  // compute Cos(Phi)
  // Phi is defined in protein convention Phi(trans)=Pi
  cos_phi = double3::dot(dr, ds);

  // Ensure CosPhi is between -1 and 1.
  cos_phi = std::clamp(cos_phi, -1.0, 1.0);
  cos_phi2 = cos_phi * cos_phi;

  switch (type)
  {
    case TorsionType::Fixed:
      return 0.0;
    case TorsionType::Harmonic:
      // ========================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      // potential defined in terms of 'phi' and therefore contains a singularity
      // the sign of the angle-phi is positive if (Rab x Rcb) x (Rcb x Rdc) is in the
      // same direction as Rbc, and negative otherwise
      sign = double3::dot(Dcb, double3::cross(double3::cross(Dab, Dcb), double3::cross(Dcb, Ddc)));
      phi = std::copysign(std::acos(cos_phi), sign);
      return parameters[0] * (1.0 + std::cos(parameters[1] * phi - parameters[2]));
    case TorsionType::HarmonicCosine:
      // (1/2)*p_0*(cos(phi)-cos(p_1))^2
      // ===============================================
      // p_0/k_B [K]
      // p_1     [degrees]
      temp = cos_phi - parameters[1];
      temp2 = temp * temp;
      return 0.5 * parameters[0] * temp2;
    case TorsionType::ThreeCosine:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      return 0.5 * parameters[0] * (1.0 + cos_phi) + parameters[1] * (1.0 - cos_phi2) +
             0.5 * parameters[2] * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2);
    case TorsionType::RyckaertBellemans:
      // Prod_i=0^5 p_i*cos(phi)^i
      // =========================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      // the Ryckaert-Bellemans potentials is often used for alkanes, the use implies exclusion of VDW-interactions
      // between the first and last atoms of the dihedral, and Phi'=Phi-Pi is defined accoording to the
      // polymer convention Phi'(trans)=0.
      return parameters[0] - parameters[1] * cos_phi + parameters[2] * cos_phi2 - parameters[3] * cos_phi * cos_phi2 +
             parameters[4] * cos_phi2 * cos_phi2 - parameters[5] * cos_phi2 * cos_phi2 * cos_phi;
    case TorsionType::TraPPE:
      // p_0+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
      // ==========================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      return parameters[0] +
             (1.0 + cos_phi) * (parameters[1] + parameters[3] -
                                2.0 * (cos_phi - 1.0) * (parameters[2] - 2.0 * parameters[3] * cos_phi));
    case TorsionType::TraPPE_Extended:
      // p_0+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
      // =============================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      return parameters[0] - parameters[2] + parameters[4] + (parameters[1] - 3.0 * parameters[3]) * cos_phi +
             (2.0 * parameters[2] - 8.0 * parameters[4]) * cos_phi2 + 4.0 * parameters[3] * cos_phi2 * cos_phi +
             8.0 * parameters[4] * cos_phi2 * cos_phi2;
    case TorsionType::ModifiedTraPPE:
      /* Salvador modification: 16/08/2016
      add phase in cos function:
      p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
     */
      sign = double3::dot(Dcb, double3::cross(double3::cross(Dab, Dcb), double3::cross(Dcb, Ddc)));
      phi = std::copysign(std::acos(cos_phi), sign);
      phi -= parameters[4];  // shift Phi as Phi+parameters[4]
      phi -= std::rint(phi / (2.0 * std::numbers::pi)) * 2.0 * std::numbers::pi;
      shifted_cos_phi = std::cos(phi);
      shifted_cos_phi2 = shifted_cos_phi;
      return parameters[0] + parameters[1] + parameters[3] + (parameters[1] - 3.0 * parameters[3]) * shifted_cos_phi -
             2.0 * parameters[2] * shifted_cos_phi2 + 4.0 * parameters[3] * shifted_cos_phi * shifted_cos_phi2;
    case TorsionType::CVFF:
      // p_0*(1+cos(p_1*phi-p_2))
      // ========================
      // p_0/k_B [K]
      // p_1     [-]
      // p_2     [degrees]
      // potential defined in terms of 'phi' and therefore contains a singularity
      // the sign of the angle-phi is positive if (Rab x Rbc) x (Rbc x Rcd) is in the
      // same direction as Rbc, and negative otherwise
      sign = double3::dot(Dcb, double3::cross(double3::cross(Dab, Dcb), double3::cross(Dcb, Ddc)));
      phi = std::copysign(std::acos(cos_phi), sign);
      return parameters[0] * (1.0 + std::cos(parameters[1] * phi - parameters[2]));
    case TorsionType::CFF:
      // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      return parameters[0] * (1.0 - cos_phi) + 2.0 * parameters[1] * (1.0 - cos_phi2) +
             parameters[2] * (1.0 + 3.0 * cos_phi - 4.0 * cos_phi * cos_phi2);
    case TorsionType::CFF2:
      // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
      // ======================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      return parameters[0] * (1.0 + cos_phi) + parameters[2] +
             cos_phi * (-3.0 * parameters[2] + 2.0 * cos_phi * (parameters[1] + 2.0 * parameters[2] * cos_phi));
    case TorsionType::OPLS:
      // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
      // =================================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      return 0.5 * (parameters[0] +
                    (1.0 + cos_phi) * (parameters[1] + parameters[3] -
                                       2.0 * (cos_phi - 1.0) * (parameters[2] - 2.0 * parameters[3] * cos_phi)));
    case TorsionType::MM3:
      // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
      // ========================================================================
      // p_0     [kcal/mol]
      // p_1     [kcal/mol]
      // p_2     [kcal/mol]
      return 0.5 * parameters[0] * (1.0 + cos_phi) + parameters[1] * (1.0 - cos_phi2) +
             0.5 * parameters[2] * (1.0 - 3.0 * cos_phi + 4.0 * cos_phi * cos_phi2);
    case TorsionType::FourierSeries:
      // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
      // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
      // =======================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      return 0.5 * (parameters[0] + 2.0 * parameters[1] + parameters[2] + parameters[4] + 2.0 * parameters[5] +
                    (parameters[0] - 3.0 * parameters[2] + 5.0 * parameters[4]) * cos_phi -
                    2.0 * (parameters[1] - 4.0 * parameters[3] + 9.0 * parameters[5]) * cos_phi2 +
                    4.0 * (parameters[2] - 5.0 * parameters[4]) * cos_phi2 * cos_phi -
                    8.0 * (parameters[3] - 6.0 * parameters[5]) * cos_phi2 * cos_phi2 +
                    16.0 * parameters[4] * cos_phi2 * cos_phi2 * cos_phi -
                    32.0 * parameters[5] * cos_phi2 * cos_phi2 * cos_phi2);
    case TorsionType::FourierSeries2:
      // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
      // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
      // =======================================================================
      // p_0/k_B [K]
      // p_1/k_B [K]
      // p_2/k_B [K]
      // p_3/k_B [K]
      // p_4/k_B [K]
      // p_5/k_B [K]
      return 0.5 * (parameters[2] + 2.0 * parameters[3] + parameters[4] - 3.0 * parameters[2] * cos_phi +
                    5.0 * parameters[4] * cos_phi + parameters[0] * (1.0 + cos_phi) +
                    2.0 * (parameters[1] - parameters[1] * cos_phi2 +
                           cos_phi2 * (parameters[5] * (3.0 - 4.0 * cos_phi2) * (3.0 - 4.0 * cos_phi2) +
                                       4.0 * parameters[3] * (cos_phi2 - 1.0) +
                                       2.0 * cos_phi * (parameters[2] + parameters[4] * (4.0 * cos_phi2 - 5.0)))));
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const TorsionPotential &b)
{
  archive << b.versionNumber;

  archive << b.type;
  archive << b.identifiers;
  archive << b.parameters;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, TorsionPotential &b)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > b.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'TorsionPotential' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> b.type;
  archive >> b.identifiers;
  archive >> b.parameters;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("TorsionPotential: Error in binary restart\n"));
  }
#endif

  return archive;
}
