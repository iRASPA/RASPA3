module;

#ifdef USE_LEGACY_HEADERS
#include <ostream>
#include <sstream>
#include <string>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#endif

module units;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <sstream>;
import <ostream>;
#if defined(__has_include) && __has_include(<print>)
import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
import print;
#endif

import stringutils;
import hdf5;

std::string Units::printStatus()
{
  std::ostringstream stream;

  std::print(stream, "Mutual consistent basic set of units:\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Unit of temperature: {}\n", "Kelvin");
  std::print(stream, "Unit of length:      {} [m]\n", Units::LengthUnit);
  std::print(stream, "Unit of time:        {} [s]\n", Units::TimeUnit);
  std::print(stream, "Unit of mass:        {} [kg]\n", Units::MassUnit);
  std::print(stream, "Unit of charge:      {} [C/particle]\n", Units::ChargeUnit);
  std::print(stream, "\n\n");

  std::print(stream, "Derived units and their conversion factors:\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Boltzmann constant:          {} [-]\n", Units::KB);
  std::print(stream, "Unit of energy:              {} [J]\n", Units::EnergyConversionFactor);
  std::print(stream, "Unit of force:               {} [N]\n", Units::ForceConversionFactor);
  std::print(stream, "Unit of pressure:            {} [Pa]\n", Units::PressureConversionFactor);
  std::print(stream, "Unit of velocity:            {} [m/s]\n", Units::VelocityConversionFactor);
  std::print(stream, "Unit of acceleration:        {} [m²/s]\n", Units::AccelerationConversionFactor);
  std::print(stream, "Unit of diffusion:           {} [m²/s]\n", Units::DiffusionConversionFactor);
  std::print(stream, "Unit of dipole moment:       {} [C.m]\n", Units::DipoleMomentConversionFactor);
  std::print(stream, "Unit of electric potential:  {} [V]\n", Units::ElectricPotentialConversionFactor);
  std::print(stream, "Unit of electric field:      {} [V]\n", Units::ElectricFieldConversionFactor);
  std::print(stream, "Unit of polarizability:      {} [-]\n", Units::PolarizilibityConversionFactor);
  std::print(stream, "Unit of Coulomb potential:   {} [K]\n", Units::CoulombicConversionFactor * Units::EnergyToKelvin);
  std::print(stream, "Unit of dielectric constant: {} [s² C²/(kg m³)]\n", Units::DielectricConstantConversionFactor);
  std::print(stream, "\n\n");

  std::print(stream, "Internal conversion factors:\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Energy to Kelvin:      {} [-]\n", Units::EnergyToKelvin);
  std::print(stream, "Kelvin to energy:      {} [-]\n", Units::KelvinToEnergy);
  std::print(stream, "\n\n");

  return stream.str();
}

void Units::logStatus(HDF5Writer& hdf5)
{
  hdf5.createGroup("units");
  hdf5.writeMetaInfo("units", "temperature", "Kelvin");
  hdf5.writeMetaInfo("units", "length [m]", Units::LengthUnit);
  hdf5.writeMetaInfo("units", "time [s]", Units::TimeUnit);
  hdf5.writeMetaInfo("units", "mass [kg]", Units::MassUnit);
  hdf5.writeMetaInfo("units", "charge [C/particle]", Units::ChargeUnit);
  hdf5.writeMetaInfo("units", "Boltzmann constant [-]", Units::KB);
  hdf5.writeMetaInfo("units", "energy [J]", Units::EnergyConversionFactor);
  hdf5.writeMetaInfo("units", "force [N]", Units::ForceConversionFactor);
  hdf5.writeMetaInfo("units", "pressure [Pa]", Units::PressureConversionFactor);
  hdf5.writeMetaInfo("units", "velocity [m/s]", Units::VelocityConversionFactor);
  hdf5.writeMetaInfo("units", "acceleration [m²/s]", Units::AccelerationConversionFactor);
  hdf5.writeMetaInfo("units", "diffusion [m²/s]", Units::DiffusionConversionFactor);
  hdf5.writeMetaInfo("units", "dipole moment [C.m]", Units::DipoleMomentConversionFactor);
  hdf5.writeMetaInfo("units", "electric potential [V]", Units::ElectricPotentialConversionFactor);
  hdf5.writeMetaInfo("units", "electric field [V]", Units::ElectricFieldConversionFactor);
  hdf5.writeMetaInfo("units", "polarizability [-]", Units::PolarizilibityConversionFactor);
  hdf5.writeMetaInfo("units", "Coulomb potential [K]", Units::CoulombicConversionFactor * Units::EnergyToKelvin);
  hdf5.writeMetaInfo("units", "Unit of dielectric constant:[s² C²/(kg m³)]", Units::DielectricConstantConversionFactor);
  hdf5.writeMetaInfo("units", "Energy to Kelvin [-]", Units::EnergyToKelvin);
  hdf5.writeMetaInfo("units", "Kelvin to energy [-]", Units::KelvinToEnergy);
}