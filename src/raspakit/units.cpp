module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <sstream>
#include <ostream>
#include <print>
#endif

module units;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <sstream>;
import <ostream>;
import <print>;
#endif


import stringutils;
import json;

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

nlohmann::json Units::jsonStatus()
{
  nlohmann::json units;

  units["temperature"] = "Kelvin";
  units["length [m]"] = Units::LengthUnit;
  units["time [s]"] = Units::TimeUnit;
  units["mass [kg]"] = Units::MassUnit;
  units["charge [C/particle]"] = Units::ChargeUnit;
  units["Boltzmann constant [-]"] = Units::KB;
  units["energy [J]"] = Units::EnergyConversionFactor;
  units["force [N]"] = Units::ForceConversionFactor;
  units["pressure [Pa]"] = Units::PressureConversionFactor;
  units["velocity [m/s]"] = Units::VelocityConversionFactor;
  units["acceleration [m²/s]"] = Units::AccelerationConversionFactor;
  units["diffusion [m²/s]"] = Units::DiffusionConversionFactor;
  units["dipole moment [C.m]"] = Units::DipoleMomentConversionFactor;
  units["electric potential [V]"] = Units::ElectricPotentialConversionFactor;
  units["electric field [V]"] = Units::ElectricFieldConversionFactor;
  units["polarizability [-]"] = Units::PolarizilibityConversionFactor;
  units["dielectric constant [s² C²/(kg m³)]"] = Units::DielectricConstantConversionFactor;
  units["Boltzmann constant [-]"] = Units::KB;
  units["energy [J]"] = Units::EnergyConversionFactor;
  units["force [N]"] = Units::ForceConversionFactor;
  units["pressure [Pa]"] = Units::PressureConversionFactor;
  units["velocity [m/s]"] = Units::VelocityConversionFactor;
  units["Coulomb potential [K]"] = Units::CoulombicConversionFactor * Units::EnergyToKelvin;
  units["Energy to Kelvin [-]"] = Units::EnergyToKelvin;
  units["Kelvin to energy [-]"] = Units::KelvinToEnergy;

  return units;
}
