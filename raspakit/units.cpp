module;

module units;

import <string>;
import <sstream>;
import <ostream>;
import <print>;

import stringutils;

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
  std::print(stream, "Unit of acceleration:        {} [m^2/s]\n", Units::AccelerationConversionFactor);
  std::print(stream, "Unit of diffusion:           {} [m^2/s]\n", Units::DiffusionConversionFactor);
  std::print(stream, "Unit of dipole moment:       {} [C.m]\n", Units::DipoleMomentConversionFactor);
  std::print(stream, "Unit of electric potential:  {} [V]\n", Units::ElectricPotentialConversionFactor);
  std::print(stream, "Unit of electric field:      {} [V]\n", Units::ElectricFieldConversionFactor);
  std::print(stream, "Unit of polarizability:      {} [-]\n", Units::PolarizilibityConversionFactor);
  std::print(stream, "Unit of Coulomb potential:   {} [K]\n", Units::CoulombicConversionFactor * Units::EnergyToKelvin);
  std::print(stream, "Unit of dielectric constant: {} [s^2 C^2/(kg m^3)]\n", Units::DielectricConstantConversionFactor);
  std::print(stream, "\n\n");

  std::print(stream, "Derived units and their conversion factors:\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Boltzmann constant:          {} [-]\n", Units::KB);
  std::print(stream, "Unit of energy:              {} [J]\n", Units::EnergyConversionFactor);
  std::print(stream, "Unit of force:               {} [N]\n", Units::ForceConversionFactor);
  std::print(stream, "Unit of pressure:            {} [Pa]\n", Units::PressureConversionFactor);
  std::print(stream, "Unit of velocity:            {} [m/s]\n", Units::VelocityConversionFactor);
  std::print(stream, "Unit of acceleration:        {} [m^2/s]\n", Units::AccelerationConversionFactor);
  std::print(stream, "Unit of diffusion:           {} [m^2/s]\n", Units::DiffusionConversionFactor);
  std::print(stream, "Unit of dipole moment:       {} [C.m]\n", Units::DipoleMomentConversionFactor);
  std::print(stream, "Unit of electric potential:  {} [V]\n", Units::ElectricPotentialConversionFactor);
  std::print(stream, "Unit of electric field:      {} [V]\n", Units::ElectricFieldConversionFactor);
  std::print(stream, "Unit of polarizability:      {} [-]\n", Units::PolarizilibityConversionFactor);
  std::print(stream, "Unit of Coulomb potential:   {} [K]\n", Units::CoulombicConversionFactor * Units::EnergyToKelvin);
  std::print(stream, "Unit of dielectric constant: {} [s^2 C^2/(kg m^3)]\n", Units::DielectricConstantConversionFactor);
  std::print(stream, "\n\n");

  std::print(stream, "Internal conversion factors:\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Energy to Kelvin:      {} [-]\n", Units::EnergyToKelvin);
  std::print(stream, "Kelvin to energy:      {} [-]\n", Units::KelvinToEnergy);
  std::print(stream, "\n\n");

  return stream.str();
}
