module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <ostream>
#include <print>
#include <sstream>
#include <string>
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

  switch (unitSystem)
  {
    case System::RASPA:
      std::print(stream, "Mutual consistent basic set of units:\n");
      std::print(stream, "===============================================================================\n\n");
      std::print(stream, "Unit of temperature: {} [{}]\n", 1, unitOfTemperatureString);
      std::print(stream, "Unit of length:      {} [{}]\n", Units::LengthUnit, unitOfLengthString);
      std::print(stream, "Unit of time:        {} [{}]\n", Units::TimeUnit, unitOfTimeString);
      std::print(stream, "Unit of mass:        {} [{}]\n", Units::MassUnit, unitOfMassString);
      std::print(stream, "Unit of charge:      {} [{}]\n", Units::ChargeUnit, unitOfChargeString);
      std::print(stream, "\n\n");

      std::print(stream, "Boltzmann constant (internal units):          {} [-]\n", Units::KB);
      std::print(stream, "represents conversion from energy (internal units) to temperature (in Kelvin)\n\n");

      std::print(stream, "Derived units and their conversion factors:\n");
      std::print(stream, "===============================================================================\n\n");
      std::print(stream, "Unit of energy:                           {} [{}]\n", Units::EnergyConversionFactor,
                 unitOfEnergyString);
      std::print(stream, "Unit of energy:                           {} [{}]\n",
                 Units::EnergyConversionFactor * AvogadroConstant, unitOfEnergyPerMolString);
      std::print(stream, "Unit of velocity:                         {} [{}]\n", Units::VelocityConversionFactor,
                 unitOfVelocityString);
      std::print(stream, "Unit of force:                            {} [{}]\n", Units::ForceConversionFactor,
                 unitOfForceString);
      std::print(stream, "Unit of diffusion:                        {} [{}]\n", Units::DiffusionConversionFactor,
                 unitOfDiffusionString);
      std::print(stream, "Unit of acceleration:                     {} [{}]\n", Units::AccelerationConversionFactor,
                 unitOfAccelerationString);
      std::print(stream, "Unit of torque:                           {} [{}]\n", Units::TorqueConversionFactor,
                 unitOfTorqueString);
      std::print(stream, "Unit of pressure:                         {} [{}]\n", Units::PressureConversionFactor,
                 unitOfPressureString);
      std::print(stream, "Unit of volume:                           {} [{}]\n", Units::VolumeConversionFactor,
                 unitOfVolumeString);
      std::print(stream, "Unit of density:                          {} [{}]\n", Units::DensityConversionFactor,
                 unitOfDensityString);
      std::print(stream, "Unit of dynamic viscosity:                {} [{}]\n", Units::DynamicViscosityConversionFactor,
                 unitOfDynamicViscosityString);
      std::print(stream, "Unit of enthalpy:                         {} [{}]\n", Units::EnthalpyConversionFactor,
                 unitOfEnthalpyString);
      std::print(stream, "Unit of polarizability:                   {} [{}]\n", Units::PolarizilibityConversionFactor,
                 unitOfPolarizabilityString);
      std::print(stream, "Unit of Coulomb potential:                {} [{}]\n",
                 Units::CoulombicConversionFactor * Units::EnergyToKelvin, unitOfCoulombPotentialString);
      std::print(stream, "Unit of dielectric constant:              {} [{}]\n",
                 Units::DielectricConstantConversionFactor, unitOfDielectricConstantString);
      std::print(stream, "Unit of dipole moment:                    {} [{}]\n", Units::DipoleMomentConversionFactor,
                 unitOfDipoleMomentString);
      std::print(stream, "Unit of Debye:                            {} [{}]\n", Units::DebyeConversionFactor,
                 unitOfDebyeString);
      std::print(stream, "Unit of electric potential:               {} [{}]\n",
                 Units::ElectricPotentialConversionFactor, unitOfElectricPotentialString);
      std::print(stream, "Unit of electric field:                   {} [{}]\n", Units::ElectricFieldConversionFactor,
                 unitOfElectricFieldString);
      std::print(stream, "Unit of isothermal compressibility:       {} [{}]\n",
                 Units::IsothermalCompressibilityConversionFactor, unitOfIsothermalCompressibilityString);
      std::print(stream, "Unit of heat capacity:                    {} [{}]\n", Units::HeatCapacityConversionFactor,
                 unitOfHeatCapacityString);
      std::print(stream, "Unit of volumetric expansion coefficient: {} [{}]\n",
                 Units::VolumetricExpansionCoefficientConversionFactor, unitOfVolumetricExpansionCoefficientString);
      std::print(stream, "\n\n");

      std::print(stream, "Internal conversion factors:\n");
      std::print(stream, "===============================================================================\n\n");
      std::print(stream, "Energy (internal units) to Kelvin:      {} [-]\n", Units::EnergyToKelvin);
      std::print(stream, "Kelvin to energy (internal units):      {} [-]\n", Units::KelvinToEnergy);
      std::print(stream, "\n\n");
      break;
    case System::ReducedUnits:
      std::print(stream, "Mutual consistent basic set of units:\n");
      std::print(stream, "===============================================================================\n\n");
      std::print(stream, "Unit of temperature: {} [{}]\n", 1, unitOfTemperatureString);
      std::print(stream, "Unit of energy:      {} [{}]\n", Units::EnergyConversionFactor, unitOfEnergyString);
      std::print(stream, "Unit of length:      {} [{}]\n", Units::LengthUnit, unitOfLengthString);
      std::print(stream, "Unit of mass:        {} [{}]\n", Units::MassUnit, unitOfMassString);
      std::print(stream, "Unit of charge:      {} [{}]\n", Units::ChargeUnit, unitOfChargeString);
      std::print(stream, "\n\n");

      std::print(stream, "Boltzmann constant (internal units):          {} [-]\n", Units::KB);
      std::print(stream, "represents conversion from energy (internal units) to temperature (in Kelvin)\n\n");

      std::print(stream, "Derived units and their conversion factors:\n");
      std::print(stream, "===============================================================================\n\n");
      std::print(stream, "Unit of time:                             {} [{}]\n", Units::TimeUnit, unitOfTimeString);
      std::print(stream, "Unit of velocity:                         {} [{}]\n", Units::VelocityConversionFactor,
                 unitOfVelocityString);
      std::print(stream, "Unit of force:                            {} [{}]\n", Units::ForceConversionFactor,
                 unitOfForceString);
      std::print(stream, "Unit of diffusion:                        {} [{}]\n", Units::DiffusionConversionFactor,
                 unitOfDiffusionString);
      std::print(stream, "Unit of acceleration:                     {} [{}]\n", Units::AccelerationConversionFactor,
                 unitOfAccelerationString);
      std::print(stream, "Unit of torque:                           {} [{}]\n", Units::TorqueConversionFactor,
                 unitOfTorqueString);
      std::print(stream, "Unit of pressure:                         {} [{}]\n", Units::PressureConversionFactor,
                 unitOfPressureString);
      std::print(stream, "Unit of volume:                           {} [{}]\n", Units::VolumeConversionFactor,
                 unitOfVolumeString);
      std::print(stream, "Unit of density:                          {} [{}]\n", Units::DensityConversionFactor,
                 unitOfDensityString);
      std::print(stream, "Unit of dynamic viscosity:                {} [{}]\n", Units::DynamicViscosityConversionFactor,
                 unitOfDynamicViscosityString);
      std::print(stream, "Unit of enthalpy:                         {} [{}]\n", Units::EnthalpyConversionFactor,
                 unitOfEnthalpyString);
      std::print(stream, "Unit of polarizability:                   {} [{}]\n", Units::PolarizilibityConversionFactor,
                 unitOfPolarizabilityString);
      std::print(stream, "Unit of Coulomb potential:                {} [{}]\n",
                 Units::CoulombicConversionFactor * Units::EnergyToKelvin, unitOfCoulombPotentialString);
      std::print(stream, "Unit of dielectric constant:              {} [{}]\n",
                 Units::DielectricConstantConversionFactor, unitOfDielectricConstantString);
      std::print(stream, "Unit of dipole moment:                    {} [{}]\n", Units::DipoleMomentConversionFactor,
                 unitOfDipoleMomentString);
      std::print(stream, "Unit of Debye:                            {} [{}]\n", Units::DebyeConversionFactor,
                 unitOfDebyeString);
      std::print(stream, "Unit of electric potential:               {} [{}]\n",
                 Units::ElectricPotentialConversionFactor, unitOfElectricPotentialString);
      std::print(stream, "Unit of electric field:                   {} [{}]\n", Units::ElectricFieldConversionFactor,
                 unitOfElectricFieldString);
      std::print(stream, "Unit of isothermal compressibility:       {} [{}]\n",
                 Units::IsothermalCompressibilityConversionFactor, unitOfIsothermalCompressibilityString);
      std::print(stream, "Unit of heat capacity:                    {} [{}]\n", Units::HeatCapacityConversionFactor,
                 unitOfHeatCapacityString);
      std::print(stream, "Unit of volumetric expansion coefficient: {} [{}]\n",
                 Units::VolumetricExpansionCoefficientConversionFactor, unitOfVolumetricExpansionCoefficientString);
      std::print(stream, "\n\n");
      break;
  }

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
