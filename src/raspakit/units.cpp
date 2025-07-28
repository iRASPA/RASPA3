module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <numbers>
#include <ostream>
#include <print>
#include <sstream>
#include <string>
#endif

module units;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import json;

Units::System Units::unitSystem = Units::System::RASPA;

// Set of base units
/// The simulation's base unit of length in meters.
double Units::LengthUnit = Angstrom;
/// The simulation's base unit of time in seconds.
double Units::TimeUnit = PicoSecond;
/// The simulation's base unit of mass in kilograms.
double Units::MassUnit = AtomicMassUnit;
/// The simulation's base unit of charge in coulombs per particle.
double Units::ChargeUnit = ElectronicChargeUnit;

/// Conversion factor for energy from simulation units to J.
double Units::EnergyUnit =
    Units::MassUnit * Units::LengthUnit * Units::LengthUnit / (Units::TimeUnit * Units::TimeUnit);

/// Derived units and their conversion factors
///
/// Conversion factor for length from simulation units to meter.
double Units::LengthConversionFactor = Units::LengthUnit;  // m
/// Conversion factor for energy from simulation units to joules.
double Units::EnergyConversionFactor =
    Units::MassUnit * Units::LengthUnit * Units::LengthUnit / (Units::TimeUnit * Units::TimeUnit);  // J
/// Conversion factor for mass from simulation units to kilograms.
double Units::MassConversionFactor = Units::AtomicMassUnit;  // kg
/// Conversion factor for charge from simulation units to Coulomb per particle.
double Units::ChargeConversionFactor = Units::ChargeUnit;  // C/particle
/// Conversion factor for time from simulation units to seconds.
double Units::TimeConversionFactor = Units::TimeUnit;  // s
/// Conversion factor for velocity from simulation units to meters per second.
double Units::VelocityConversionFactor = Units::LengthUnit / Units::TimeUnit;  // m/s
/// Conversion factor for force from simulation units to newtons.
double Units::ForceConversionFactor = Units::EnergyConversionFactor / Units::LengthUnit;  // N (=m/s^2)
/// Conversion factor for diffusion from simulation units to square meters per second.
double Units::DiffusionConversionFactor = Units::LengthUnit * Units::LengthUnit / Units::TimeUnit;  // m^2/s
/// Conversion factor for acceleration from simulation units to meters per square second.
double Units::AccelerationConversionFactor = Units::LengthUnit * Units::LengthUnit / Units::TimeUnit;  // m^2/s
/// Conversion factor for torque from simulation units to Newton meters.
double Units::TorqueConversionFactor =
    Units::MassUnit * Units::LengthUnit * Units::LengthUnit / (Units::TimeUnit * Units::TimeUnit);  // kg.m^2.s^-2
/// Conversion factor for temperature from simulation units to K
double Units::TemperatureConversionFactor = 1.0;

double Units::KB = Units::BoltzmannConstant / Units::EnergyConversionFactor;
/// Conversion factor for pressure from simulation units to pascals.
double Units::PressureConversionFactor =
    Units::MassUnit / (Units::LengthUnit * Units::TimeUnit * Units::TimeUnit);  // Pa(=N/m^2)
/// Conversion factor for volume from simulation units to cubic meters.
double Units::VolumeConversionFactor = Units::LengthUnit * Units::LengthUnit * Units::LengthUnit;  // m^-3
/// Conversion factor for density from simulation units to kilograms per cubic meter.
double Units::DensityConversionFactor = Units::MassConversionFactor / Units::VolumeConversionFactor;  // kg/m^3
/// Conversion factor for dynamic viscosity from simulation units to Pascal second.
double Units::DynamicViscosityConversionFactor =
    Units::MassConversionFactor / Units::LengthUnit / Units::TimeUnit;  // kg.m^-1.s^-1
/// Conversion factor for enthalpy from simulation units to J
double Units::EnthalpyConversionFactor =
    Units::MassUnit * Units::LengthUnit * Units::LengthUnit / (Units::TimeUnit * Units::TimeUnit);  // J
/// Conversion factor for polarizability from simulation units to SI units.
double Units::PolarizilibityConversionFactor =
    Units::ChargeUnit * Units::ChargeUnit * Units::LengthUnit * Units::LengthUnit / Units::EnergyConversionFactor;
/// Conversion factor for dielectric constant from simulation units to SI units.
double Units::DielectricConstantConversionFactor =
    Units::TimeUnit * Units::TimeUnit * Units::ChargeUnit * Units::ChargeUnit /
    (Units::AtomicMassUnit * Units::LengthUnit * Units::LengthUnit * Units::LengthUnit);
/// Conversion factor for Coulomb potential from simulation units to SI units.
double Units::CoulombicConversionFactor =
    Units::ChargeUnit * Units::ChargeUnit /
    (4.0 * std::numbers::pi * Units::DielectricConstantVacuum * Units::LengthUnit * Units::EnergyConversionFactor);
/// Conversion factor for position from simulation units to meters.
double Units::DipoleMomentConversionFactor = Units::ChargeUnit * Units::LengthUnit;  // C.m
/// Conversion factor from simulation units to Debye units.
double Units::DebyeConversionFactor = Units::ChargeUnit * Units::LengthUnit / Units::Debye;  // C.m
/// Conversion factor for electric potential from simulation units to volts.
double Units::ElectricPotentialConversionFactor = (Units::EnergyConversionFactor / Units::ChargeUnit);  // V
/// Conversion factor for electric field from simulation units to volts per meter.
double Units::ElectricFieldConversionFactor = Units::ElectricPotentialConversionFactor / Units::LengthUnit;  // V
/// Conversion factor for isothermal compressibility from simulation units to per Pascal
double Units::IsothermalCompressibilityConversionFactor = 1.0 / Units::PressureConversionFactor;  // Pa^-1
/// Conversion factor for heat capacity from simulation units to joules per mole per kelvin.
double Units::HeatCapacityConversionFactor = (Units::EnergyConversionFactor * Units::AvogadroConstant);  // J/mol/K
/// Conversion factor for volumetric expansion coefficient from simulation units to SI units.
double Units::VolumetricExpansionCoefficientConversionFactor = 1.0;  // 1/K
/// Conversion factor for Feynman-Hibbs correction from simulation units to SI units.
double Units::FeymannHibbsConversionFactor =
    (Units::PlanckConstant * Units::PlanckConstant / (2.0 * std::numbers::pi)) /
    (24.0 * Units::AtomicMassUnit * Units::BoltzmannConstant * Units::LengthUnit * Units::LengthUnit);

/// Conversion factor from simulation energy units to kelvin.
double Units::EnergyToKelvin = Units::EnergyConversionFactor * Units::AvogadroConstant / Units::MolarGasConstant;
/// Conversion factor from kelvin to simulation energy units.
double Units::KelvinToEnergy = Units::MolarGasConstant / (Units::EnergyConversionFactor * Units::AvogadroConstant);
/// Conversion factor from simulation energy units to kilojoules per mole.
double Units::EnergyToKJPerMol = Units::EnergyConversionFactor * Units::AvogadroConstant / 1000.0;
/// Conversion factor from simulation energy units to electron volts.
double Units::EnergyToEV = Units::EnergyToKelvin / 11604.23;
/// Conversion factor from simulation energy units to kilocalories per mole.
double Units::EnergyToKCalPerMol = Units::EnergyConversionFactor * Units::AvogadroConstant / 4184.0;
/// Conversion factor from kilocalories per mole to simulation energy units.
double Units::KCalPerMolToEnergy = 4184.0 / (Units::EnergyConversionFactor * Units::AvogadroConstant);

std::string Units::unitOfLengthString{"m"};
std::string Units::unitOfEnergyString{"J"};
std::string Units::unitOfEnergyPerMolString{"J.mol⁻¹"};
std::string Units::unitOfMassString{"kg"};
std::string Units::unitOfChargeString{"C/particle"};
std::string Units::unitOfTimeString{"s"};

std::string Units::displayedUnitOfLengthString{"Å"};
std::string Units::displayedUnitOfEnergyString{"K"};
std::string Units::displayedUnitOfEnergyConversionString = "/kʙ";

std::string Units::unitOfVelocityString{"m.s⁻¹"};
std::string Units::unitOfForceString{"N"};
std::string Units::unitOfDiffusionString{"m².s⁻¹"};
std::string Units::unitOfAccelerationString{"m².s⁻¹"};
std::string Units::unitOfTorqueString{"N.m"};
std::string Units::unitOfTemperatureString{"K"};
std::string Units::unitOfBoltzmannFactorString{"m².kg.s⁻².K⁻¹"};
std::string Units::unitOfPressureString{"Pa"};
std::string Units::unitOfVolumeString{"m³"};
std::string Units::unitOfDensityString{"kg.m⁻³"};
std::string Units::unitOfDynamicViscosityString{"N.s.m⁻²"};
std::string Units::unitOfEnthalpyString{"J"};
std::string Units::unitOfEnthalpyPerMolString{"J.mol⁻¹"};
std::string Units::unitOfPolarizabilityString{"-"};
std::string Units::unitOfCoulombPotentialString{"K"};
std::string Units::unitOfDielectricConstantString{"s².C².kg⁻¹.m⁻³"};
std::string Units::unitOfDipoleMomentString{"C.m"};
std::string Units::unitOfDebyeString{"C.m"};
std::string Units::unitOfElectricPotentialString{"V"};
std::string Units::unitOfElectricFieldString{"V.m⁻¹"};
std::string Units::unitOfIsothermalCompressibilityString{"Pa⁻¹"};
std::string Units::unitOfHeatCapacityString{"J.mol⁻¹.K⁻¹"};
std::string Units::unitOfVolumetricExpansionCoefficientString{"K⁻¹"};

void Units::setUnits(Units::System unit_system)
{
  unitSystem = unit_system;
  switch (unit_system)
  {
    case System::RASPA:
      LengthUnit = Angstrom;
      TimeUnit = PicoSecond;
      MassUnit = AtomicMassUnit;
      ChargeUnit = ElectronicChargeUnit;

      EnergyUnit = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit);

      LengthConversionFactor = LengthUnit;                                                  // m
      EnergyConversionFactor = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit);  // J
      MassConversionFactor = AtomicMassUnit;                                                // kg
      ChargeConversionFactor = ChargeUnit;                                                  // C/particle
      TimeConversionFactor = TimeUnit;                                                      // s
      VelocityConversionFactor = LengthUnit / TimeUnit;                                     // m/s
      ForceConversionFactor = EnergyConversionFactor / LengthUnit;                          // N (=m/s^2)
      DiffusionConversionFactor = LengthUnit * LengthUnit / TimeUnit;                       // m^2/s
      AccelerationConversionFactor = LengthUnit * LengthUnit / TimeUnit;                    // m^2/s
      TorqueConversionFactor = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit);  // kg.m^2.s^-2
      TemperatureConversionFactor = 1.0;

      KB = BoltzmannConstant / EnergyConversionFactor;
      PressureConversionFactor = MassUnit / (LengthUnit * TimeUnit * TimeUnit);               // Pa(=N/m^2)
      VolumeConversionFactor = LengthUnit * LengthUnit * LengthUnit;                          // m^-3
      DensityConversionFactor = MassConversionFactor / VolumeConversionFactor;                // kg/m^3
      DynamicViscosityConversionFactor = MassConversionFactor / LengthUnit / TimeUnit;        // kg.m^-1.s^-1
      EnthalpyConversionFactor = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit);  // J
      PolarizilibityConversionFactor = ChargeUnit * ChargeUnit * LengthUnit * LengthUnit / EnergyConversionFactor;
      DielectricConstantConversionFactor =
          TimeUnit * TimeUnit * ChargeUnit * ChargeUnit / (AtomicMassUnit * LengthUnit * LengthUnit * LengthUnit);
      CoulombicConversionFactor =
          ChargeUnit * ChargeUnit /
          (4.0 * std::numbers::pi * DielectricConstantVacuum * LengthUnit * EnergyConversionFactor);
      DipoleMomentConversionFactor = ChargeUnit * LengthUnit;                          // C.m
      DebyeConversionFactor = ChargeUnit * LengthUnit / Debye;                         // C.m
      ElectricPotentialConversionFactor = (EnergyConversionFactor / ChargeUnit);       // V
      ElectricFieldConversionFactor = ElectricPotentialConversionFactor / LengthUnit;  // V
      IsothermalCompressibilityConversionFactor = 1.0 / PressureConversionFactor;      // Pa^-1
      HeatCapacityConversionFactor = (EnergyConversionFactor * AvogadroConstant);      // J/mol/K
      VolumetricExpansionCoefficientConversionFactor = 1.0;                            // 1/K
      FeymannHibbsConversionFactor = (PlanckConstant * PlanckConstant / (2.0 * std::numbers::pi)) /
                                     (24.0 * AtomicMassUnit * BoltzmannConstant * LengthUnit * LengthUnit);

      EnergyToKelvin = EnergyConversionFactor * AvogadroConstant / MolarGasConstant;
      KelvinToEnergy = MolarGasConstant / (EnergyConversionFactor * AvogadroConstant);
      EnergyToKJPerMol = EnergyConversionFactor * AvogadroConstant / 1000.0;
      EnergyToEV = EnergyToKelvin / 11604.23;
      EnergyToKCalPerMol = EnergyConversionFactor * AvogadroConstant / 4184.0;
      KCalPerMolToEnergy = 4184.0 / (EnergyConversionFactor * AvogadroConstant);

      unitOfLengthString = "m";
      unitOfEnergyString = "J";
      unitOfEnergyPerMolString = "J.mol⁻¹";
      unitOfMassString = "kg";
      unitOfChargeString = "C/particle";
      unitOfTimeString = "s";

      displayedUnitOfLengthString = "Å";
      displayedUnitOfEnergyString = "K";
      displayedUnitOfEnergyConversionString = "/kʙ";

      unitOfVelocityString = "m.s⁻¹";
      unitOfForceString = "N";
      unitOfDiffusionString = "m².s⁻¹";
      unitOfAccelerationString = "m².s⁻¹";
      unitOfTorqueString = "N.m";
      unitOfTemperatureString = "K";
      unitOfBoltzmannFactorString = "m².kg.s⁻².K⁻¹";
      unitOfPressureString = "Pa";
      unitOfVolumeString = "m³";
      unitOfDensityString = "kg.m⁻³";
      unitOfDynamicViscosityString = "N.s.m⁻²";
      unitOfEnthalpyString = "J";
      unitOfEnthalpyPerMolString = "J.mol⁻¹";
      unitOfPolarizabilityString = "-";
      unitOfCoulombPotentialString = "K";
      unitOfDielectricConstantString = "s².C².kg⁻¹.m⁻³";
      unitOfDipoleMomentString = "C.m";
      unitOfDebyeString = "C.m";
      unitOfElectricPotentialString = "V";
      unitOfElectricFieldString = "V.m⁻¹";
      unitOfIsothermalCompressibilityString = "Pa⁻¹";
      unitOfHeatCapacityString = "J.mol⁻¹.K⁻¹";
      unitOfVolumetricExpansionCoefficientString = "K⁻¹";
      break;
    case System::ReducedUnits:
      LengthUnit = 1.0;
      TimeUnit = 1.0;
      MassUnit = 1.0;
      ChargeUnit = 1.0;
      EnergyUnit = 1.0;

      LengthConversionFactor = 1.0;
      EnergyConversionFactor = 1.0;
      MassConversionFactor = 1.0;
      ChargeConversionFactor = 1.0;
      TimeConversionFactor = 1.0;
      VelocityConversionFactor = 1.0;
      ForceConversionFactor = 1.0;
      DiffusionConversionFactor = 1.0;
      AccelerationConversionFactor = 1.0;
      TorqueConversionFactor = 1.0;
      TemperatureConversionFactor = 1.0;

      KB = 1.0;
      PressureConversionFactor = 1.0;
      VolumeConversionFactor = 1.0;
      DensityConversionFactor = 1.0;
      DynamicViscosityConversionFactor = 1.0;
      EnthalpyConversionFactor = 1.0;
      PolarizilibityConversionFactor = 1.0;
      DielectricConstantConversionFactor = 1.0;
      CoulombicConversionFactor = 1.0;
      DipoleMomentConversionFactor = 1.0;
      DebyeConversionFactor = 1.0;
      ElectricPotentialConversionFactor = 1.0;
      ElectricFieldConversionFactor = 1.0;
      IsothermalCompressibilityConversionFactor = 1.0;
      HeatCapacityConversionFactor = 1.0;
      VolumetricExpansionCoefficientConversionFactor = 1.0;
      FeymannHibbsConversionFactor = 1.0;

      EnergyToKelvin = 1.0;
      KelvinToEnergy = 1.0;
      EnergyToKJPerMol = 0.0;
      EnergyToEV = 0.0;
      EnergyToKCalPerMol = 0.0;
      KCalPerMolToEnergy = 0.0;

      unitOfLengthString = "σ";
      unitOfEnergyString = "ε";
      unitOfMassString = "m";
      unitOfChargeString = "q′";
      unitOfTimeString = "𝜏=σ√(m/ε)";

      displayedUnitOfLengthString = "σ";
      displayedUnitOfEnergyString = "ε";
      displayedUnitOfEnergyConversionString = "   ";

      unitOfVelocityString = "σ.𝜏⁻¹";
      unitOfForceString = "ε.σ⁻¹";
      unitOfDiffusionString = "ε.σ⁻¹.m⁻¹";
      unitOfAccelerationString = "ε.σ⁻¹.m⁻¹";
      unitOfTorqueString = "ε";
      unitOfTemperatureString = "ε.kʙ⁻¹";
      unitOfBoltzmannFactorString = "kʙ";
      unitOfPressureString = "ε.σ⁻³";
      unitOfVolumeString = "σ³";
      unitOfDensityString = "m.σ⁻³";
      unitOfDynamicViscosityString = "ε.𝜏.σ⁻³";
      unitOfEnthalpyString = "ε";
      unitOfEnthalpyPerMolString = "";
      unitOfPolarizabilityString = "-";
      unitOfCoulombPotentialString = "ε";
      unitOfDielectricConstantString = "q².ε⁻¹.σ⁻¹";
      unitOfDipoleMomentString = "√(4πε₀σ³ε)";
      unitOfDebyeString = "q.σ";
      unitOfElectricPotentialString = "ε.q⁻¹";
      unitOfElectricFieldString = "ε.σ⁻¹.(√(4πε₀σε))⁻¹";
      unitOfIsothermalCompressibilityString = "σ³.ε⁻¹";
      unitOfHeatCapacityString = "kʙ";
      unitOfVolumetricExpansionCoefficientString = "kʙ.ε⁻¹";
      break;
  }
}

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
