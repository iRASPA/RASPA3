module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <numbers>
#include <print>
#include <sstream>
#include <string>
#endif

export module units;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import json;
/**
 * \brief Contains fundamental constants and unit conversion factors used throughout the simulation.
 *
 * The Units namespace provides definitions for various physical constants, base units, and derived units.
 * It includes fundamental constants such as the speed of light, Boltzmann constant, and Planck constant.
 * Additionally, it defines the simulation's base units for length, time, mass, and charge, and calculates
 * derived units and their conversion factors, which are used to convert between simulation units and SI units.
 * The namespace also provides functions to output the unit settings in human-readable or JSON format.
 */
export struct Units
{
  enum class System : std::size_t
  {
    RASPA = 0,
    ReducedUnits = 1
  };

  static Units::System unitSystem;

  // Fundamental constants of nature
  /// The length of one angstrom in meters.
  static constexpr double Angstrom = 1e-10;
  /// The volume of one cubic angstrom in cubic meters.
  static constexpr double AngstromCubed = 1e-30;
  /// The length of one nanosecond in seconds.
  static constexpr double NanoSecond = 1e-9;
  /// The length of one picosecond in seconds.
  static constexpr double PicoSecond = 1e-12;
  /// The length of one femtosecond in seconds.
  static constexpr double FemtoSecond = 1e-15;
  /// The mass of one atomic mass unit in kilograms.
  static constexpr double AtomicMassUnit = 1.6605402e-27;  // kg
  /// The charge of an electron in coulombs per particle.
  static constexpr double ElectronicChargeUnit = 1.60217733e-19;  // C/particle
  /// The molar gas constant in joules per mole per kelvin.
  static constexpr double MolarGasConstant = 8.314464919;  // J mol^-1 K^-1
  /// The Boltzmann constant in joules per kelvin.
  static constexpr double BoltzmannConstant = 1.380650324e-23;  // J K^-1
  /// The speed of light in meters per second.
  static constexpr double SpeedOfLight = 299792458.0;  // m s^-1
  /// Avogadro's constant in particles per mole.
  static constexpr double AvogadroConstant = 6.0221419947e23;  // mol^-1
  /// The electric constant (vacuum permittivity) in coulombs squared per newton square meters.
  static constexpr double ElectricConstant = 8.8541878176e-12;  // C^2/(N.m^2)
  /// Planck's constant in joule seconds.
  static constexpr double PlanckConstant = 6.6260687652e-34;  // J.s
  /// The vacuum dielectric constant in farads per meter.
  static constexpr double DielectricConstantVacuum = 8.8541878176e-12;  // F/m
  /// The value of one Debye in coulomb meters.
  static constexpr double Debye = 3.335640952e-30;  // C.m

  // 1 calorie (International Table) = 4.1868 J
  // 1 calorie (thermochemical) = 4.184 J
  // 1 calorie (15C) = 4.1855 J
  /// Conversion factor from calories to joules.
  static constexpr double CalToJoule = 4.184;
  /// Conversion factor from joules to calories.
  static constexpr double JouleToCal = 1.0 / CalToJoule;
  /// Conversion factor from electron volts to kilojoules per mole.
  static constexpr double EvToKJPerMol = 96.48534;
  /// Conversion factor from electron volts to kelvin.
  static constexpr double EvToKelvin = EvToKJPerMol * 1000.0 / MolarGasConstant;

  // Pressure conversion factors
  /// Conversion factor from pascals to atmospheres.
  static constexpr double PaToAtm = (1.0 / 101325.0);
  /// Conversion factor from pascals to bars.
  static constexpr double PatoBar = 0.00001;
  /// Conversion factor from pascals to torr.
  static constexpr double PaToTorr = 0.0075;
  /// Conversion factor from kilopascals to atmospheres.
  static constexpr double KPaToAtm = 1.0 / 101.325;
  /// Conversion factor from kilopascals to bars.
  static constexpr double KPaToBar = 0.01;
  /// Conversion factor from kilopascals to torr.
  static constexpr double KPaToTorr = 7.5;
  /// Conversion factor from bars to torr.
  static constexpr double BarToTorr = 750.062;
  /// Conversion factor from bars to atmospheres.
  static constexpr double BarToAtm = 0.9869;
  /// Conversion factor from bars to kilopascals.
  static constexpr double BarTokPa = 100.0;
  /// Conversion factor from atmospheres to torr.
  static constexpr double AtmToTorr = 760;
  /// Conversion factor from atmospheres to bars.
  static constexpr double AtmToBar = 1.0 / 0.9869;
  /// Conversion factor from atmospheres to kilopascals.
  static constexpr double AtmTokPa = 101.325;
  /// Conversion factor from atmospheres to pascals.
  static constexpr double AtmToPa = 101325.0;

  static constexpr double DegreesToRadians = std::numbers::pi / 180.0;
  static constexpr double RadiansToDegrees = 180.0 / std::numbers::pi;

  // Set of base units
  /// The simulation's base unit of length in meters.
  static double LengthUnit;
  /// The simulation's base unit of time in seconds.
  static double TimeUnit;
  /// The simulation's base unit of mass in kilograms.
  static double MassUnit;
  /// The simulation's base unit of charge in coulombs per particle.
  static double ChargeUnit;

  /// Conversion factor for energy from simulation units to J.
  static double EnergyUnit;

  /// Derived units and their conversion factors
  ///
  /// Conversion factor for length from simulation units to meter.
  static double LengthConversionFactor;
  /// Conversion factor for energy from simulation units to joules.
  static double EnergyConversionFactor;
  /// Conversion factor for mass from simulation units to kilograms.
  static double MassConversionFactor;
  /// Conversion factor for charge from simulation units to Coulomb per particle.
  static double ChargeConversionFactor;
  /// Conversion factor for time from simulation units to seconds.
  static double TimeConversionFactor;
  /// Conversion factor for velocity from simulation units to meters per second.
  static double VelocityConversionFactor;
  /// Conversion factor for force from simulation units to newtons.
  static double ForceConversionFactor;
  /// Conversion factor for diffusion from simulation units to square meters per second.
  static double DiffusionConversionFactor;
  /// Conversion factor for acceleration from simulation units to meters per square second.
  static double AccelerationConversionFactor;
  /// Conversion factor for torque from simulation units to Newton meters.
  static double TorqueConversionFactor;
  /// Conversion factor for temperature from simulation units to K
  static double TemperatureConversionFactor;

  static double KB;
  /// Conversion factor for pressure from simulation units to pascals.
  static double PressureConversionFactor;
  /// Conversion factor for volume from simulation units to cubic meters.
  static double VolumeConversionFactor;
  /// Conversion factor for density from simulation units to kilograms per cubic meter.
  static double DensityConversionFactor;
  /// Conversion factor for dynamic viscosity from simulation units to Pascal second.
  static double DynamicViscosityConversionFactor;
  /// Conversion factor for enthalpy from simulation units to J
  static double EnthalpyConversionFactor;
  /// Conversion factor for polarizability from simulation units to SI units.
  static double PolarizilibityConversionFactor;
  /// Conversion factor for dielectric constant from simulation units to SI units.
  static double DielectricConstantConversionFactor;
  /// Conversion factor for Coulomb potential from simulation units to SI units.
  static double CoulombicConversionFactor;
  /// Conversion factor for position from simulation units to meters.
  static double DipoleMomentConversionFactor;
  /// Conversion factor from simulation units to Debye units.
  static double DebyeConversionFactor;
  /// Conversion factor for electric potential from simulation units to volts.
  static double ElectricPotentialConversionFactor;
  /// Conversion factor for electric field from simulation units to volts per meter.
  static double ElectricFieldConversionFactor;
  /// Conversion factor for isothermal compressibility from simulation units to per Pascal
  static double IsothermalCompressibilityConversionFactor;
  /// Conversion factor for heat capacity from simulation units to joules per mole per kelvin.
  static double HeatCapacityConversionFactor;
  /// Conversion factor for volumetric expansion coefficient from simulation units to SI units.
  static double VolumetricExpansionCoefficientConversionFactor;
  /// Conversion factor for Feynman-Hibbs correction from simulation units to SI units.
  static double FeymannHibbsConversionFactor;

  /// Conversion factor from simulation energy units to kelvin.
  static double EnergyToKelvin;
  /// Conversion factor from kelvin to simulation energy units.
  static double KelvinToEnergy;
  /// Conversion factor from simulation energy units to kilojoules per mole.
  static double EnergyToKJPerMol;
  /// Conversion factor from simulation energy units to electron volts.
  static double EnergyToEV;
  /// Conversion factor from simulation energy units to kilocalories per mole.
  static double EnergyToKCalPerMol;
  /// Conversion factor from kilocalories per mole to simulation energy units.
  static double KCalPerMolToEnergy;

  static std::string unitOfLengthString;
  static std::string unitOfEnergyString;
  static std::string unitOfEnergyPerMolString;
  static std::string unitOfMassString;
  static std::string unitOfChargeString;
  static std::string unitOfTimeString;

  static std::string displayedUnitOfLengthString;
  static std::string displayedUnitOfEnergyString;
  static std::string displayedUnitOfEnergyConversionString;

  static std::string unitOfVelocityString;
  static std::string unitOfForceString;
  static std::string unitOfDiffusionString;
  static std::string unitOfAccelerationString;
  static std::string unitOfTorqueString;
  static std::string unitOfTemperatureString;
  static std::string unitOfBoltzmannFactorString;
  static std::string unitOfPressureString;
  static std::string unitOfVolumeString;
  static std::string unitOfDensityString;
  static std::string unitOfDynamicViscosityString;
  static std::string unitOfEnthalpyString;
  static std::string unitOfEnthalpyPerMolString;
  static std::string unitOfPolarizabilityString;
  static std::string unitOfCoulombPotentialString;
  static std::string unitOfDielectricConstantString;
  static std::string unitOfDipoleMomentString;
  static std::string unitOfDebyeString;
  static std::string unitOfElectricPotentialString;
  static std::string unitOfElectricFieldString;
  static std::string unitOfIsothermalCompressibilityString;
  static std::string unitOfHeatCapacityString;
  static std::string unitOfVolumetricExpansionCoefficientString;

  static void setUnits(Units::System unit_system);

  /**
   * \brief Generates a string representation of the current unit settings.
   *
   * Returns a formatted string that lists the base units, derived units, and their conversion factors
   * used in the simulation. This can be used to output the unit settings to logs or console.
   *
   * \return A string containing the unit settings.
   */
  static std::string printStatus();

  /**
   * \brief Generates a JSON object containing the current unit settings.
   *
   * Returns a JSON object that includes the base units, derived units, and their conversion factors.
   * This can be used to serialize the unit settings for output to files or for inter-process communication.
   *
   * \return A JSON object with the unit settings.
   */
  static nlohmann::json jsonStatus();
};  // namespace Units
