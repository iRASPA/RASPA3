module;

#ifdef USE_LEGACY_HEADERS
#include <numbers>
#include <print>
#include <sstream>
#include <string>
#endif

export module units;

#ifndef USE_LEGACY_HEADERS
import <numbers>;
import <string>;
import <sstream>;
import <print>;
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
export namespace Units
{
// Fundamental constants of nature
/// The length of one angstrom in meters.
constexpr double Angstrom = 1e-10;
/// The volume of one cubic angstrom in cubic meters.
constexpr double AngstromCubed = 1e-30;
/// The length of one nanosecond in seconds.
constexpr double NanoSecond = 1e-9;
/// The length of one picosecond in seconds.
constexpr double PicoSecond = 1e-12;
/// The length of one femtosecond in seconds.
constexpr double FemtoSecond = 1e-15;
/// The mass of one atomic mass unit in kilograms.
constexpr double AtomicMassUnit = 1.6605402e-27;  // kg
/// The charge of an electron in coulombs per particle.
constexpr double ElectronicChargeUnit = 1.60217733e-19;  // C/particle
/// The molar gas constant in joules per mole per kelvin.
constexpr double MolarGasConstant = 8.314464919;  // J mol^-1 K^-1
/// The Boltzmann constant in joules per kelvin.
constexpr double BoltzmannConstant = 1.380650324e-23;  // J K^-1
/// The speed of light in meters per second.
constexpr double SpeedOfLight = 299792458.0;  // m s^-1
/// Avogadro's constant in particles per mole.
constexpr double AvogadroConstant = 6.0221419947e23;  // mol^-1
/// The electric constant (vacuum permittivity) in coulombs squared per newton square meters.
constexpr double ElectricConstant = 8.8541878176e-12;  // C^2/(N.m^2)
/// Planck's constant in joule seconds.
constexpr double PlanckConstant = 6.6260687652e-34;  // J.s
/// The vacuum dielectric constant in farads per meter.
constexpr double DielectricConstantVacuum = 8.8541878176e-12;  // F/m
/// The value of one Debye in coulomb meters.
constexpr double Debye = 3.335640952e-30;  // C.m

// Set of base units
/// The simulation's base unit of length in meters.
constexpr double LengthUnit = Angstrom;
/// The simulation's base unit of time in seconds.
constexpr double TimeUnit = PicoSecond;
/// The simulation's base unit of mass in kilograms.
constexpr double MassUnit = AtomicMassUnit;
/// The simulation's base unit of charge in coulombs per particle.
constexpr double ChargeUnit = ElectronicChargeUnit;

// Derived units and their conversion factors
/// Conversion factor for mass from simulation units to kilograms.
constexpr double MassConversionFactor = AtomicMassUnit;  // kg
/// Conversion factor for volume from simulation units to cubic meters.
constexpr double VolumeConversionFactor = LengthUnit * LengthUnit * LengthUnit;  // m^-3
/// Conversion factor for density from simulation units to kilograms per cubic meter.
constexpr double DensityConversionfactor = MassConversionFactor / VolumeConversionFactor;  // kg/m^3
/// Conversion factor for energy from simulation units to joules.
constexpr double EnergyConversionFactor = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit);  // J
/// Conversion factor for force from simulation units to newtons.
constexpr double ForceConversionFactor = EnergyConversionFactor / LengthUnit;  // N (=m/s^2)
/// Conversion factor for pressure from simulation units to pascals.
constexpr double PressureConversionFactor = MassUnit / (LengthUnit * TimeUnit * TimeUnit);  // Pa(=N/m^2)
/// Conversion factor for position from simulation units to meters.
constexpr double PositionConversionFactor = LengthUnit;  // m
/// Conversion factor for time from simulation units to seconds.
constexpr double TimeConversionFactor = TimeUnit;  // s
/// Conversion factor for diffusion from simulation units to square meters per second.
constexpr double DiffusionConversionFactor = LengthUnit * LengthUnit / TimeUnit;  // m^2/s
/// Conversion factor for velocity from simulation units to meters per second.
constexpr double VelocityConversionFactor = LengthUnit / TimeUnit;  // m/s
/// Conversion factor for acceleration from simulation units to meters per square second.
constexpr double AccelerationConversionFactor = LengthUnit * LengthUnit / TimeUnit;  // m^2/s
/// Conversion factor for dipole moment from simulation units to coulomb meters.
constexpr double DipoleMomentConversionFactor = ChargeUnit * LengthUnit;  // C.m
/// Conversion factor from simulation units to Debye units.
constexpr double DebyeConversionFactor = ChargeUnit * LengthUnit / Debye;  // C.m
/// Conversion factor for electric potential from simulation units to volts.
constexpr double ElectricPotentialConversionFactor = (EnergyConversionFactor / ChargeUnit);  // V
/// Conversion factor for polarizability from simulation units to SI units.
constexpr double PolarizilibityConversionFactor =
    ChargeUnit * ChargeUnit * LengthUnit * LengthUnit / EnergyConversionFactor;
/// Conversion factor for electric field from simulation units to volts per meter.
constexpr double ElectricFieldConversionFactor = ElectricPotentialConversionFactor / LengthUnit;  // V
/// Boltzmann constant in simulation units.
constexpr double KB = BoltzmannConstant / EnergyConversionFactor;
/// Conversion factor for Coulomb potential from simulation units to SI units.
constexpr double CoulombicConversionFactor =
    ChargeUnit * ChargeUnit / (4.0 * std::numbers::pi * DielectricConstantVacuum * LengthUnit * EnergyConversionFactor);
/// Conversion factor for dielectric constant from simulation units to SI units.
constexpr double DielectricConstantConversionFactor =
    TimeUnit * TimeUnit * ChargeUnit * ChargeUnit / (AtomicMassUnit * LengthUnit * LengthUnit * LengthUnit);
/// Conversion factor for isothermal compressibility from simulation units to SI units.
constexpr double IsothermalCompressibilityConversionFactor = 1.0 / PressureConversionFactor;
/// Conversion factor for heat capacity from simulation units to joules per mole per kelvin.
constexpr double HeatCapacityConversionFactor = (EnergyConversionFactor * AvogadroConstant);  // J/mol/K
/// Conversion factor for volumetric expansion coefficient from simulation units to SI units.
constexpr double VolumetricExpansionCoefficientConversionFactor = 1.0;  // 1/K
/// Conversion factor for Feynman-Hibbs correction from simulation units to SI units.
constexpr double FeymannHibbsConversionFactor = (PlanckConstant * PlanckConstant / (2.0 * std::numbers::pi)) /
                                                (24.0 * AtomicMassUnit * BoltzmannConstant * LengthUnit * LengthUnit);

// 1 calorie (International Table) = 4.1868 J
// 1 calorie (thermochemical) = 4.184 J
// 1 calorie (15C) = 4.1855 J
/// Conversion factor from calories to joules.
constexpr double CalToJoule = 4.184;
/// Conversion factor from joules to calories.
constexpr double JouleToCal = 1.0 / CalToJoule;
/// Conversion factor from electron volts to kilojoules per mole.
constexpr double EvToKJPerMol = 96.48534;
/// Conversion factor from electron volts to kelvin.
constexpr double EvToKelvin = EvToKJPerMol * 1000.0 / MolarGasConstant;

/// Conversion factor from simulation energy units to kelvin.
constexpr double EnergyToKelvin = EnergyConversionFactor * AvogadroConstant / MolarGasConstant;
/// Conversion factor from kelvin to simulation energy units.
constexpr double KelvinToEnergy = MolarGasConstant / (EnergyConversionFactor * AvogadroConstant);
/// Conversion factor from simulation energy units to kilojoules per mole.
constexpr double EnergyToKJPerMol = EnergyConversionFactor * AvogadroConstant / 1000.0;
/// Conversion factor from simulation energy units to electron volts.
constexpr double EnergyToEV = EnergyToKelvin / 11604.23;
/// Conversion factor from simulation energy units to kilocalories per mole.
constexpr double EnergyToKCalPerMol = EnergyConversionFactor * AvogadroConstant / 4184.0;
/// Conversion factor from kilocalories per mole to simulation energy units.
constexpr double KCalPerMolToEnergy = 4184.0 / (EnergyConversionFactor * AvogadroConstant);

// Pressure conversion factors
/// Conversion factor from pascals to atmospheres.
constexpr double PaToAtm = (1.0 / 101325.0);
/// Conversion factor from pascals to bars.
constexpr double PatoBar = 0.00001;
/// Conversion factor from pascals to torr.
constexpr double PaToTorr = 0.0075;
/// Conversion factor from kilopascals to atmospheres.
constexpr double KPaToAtm = 1.0 / 101.325;
/// Conversion factor from kilopascals to bars.
constexpr double KPaToBar = 0.01;
/// Conversion factor from kilopascals to torr.
constexpr double KPaToTorr = 7.5;
/// Conversion factor from bars to torr.
constexpr double BarToTorr = 750.062;
/// Conversion factor from bars to atmospheres.
constexpr double BarToAtm = 0.9869;
/// Conversion factor from bars to kilopascals.
constexpr double BarTokPa = 100.0;
/// Conversion factor from atmospheres to torr.
constexpr double AtmToTorr = 760;
/// Conversion factor from atmospheres to bars.
constexpr double AtmToBar = 1.0 / 0.9869;
/// Conversion factor from atmospheres to kilopascals.
constexpr double AtmTokPa = 101.325;
/// Conversion factor from atmospheres to pascals.
constexpr double AtmToPa = 101325.0;

/**
 * \brief Generates a string representation of the current unit settings.
 *
 * Returns a formatted string that lists the base units, derived units, and their conversion factors
 * used in the simulation. This can be used to output the unit settings to logs or console.
 *
 * \return A string containing the unit settings.
 */
std::string printStatus();

/**
 * \brief Generates a JSON object containing the current unit settings.
 *
 * Returns a JSON object that includes the base units, derived units, and their conversion factors.
 * This can be used to serialize the unit settings for output to files or for inter-process communication.
 *
 * \return A JSON object with the unit settings.
 */
nlohmann::json jsonStatus();
};  // namespace Units
