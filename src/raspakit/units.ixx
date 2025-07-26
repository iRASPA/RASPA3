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
export namespace Units
{
enum class System : size_t
{
  RASPA = 0,
  ReducedUnits = 1
};

Units::System unitSystem{Units::System::RASPA};

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
double BoltzmannConstant = 1.380650324e-23;  // J K^-1
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

// Set of base units
/// The simulation's base unit of length in meters.
inline double LengthUnit = Angstrom;
/// The simulation's base unit of time in seconds.
inline double TimeUnit = PicoSecond;
/// The simulation's base unit of mass in kilograms.
inline double MassUnit = AtomicMassUnit;
/// The simulation's base unit of charge in coulombs per particle.
inline double ChargeUnit = ElectronicChargeUnit;

/// Conversion factor for energy from simulation units to J.
inline double EnergyUnit = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit);

/// Derived units and their conversion factors
///
/// Conversion factor for length from simulation units to meter.
inline double LengthConversionFactor = LengthUnit;  // m
/// Conversion factor for energy from simulation units to joules.
inline double EnergyConversionFactor = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit);  // J
/// Conversion factor for mass from simulation units to kilograms.
inline double MassConversionFactor = AtomicMassUnit;  // kg
/// Conversion factor for charge from simulation units to Coulomb per particle.
inline double ChargeConversionFactor = ChargeUnit;  // C/particle
/// Conversion factor for time from simulation units to seconds.
inline double TimeConversionFactor = TimeUnit;  // s
/// Conversion factor for velocity from simulation units to meters per second.
inline double VelocityConversionFactor = LengthUnit / TimeUnit;  // m/s
/// Conversion factor for force from simulation units to newtons.
inline double ForceConversionFactor = EnergyConversionFactor / LengthUnit;  // N (=m/s^2)
/// Conversion factor for diffusion from simulation units to square meters per second.
inline double DiffusionConversionFactor = LengthUnit * LengthUnit / TimeUnit;  // m^2/s
/// Conversion factor for acceleration from simulation units to meters per square second.
inline double AccelerationConversionFactor = LengthUnit * LengthUnit / TimeUnit;  // m^2/s
/// Conversion factor for torque from simulation units to Newton meters.
inline double TorqueConversionFactor = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit);  // kg.m^2.s^-2
/// Conversion factor for temperature from simulation units to K
inline double TemperatureConversionFactor = 1.0;

inline double KB = BoltzmannConstant / EnergyConversionFactor;
/// Conversion factor for pressure from simulation units to pascals.
inline double PressureConversionFactor = MassUnit / (LengthUnit * TimeUnit * TimeUnit);  // Pa(=N/m^2)
/// Conversion factor for volume from simulation units to cubic meters.
inline double VolumeConversionFactor = LengthUnit * LengthUnit * LengthUnit;  // m^-3
/// Conversion factor for density from simulation units to kilograms per cubic meter.
inline double DensityConversionFactor = MassConversionFactor / VolumeConversionFactor;  // kg/m^3
/// Conversion factor for dynamic viscosity from simulation units to Pascal second.
inline double DynamicViscosityConversionFactor = MassConversionFactor / LengthUnit / TimeUnit;  // kg.m^-1.s^-1
/// Conversion factor for enthalpy from simulation units to J
inline double EnthalpyConversionFactor = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit);  // J
/// Conversion factor for polarizability from simulation units to SI units.
inline double PolarizilibityConversionFactor = ChargeUnit * ChargeUnit * LengthUnit * LengthUnit / EnergyConversionFactor;
/// Conversion factor for dielectric constant from simulation units to SI units.
inline double DielectricConstantConversionFactor =
    TimeUnit * TimeUnit * ChargeUnit * ChargeUnit / (AtomicMassUnit * LengthUnit * LengthUnit * LengthUnit);
/// Conversion factor for Coulomb potential from simulation units to SI units.
inline double CoulombicConversionFactor =
    ChargeUnit * ChargeUnit / (4.0 * std::numbers::pi * DielectricConstantVacuum * LengthUnit * EnergyConversionFactor);
/// Conversion factor for position from simulation units to meters.
inline double DipoleMomentConversionFactor = ChargeUnit * LengthUnit;  // C.m
/// Conversion factor from simulation units to Debye units.
inline double DebyeConversionFactor = ChargeUnit * LengthUnit / Debye;  // C.m
/// Conversion factor for electric potential from simulation units to volts.
inline double ElectricPotentialConversionFactor = (EnergyConversionFactor / ChargeUnit);  // V
/// Conversion factor for electric field from simulation units to volts per meter.
inline double ElectricFieldConversionFactor = ElectricPotentialConversionFactor / LengthUnit;  // V
/// Conversion factor for isothermal compressibility from simulation units to per Pascal
inline double IsothermalCompressibilityConversionFactor = 1.0 / PressureConversionFactor;  // Pa^-1
/// Conversion factor for heat capacity from simulation units to joules per mole per kelvin.
inline double HeatCapacityConversionFactor = (EnergyConversionFactor * AvogadroConstant);  // J/mol/K
/// Conversion factor for volumetric expansion coefficient from simulation units to SI units.
inline double VolumetricExpansionCoefficientConversionFactor = 1.0;  // 1/K
/// Conversion factor for Feynman-Hibbs correction from simulation units to SI units.
inline double FeymannHibbsConversionFactor = (PlanckConstant * PlanckConstant / (2.0 * std::numbers::pi)) /
                                      (24.0 * AtomicMassUnit * BoltzmannConstant * LengthUnit * LengthUnit);

/// Conversion factor from simulation energy units to kelvin.
inline double EnergyToKelvin = EnergyConversionFactor * AvogadroConstant / MolarGasConstant;
/// Conversion factor from kelvin to simulation energy units.
inline double KelvinToEnergy = MolarGasConstant / (EnergyConversionFactor * AvogadroConstant);
/// Conversion factor from simulation energy units to kilojoules per mole.
inline double EnergyToKJPerMol = EnergyConversionFactor * AvogadroConstant / 1000.0;
/// Conversion factor from simulation energy units to electron volts.
inline double EnergyToEV = EnergyToKelvin / 11604.23;
/// Conversion factor from simulation energy units to kilocalories per mole.
inline double EnergyToKCalPerMol = EnergyConversionFactor * AvogadroConstant / 4184.0;
/// Conversion factor from kilocalories per mole to simulation energy units.
inline double KCalPerMolToEnergy = 4184.0 / (EnergyConversionFactor * AvogadroConstant);

inline double DegreesToRadians = std::numbers::pi / 180.0;
inline double RadiansToDegrees = 180.0 / std::numbers::pi;

inline  std::string unitOfLengthString{"m"};
inline  std::string unitOfEnergyString{"J"};
inline  std::string unitOfEnergyPerMolString{"J.mol⁻¹"};
inline  std::string unitOfMassString{"kg"};
inline  std::string unitOfChargeString{"C/particle"};
inline  std::string unitOfTimeString{"s"};

inline std::string displayedUnitOfLengthString{"Å"};
inline std::string displayedUnitOfEnergyString{"K"};
inline std::string displayedUnitOfEnergyConversionString = "/kʙ";

inline std::string unitOfVelocityString{"m.s⁻¹"};
inline std::string unitOfForceString{"N"};
inline std::string unitOfDiffusionString{"m².s⁻¹"};
inline std::string unitOfAccelerationString{"m².s⁻¹"};
inline std::string unitOfTorqueString{"N.m"};
inline std::string unitOfTemperatureString{"K"};
inline std::string unitOfBoltzmannFactorString{"m².kg.s⁻².K⁻¹"};
inline std::string unitOfPressureString{"Pa"};
inline std::string unitOfVolumeString{"m³"};
inline std::string unitOfDensityString{"kg.m⁻³"};
inline std::string unitOfDynamicViscosityString{"N.s.m⁻²"};
inline std::string unitOfEnthalpyString{"J"};
inline std::string unitOfEnthalpyPerMolString{"J.mol⁻¹"};
inline std::string unitOfPolarizabilityString{"-"};
inline std::string unitOfCoulombPotentialString{"K"};
inline std::string unitOfDielectricConstantString{"s².C².kg⁻¹.m⁻³"};
inline std::string unitOfDipoleMomentString{"C.m"};
inline std::string unitOfDebyeString{"C.m"};
inline std::string unitOfElectricPotentialString{"V"};
inline std::string unitOfElectricFieldString{"V.m⁻¹"};
inline std::string unitOfIsothermalCompressibilityString{"Pa⁻¹"};
inline std::string unitOfHeatCapacityString{"J.mol⁻¹.K⁻¹"};
inline std::string unitOfVolumetricExpansionCoefficientString{"K⁻¹"};

void setUnits(Units::System unit_system)
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
