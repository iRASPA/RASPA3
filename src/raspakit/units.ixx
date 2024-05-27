module;

#ifdef USE_LEGACY_HEADERS
#include <numbers>
#include <string>
#include <sstream>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

export module units;

#ifndef USE_LEGACY_HEADERS
import <numbers>;
import <string>;
import <sstream>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import hdf5;

  export namespace Units
  {
  // Fundamental constants of nature
  constexpr double Angstrom = 1e-10;
  constexpr double AngstromCubed = 1e-30;
  constexpr double NanoSecond = 1e-9;
  constexpr double PicoSecond = 1e-12;
  constexpr double FemtoSecond = 1e-15;
  constexpr double AtomicMassUnit = 1.6605402e-27;        // kg
  constexpr double ElectronicChargeUnit = 1.60217733e-19; // C/particle
  constexpr double MolarGasConstant = 8.314464919;        // J mol^-1 K^-1
  constexpr double BoltzmannConstant = 1.380650324e-23;   // J K^-1
  constexpr double SpeedOfLight = 299792458.0;            // m s^-1
  constexpr double AvogadroConstant = 6.0221419947e23;    // mol^-1
  constexpr double ElectricConstant = 8.8541878176e-12;   // C^2/(N.m^2)
  constexpr double PlanckConstant = 6.6260687652e-34;     // J.s
  constexpr double DielectricConstantVacuum = 8.8541878176e-12; // F/m
  constexpr double Debye = 3.335640952e-30;                     // C.m

  // Set of basis units
  constexpr double LengthUnit = Angstrom;
  constexpr double TimeUnit = PicoSecond;
  constexpr double MassUnit = AtomicMassUnit;
  constexpr double ChargeUnit = ElectronicChargeUnit;

  // Derived units and their conversion factors
  constexpr double MassConversionFactor = AtomicMassUnit; // kg
  constexpr double VolumeConversionFactor = LengthUnit * LengthUnit * LengthUnit; // m^-3
  constexpr double DensityConversionfactor = MassConversionFactor / VolumeConversionFactor; // kg/m^3
  constexpr double EnergyConversionFactor = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit); // J
  constexpr double ForceConversionFactor = EnergyConversionFactor / LengthUnit; // N (=m/s^2)
  constexpr double PressureConversionFactor = MassUnit / (LengthUnit * TimeUnit * TimeUnit); // Pa(=N/m^2)
  constexpr double PositionConversionFactor = LengthUnit; // m
  constexpr double TimeConversionFactor = TimeUnit; // s
  constexpr double DiffusionConversionFactor = LengthUnit * LengthUnit / TimeUnit; // m^2/s
  constexpr double VelocityConversionFactor = LengthUnit / TimeUnit; // m/s
  constexpr double AccelerationConversionFactor = LengthUnit * LengthUnit / TimeUnit; // m^2/s
  constexpr double DipoleMomentConversionFactor = ChargeUnit * LengthUnit; // C.m
  constexpr double DebyeConversionFactor = ChargeUnit * LengthUnit / Debye; // C.m
  constexpr double ElectricPotentialConversionFactor = (EnergyConversionFactor / ChargeUnit); // V
  constexpr double PolarizilibityConversionFactor = ChargeUnit * ChargeUnit * LengthUnit * LengthUnit / EnergyConversionFactor;
  constexpr double ElectricFieldConversionFactor = ElectricPotentialConversionFactor / LengthUnit; // V
  constexpr double KB = BoltzmannConstant / EnergyConversionFactor;
  constexpr double CoulombicConversionFactor = ChargeUnit * ChargeUnit / (4.0 * std::numbers::pi * DielectricConstantVacuum * LengthUnit * EnergyConversionFactor);
  constexpr double DielectricConstantConversionFactor = TimeUnit * TimeUnit * ChargeUnit * ChargeUnit / (AtomicMassUnit * LengthUnit * LengthUnit * LengthUnit);
  constexpr double IsothermalCompressibilityConversionFactor = 1.0 / PressureConversionFactor;
  constexpr double HeatCapacityConversionFactor = (EnergyConversionFactor * AvogadroConstant); // J/mol/K
  constexpr double VolumetricExpansionCoefficientConversionFactor = 1.0; // 1/K
  constexpr double FeymannHibbsConversionFactor = (PlanckConstant * PlanckConstant / (2.0 * std::numbers::pi)) / 
                    (24.0 * AtomicMassUnit * BoltzmannConstant * LengthUnit * LengthUnit);

  // 1 calorie (International Table) = 4.1868 J
  // 1 calorie (thermochemical) = 4.184 J
  // 1 calorie (15C) = 4.1855 J
  constexpr double CalToJoule = 4.184;
  constexpr double JouleToCal = 1.0 / CalToJoule;
  constexpr double EvToKJPerMol = 96.48534;
  constexpr double EvToKelvin = EvToKJPerMol * 1000.0 / MolarGasConstant;
  
  constexpr double EnergyToKelvin = EnergyConversionFactor * AvogadroConstant / MolarGasConstant;
  constexpr double KelvinToEnergy = MolarGasConstant / (EnergyConversionFactor * AvogadroConstant);
  constexpr double EnergyToKJPerMol = EnergyConversionFactor * AvogadroConstant / 1000.0;
  constexpr double EnergyToEV = EnergyToKelvin / 11604.23;
  constexpr double EnergyToKCalPerMol = EnergyConversionFactor * AvogadroConstant / 4184.0;
  constexpr double KCalPerMolToEnergy = 4184.0 / (EnergyConversionFactor * AvogadroConstant);
  
  constexpr double PaToAtm = (1.0 / 101325.0);
  constexpr double PatoBar = 0.00001;
  constexpr double PaToTorr = 0.0075;
  constexpr double KPaToAtm = 1.0 / 101.325;
  constexpr double KPaToBar = 0.01;
  constexpr double KPaToTorr = 7.5;
  constexpr double BarToTorr = 750.062;
  constexpr double BarToAtm = 0.9869;
  constexpr double BarTokPa = 100.0;
  constexpr double AtmToTorr = 760;
  constexpr double AtmToBar = 1.0 / 0.9869;
  constexpr double AtmTokPa = 101.325;
  constexpr double AtmToPa = 101325.0;

  
  std::string printStatus();
  void logStatus(HDF5Writer& hdf5);
};
