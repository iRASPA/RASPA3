export module units;

import <numbers>;
import <string>;
import <sstream>;

import print;

export namespace Units
{
  // Fundamental constants of nature
  inline const double Angstrom = 1e-10;
  inline const double AngstromCubed = 1e-30;
  inline const double NanoSecond = 1e-9;
  inline const double PicoSecond = 1e-12;
  inline const double FemtoSecond = 1e-15;
  inline const double AtomicMassUnit = 1.6605402e-27;        // kg
  inline const double ElectronicChargeUnit = 1.60217733e-19; // C/particle
  inline const double MolarGasConstant = 8.314464919;        // J mol^-1 K^-1
  inline const double BoltzmannConstant = 1.380650324e-23;   // J K^-1
  inline const double SpeedOfLight = 299792458.0;            // m s^-1
  inline const double AvogadroConstant = 6.0221419947e23;    // mol^-1
  inline const double ElectricConstant = 8.8541878176e-12;   // C^2/(N.m^2)
  inline const double PlanckConstant = 6.6260687652e-34;     // J.s
  inline const double DielectricConstantVacuum = 8.8541878176e-12; // F/m
  inline const double Debye = 3.335640952e-30;                     // C.m

  // Set of basis units
  inline const double LengthUnit = Angstrom;
  inline const double TimeUnit = PicoSecond;
  inline const double MassUnit = AtomicMassUnit;
  inline const double ChargeUnit = ElectronicChargeUnit;

  // Derived units and their conversion factors
  inline const double MassConversionFactor = AtomicMassUnit; // kg
  inline const double VolumeConversionFactor = LengthUnit * LengthUnit * LengthUnit; // m^-3
  inline const double DensityConversionfactor = MassConversionFactor / VolumeConversionFactor; // kg/m^3
  inline const double EnergyConversionFactor = MassUnit * LengthUnit * LengthUnit / (TimeUnit * TimeUnit); // J
  inline const double ForceConversionFactor = EnergyConversionFactor / LengthUnit; // N (=m/s^2)
  inline const double PressureConversionFactor = MassUnit / (LengthUnit * TimeUnit * TimeUnit); // Pa(=N/m^2)
  inline const double PositionConversionFactor = LengthUnit; // m
  inline const double TimeConversionFactor = TimeUnit; // s
  inline const double DiffusionConversionFactor = LengthUnit * LengthUnit / TimeUnit; // m^2/s
  inline const double VelocityConversionFactor = LengthUnit / TimeUnit; // m/s
  inline const double AccelerationConversionFactor = LengthUnit * LengthUnit / TimeUnit; // m^2/s
  inline const double DipoleMomentConversionFactor = ChargeUnit * LengthUnit; // C.m
  inline const double DebyeConversionFactor = ChargeUnit * LengthUnit / Debye; // C.m
  inline const double ElectricPotentialConversionFactor = (EnergyConversionFactor / ChargeUnit); // V
  inline const double PolarizilibityConversionFactor = ChargeUnit * ChargeUnit * LengthUnit * LengthUnit / EnergyConversionFactor;
  inline const double ElectricFieldConversionFactor = ElectricPotentialConversionFactor / LengthUnit; // V
  inline const double KB = BoltzmannConstant / EnergyConversionFactor;
  inline const double CoulombicConversionFactor = ChargeUnit * ChargeUnit / (4.0 * std::numbers::pi * DielectricConstantVacuum * LengthUnit * EnergyConversionFactor);
  inline const double DielectricConstantConversionFactor = TimeUnit * TimeUnit * ChargeUnit * ChargeUnit / (AtomicMassUnit * LengthUnit * LengthUnit * LengthUnit);
  inline const double IsothermalCompressibilityConversionFactor = 1.0 / PressureConversionFactor;
  inline const double HeatCapacityConversionFactor = (EnergyConversionFactor * AvogadroConstant); // J/mol/K
  inline const double VolumetricExpansionCoefficientConversionFactor = 1.0; // 1/K
  inline const double FeymannHibbsConversionFactor = (PlanckConstant * PlanckConstant / (2.0 * std::numbers::pi)) / 
                    (24.0 * AtomicMassUnit * BoltzmannConstant * LengthUnit * LengthUnit);

  // 1 calorie (International Table) = 4.1868 J
  // 1 calorie (thermochemical) = 4.184 J
  // 1 calorie (15C) = 4.1855 J
  inline const double CalToJoule = 4.184;
  inline const double JouleToCal = 1.0 / CalToJoule;
  inline const double EvToKJPerMol = 96.48534;
  inline const double EvToKelvin = EvToKJPerMol * 1000.0 / MolarGasConstant;

  inline const double EnergyToKelvin = EnergyConversionFactor * AvogadroConstant / MolarGasConstant;
  inline const double KelvinToEnergy = MolarGasConstant / (EnergyConversionFactor * AvogadroConstant);
  inline const double EnergyToKJPerMol = EnergyConversionFactor * AvogadroConstant / 1000.0;
  inline const double EnergyToEV = EnergyToKelvin / 11604.23;
  inline const double EnergyToKCalPerMol = EnergyConversionFactor * AvogadroConstant / 4184.0;
  inline const double KCalPerMolToEnergy = 4184.0 / (EnergyConversionFactor * AvogadroConstant);

  inline const double PaToAtm = (1.0 / 101325.0);
  inline const double PatoBar = 0.00001;
  inline const double PaToTorr = 0.0075;
  inline const double KPaToAtm = 1.0 / 101.325;
  inline const double KPaToBar = 0.01;
  inline const double KPaToTorr = 7.5;
  inline const double BarToTorr = 750.062;
  inline const double BarToAtm = 0.9869;
  inline const double BarTokPa = 100.0;
  inline const double AtmToTorr = 760;
  inline const double AtmToBar = 1.0 / 0.9869;
  inline const double AtmTokPa = 101.325;
  inline const double AtmToPa = 101325.0;

  
  std::string printStatus();
};
