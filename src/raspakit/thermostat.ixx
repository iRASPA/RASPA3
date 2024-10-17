module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <print>
#include <tuple>
#include <vector>
#endif

export module thermostat;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <tuple>;
import <vector>;
import <print>;
#endif

import archive;
import units;
import randomnumbers;
import molecule;

export struct Thermostat
{
  Thermostat() {}
  Thermostat(double temperature, size_t thermostatChainLength, size_t numberOfYoshidaSuzukiSteps, double deltaT,
             size_t translationalDegreesOfFreedom, size_t rotationalDgreesOfFreedom);

  uint64_t versionNumber{1};

  double temperature;
  size_t thermostatChainLength;
  double timeScaleParameterThermostat{0.15};
  size_t numberOfRespaSteps{5};
  size_t numberOfYoshidaSuzukiSteps{5};
  double deltaT{};

  size_t translationalCenterOfMassConstraint{};
  size_t translationalDegreesOfFreedom;
  size_t rotationalDgreesOfFreedom;

  std::vector<double> thermostatDegreesOfFreedomTranslation;
  std::vector<double> thermostatForceTranslation;
  std::vector<double> thermostatVelocityTranslation;
  std::vector<double> thermostatPositionTranslation;
  std::vector<double> thermostatMassTranslation;

  std::vector<double> thermostatDegreesOfFreedomRotation;
  std::vector<double> thermostatForceRotation;
  std::vector<double> thermostatVelocityRotation;
  std::vector<double> thermostatPositionRotation;
  std::vector<double> thermostatMassRotation;

  std::vector<double> w;

  void initialize(RandomNumber &random);

  std::pair<double, double> NoseHooverNVT(double UKineticTranslation, double UKineticRotation);

  double getEnergy();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Thermostat &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Thermostat &s);
};
