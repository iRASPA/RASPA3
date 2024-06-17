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

export struct Thermostat
{
  Thermostat(double temperature, size_t thermostatChainLength, size_t numberOfYoshidaSuzukiSteps);

  void initialize(RandomNumber &random, size_t translation_degrees_of_freedom, size_t rotational_degrees_of_freedom);
  void NoseHooverNVT();

  double temperature;
  size_t thermostatChainLength;
  double timeScaleParameterThermostat{0.15};
  size_t numberOfYoshidaSuzukiSteps{5};

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
};
