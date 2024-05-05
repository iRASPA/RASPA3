module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <tuple>
#include <vector>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

export module thermostat;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <tuple>;
import <vector>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import archive;


export struct Thermostat
{

  Thermostat(size_t thermostatChainLength):
    thermostatChainLength(thermostatChainLength), 
    thermostatDegreesOfFreedomTranslation(thermostatChainLength),
    thermostatForceTranslation(thermostatChainLength),
    thermostatVelocityTranslation(thermostatChainLength),
    thermostatPositionTranslation(thermostatChainLength),
    thermostatMassTranslation(thermostatChainLength),
    thermostatDegreesOfFreedomRotation(thermostatChainLength),
    thermostatForceRotation(thermostatChainLength),
    thermostatVelocityRotation(thermostatChainLength),
    thermostatPositionRotation(thermostatChainLength),
    thermostatMassRotation(thermostatChainLength)
  {
  }


  size_t thermostatChainLength;
  double timeScaleParameterThermostat;

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
};
