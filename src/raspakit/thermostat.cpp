module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <format>
#include <print>
#include <source_location>
#include <tuple>
#include <vector>
#endif

module thermostat;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <tuple>;
import <cmath>;
import <print>;
import <source_location>;
#endif

import archive;
import units;
import randomnumbers;

Thermostat::Thermostat(double temperature, size_t thermostatChainLength, size_t numberOfYoshidaSuzukiSteps,
                       double deltaT, size_t translationalDegreesOfFreedom, size_t rotationalDgreesOfFreedom)
    : temperature(temperature),
      thermostatChainLength(thermostatChainLength),
      numberOfYoshidaSuzukiSteps(numberOfYoshidaSuzukiSteps),
      deltaT(deltaT),
      translationalDegreesOfFreedom(translationalDegreesOfFreedom),
      rotationalDgreesOfFreedom(rotationalDgreesOfFreedom),
      thermostatDegreesOfFreedomTranslation(thermostatChainLength),
      thermostatForceTranslation(thermostatChainLength),
      thermostatVelocityTranslation(thermostatChainLength),
      thermostatPositionTranslation(thermostatChainLength),
      thermostatMassTranslation(thermostatChainLength),
      thermostatDegreesOfFreedomRotation(thermostatChainLength),
      thermostatForceRotation(thermostatChainLength),
      thermostatVelocityRotation(thermostatChainLength),
      thermostatPositionRotation(thermostatChainLength),
      thermostatMassRotation(thermostatChainLength),
      w(numberOfYoshidaSuzukiSteps)
{
  switch (numberOfYoshidaSuzukiSteps)
  {
    case 1:
      w[0] = 1.0;
      break;
    case 3:
      w[0] = 1.0 / (2.0 - std::pow(2.0, 1.0 / 3.0));
      w[1] = 1.0 - 2.0 * w[0];
      w[2] = w[0];
      break;
    case 5:
      w[0] = 1.0 / (4.0 - std::pow(4.0, 1.0 / 3.0));
      w[1] = w[0];
      w[2] = 1.0 - 4.0 * w[0];
      w[3] = w[0];
      w[4] = w[0];
      break;
    case 7:
      w[0] = 0.784513610477560;
      w[1] = 0.235573213359357;
      w[2] = -1.17767998417887;
      w[3] = 1.0 - 2.0 * (w[0] + w[1] + w[2]);
      w[4] = -1.17767998417887;
      w[5] = 0.235573213359357;
      w[6] = 0.784513610477560;
      break;
    case 9:
      w[0] = 0.192;
      w[1] = 0.554910818409783619692725006662999;
      w[2] = 0.124659619941888644216504240951585;
      w[3] = -0.843182063596933505315033808282941;
      w[4] = 1.0 - 2.0 * (w[0] + w[1] + w[2] + w[3]);
      w[5] = -0.843182063596933505315033808282941;
      w[6] = 0.124659619941888644216504240951585;
      w[7] = 0.554910818409783619692725006662999;
      w[8] = 0.192;
      break;
    default:
      throw std::runtime_error(std::format("Error: Yoshida-Suzuki-steps should be: 1,3,5,7 or 9\n"));
      break;
  }
}

void Thermostat::initialize(RandomNumber &random)
{
  thermostatDegreesOfFreedomTranslation[0] =
      static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint) * Units::KB *
      temperature;
  for (size_t i = 1; i != thermostatChainLength; ++i)
  {
    thermostatDegreesOfFreedomTranslation[i] = Units::KB * temperature;
  }

  for (size_t i = 0; i != thermostatChainLength; ++i)
  {
    thermostatMassTranslation[i] =
        static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint) *
        timeScaleParameterThermostat * timeScaleParameterThermostat;
  }

  thermostatDegreesOfFreedomRotation[0] = static_cast<double>(rotationalDgreesOfFreedom) * Units::KB * temperature;
  for (size_t i = 1; i != thermostatChainLength; ++i)
  {
    thermostatDegreesOfFreedomRotation[i] = Units::KB * temperature;
  }

  for (size_t i = 0; i != thermostatChainLength; ++i)
  {
    thermostatMassRotation[i] =
        static_cast<double>(rotationalDgreesOfFreedom) * timeScaleParameterThermostat * timeScaleParameterThermostat;
  }

  for (size_t i = 0; i != thermostatChainLength; ++i)
  {
    thermostatVelocityTranslation[i] =
        random.Gaussian() *
        std::sqrt(static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint) /
                  thermostatMassTranslation[i]);

    thermostatVelocityRotation[i] =
        random.Gaussian() * std::sqrt(static_cast<double>(rotationalDgreesOfFreedom) / thermostatMassTranslation[i]);
  }
}

std::pair<double, double> Thermostat::NoseHooverNVT(double UKineticTranslation, double UKineticRotation)
{
  double AA;

  size_t M = thermostatChainLength;
  size_t nc = numberOfRespaSteps;
  size_t nyosh = numberOfYoshidaSuzukiSteps;

  double scale_translation = 1.0;

  if (translationalDegreesOfFreedom > 0)
  {
    thermostatForceTranslation[0] =
        (2.0 * UKineticTranslation - thermostatDegreesOfFreedomTranslation[0]) / thermostatMassTranslation[0];

    for (size_t i = 0; i < nc; i++)
    {
      for (size_t j = 0; j < nyosh; j++)
      {
        thermostatVelocityTranslation[M - 1] +=
            thermostatForceTranslation[M - 1] * w[j] * deltaT / (4.0 * static_cast<double>(nc));

        for (size_t k = 0; k < M - 1; k++)
        {
          AA = exp(-(w[j] * deltaT / (8.0 * static_cast<double>(nc))) * thermostatVelocityTranslation[M - k - 1]);
          thermostatVelocityTranslation[M - k - 2] =
              thermostatVelocityTranslation[M - k - 2] * AA * AA +
              thermostatForceTranslation[M - k - 2] * AA * w[j] * deltaT / (4.0 * static_cast<double>(nc));
        }

        AA = exp(-(w[j] * deltaT / (2.0 * static_cast<double>(nc))) * thermostatVelocityTranslation[0]);
        scale_translation *= AA;
        thermostatForceTranslation[0] = (scale_translation * scale_translation * 2.0 * UKineticTranslation -
                                         thermostatDegreesOfFreedomTranslation[0]) /
                                        thermostatMassTranslation[0];

        for (size_t k = 0; k < M; k++)
        {
          thermostatPositionTranslation[k] +=
              thermostatVelocityTranslation[k] * w[j] * deltaT / (2.0 * static_cast<double>(nc));
        }

        for (size_t k = 0; k < M - 1; k++)
        {
          AA = exp(-(w[j] * deltaT / (8.0 * static_cast<double>(nc))) * thermostatVelocityTranslation[k + 1]);
          thermostatVelocityTranslation[k] =
              thermostatVelocityTranslation[k] * AA * AA +
              thermostatForceTranslation[k] * AA * (w[j] * deltaT / (4.0 * static_cast<double>(nc)));

          thermostatForceTranslation[k + 1] =
              (thermostatMassTranslation[k] * thermostatVelocityTranslation[k] * thermostatVelocityTranslation[k] -
               thermostatDegreesOfFreedomTranslation[k + 1]) /
              thermostatMassTranslation[k + 1];
        }
        thermostatVelocityTranslation[M - 1] +=
            thermostatForceTranslation[M - 1] * w[j] * deltaT / (4.0 * static_cast<double>(nc));
      }
    }
  }

  double scale_rotation = 1.0;

  if (rotationalDgreesOfFreedom > 0)
  {
    thermostatForceRotation[0] =
        (2.0 * UKineticRotation - thermostatDegreesOfFreedomRotation[0]) / thermostatMassRotation[0];

    for (size_t i = 0; i < nc; i++)
    {
      for (size_t j = 0; j < nyosh; j++)
      {
        thermostatVelocityRotation[M - 1] +=
            thermostatForceRotation[M - 1] * w[j] * deltaT / (4.0 * static_cast<double>(nc));

        for (size_t k = 0; k < M - 1; k++)
        {
          AA = exp(-(w[j] * deltaT / (8.0 * static_cast<double>(nc))) * thermostatVelocityRotation[M - k - 1]);
          thermostatVelocityRotation[M - k - 2] =
              thermostatVelocityRotation[M - k - 2] * AA * AA +
              thermostatForceRotation[M - k - 2] * AA * w[j] * deltaT / (4.0 * static_cast<double>(nc));
        }

        AA = exp(-(w[j] * deltaT / (2.0 * static_cast<double>(nc))) * thermostatVelocityRotation[0]);
        scale_rotation *= AA;
        thermostatForceRotation[0] =
            (scale_rotation * scale_rotation * 2.0 * UKineticRotation - thermostatDegreesOfFreedomRotation[0]) /
            thermostatMassRotation[0];

        for (size_t k = 0; k < M; k++)
        {
          thermostatPositionRotation[k] +=
              thermostatVelocityRotation[k] * w[j] * deltaT / (2.0 * static_cast<double>(nc));
        }

        for (size_t k = 0; k < M - 1; k++)
        {
          AA = exp(-(w[j] * deltaT / (8.0 * static_cast<double>(nc))) * thermostatVelocityRotation[k + 1]);
          thermostatVelocityRotation[k] =
              thermostatVelocityRotation[k] * AA * AA +
              thermostatForceRotation[k] * AA * (w[j] * deltaT / (4.0 * static_cast<double>(nc)));

          thermostatForceRotation[k + 1] =
              (thermostatMassRotation[k] * thermostatVelocityRotation[k] * thermostatVelocityRotation[k] -
               thermostatDegreesOfFreedomRotation[k + 1]) /
              thermostatMassRotation[k + 1];
        }
        thermostatVelocityRotation[M - 1] +=
            thermostatForceRotation[M - 1] * w[j] * deltaT / (4.0 * static_cast<double>(nc));
      }
    }
  }

  return {scale_translation, scale_rotation};
}

double Thermostat::getEnergy()
{
  double energy{};

  if (translationalDegreesOfFreedom > 0)
  {
    energy += 0.5 * thermostatMassTranslation[0] * thermostatVelocityTranslation[0] * thermostatVelocityTranslation[0] +
              static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint) * Units::KB *
                  temperature * thermostatPositionTranslation[0];

    for (size_t i = 1; i < thermostatChainLength; i++)
    {
      energy +=
          0.5 * thermostatMassTranslation[i] * thermostatVelocityTranslation[i] * thermostatVelocityTranslation[i] +
          Units::KB * temperature * thermostatPositionTranslation[i];
    }
  }

  if (rotationalDgreesOfFreedom > 0)
  {
    energy += 0.5 * thermostatMassRotation[0] * thermostatVelocityRotation[0] * thermostatVelocityRotation[0] +
              static_cast<double>(rotationalDgreesOfFreedom) * Units::KB * temperature * thermostatPositionRotation[0];

    for (size_t i = 1; i < thermostatChainLength; i++)
    {
      energy += 0.5 * thermostatMassRotation[i] * thermostatVelocityRotation[i] * thermostatVelocityRotation[i] +
                Units::KB * temperature * thermostatPositionRotation[i];
    }
  }

  return energy;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Thermostat &t)
{
  archive << t.versionNumber;

  archive << t.temperature;
  archive << t.thermostatChainLength;
  archive << t.timeScaleParameterThermostat;
  archive << t.numberOfRespaSteps;
  archive << t.numberOfYoshidaSuzukiSteps;

  archive << t.thermostatDegreesOfFreedomTranslation;
  archive << t.thermostatForceTranslation;
  archive << t.thermostatVelocityTranslation;
  archive << t.thermostatPositionTranslation;
  archive << t.thermostatMassTranslation;

  archive << t.thermostatDegreesOfFreedomRotation;
  archive << t.thermostatForceRotation;
  archive << t.thermostatVelocityRotation;
  archive << t.thermostatPositionRotation;
  archive << t.thermostatMassRotation;

  archive << t.w;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Thermostat &t)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > t.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'EquationOfState' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> t.temperature;
  archive >> t.thermostatChainLength;
  archive >> t.timeScaleParameterThermostat;
  archive >> t.numberOfRespaSteps;
  archive >> t.numberOfYoshidaSuzukiSteps;

  archive >> t.thermostatDegreesOfFreedomTranslation;
  archive >> t.thermostatForceTranslation;
  archive >> t.thermostatVelocityTranslation;
  archive >> t.thermostatPositionTranslation;
  archive >> t.thermostatMassTranslation;

  archive >> t.thermostatDegreesOfFreedomRotation;
  archive >> t.thermostatForceRotation;
  archive >> t.thermostatVelocityRotation;
  archive >> t.thermostatPositionRotation;
  archive >> t.thermostatMassRotation;

  archive >> t.w;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Thermostat: Error in binary restart\n"));
  }
#endif

  return archive;
}
