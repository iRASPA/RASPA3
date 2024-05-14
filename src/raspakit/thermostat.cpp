module;

#ifdef USE_LEGACY_HEADERS
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <vector>
#include <tuple>
#include <cmath>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

module thermostat;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <tuple>;
import <cmath>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import archive;
import units;
import randomnumbers;

Thermostat::Thermostat(double temperature, size_t thermostatChainLength, size_t numberOfYoshidaSuzukiSteps):
    temperature(temperature),
    thermostatChainLength(thermostatChainLength),
    numberOfYoshidaSuzukiSteps(numberOfYoshidaSuzukiSteps),
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
    switch(numberOfYoshidaSuzukiSteps)
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


void Thermostat::initialize(RandomNumber &random, size_t translation_degrees_of_freedom, size_t rotational_degrees_of_freedom)
{
  thermostatDegreesOfFreedomTranslation[0] = static_cast<double>(translation_degrees_of_freedom) * Units::KB * temperature;
  for(size_t i = 1; i != thermostatChainLength; ++i)
  {
    thermostatDegreesOfFreedomTranslation[i] = Units::KB * temperature;
  }

  for(size_t i = 0; i != thermostatChainLength; ++i)
  {
    thermostatMassTranslation[i] = static_cast<double>(translation_degrees_of_freedom) * timeScaleParameterThermostat * timeScaleParameterThermostat;
  }

  thermostatDegreesOfFreedomRotation[0] = static_cast<double>(rotational_degrees_of_freedom) * Units::KB * temperature;
  for(size_t i = 1; i != thermostatChainLength; ++i)
  {
    thermostatDegreesOfFreedomRotation[i] = Units::KB * temperature;
  }

  for(size_t i = 0; i != thermostatChainLength; ++i)
  {
    thermostatMassRotation[i] = static_cast<double>(rotational_degrees_of_freedom) * timeScaleParameterThermostat * timeScaleParameterThermostat;
  }


  for(size_t i = 0; i != thermostatChainLength; ++i)
  {
    thermostatVelocityTranslation[i] = random.Gaussian() * std::sqrt(static_cast<double>(translation_degrees_of_freedom) / thermostatMassTranslation[i]);
    thermostatVelocityRotation[i] = random.Gaussian() * std::sqrt(static_cast<double>(rotational_degrees_of_freedom) / thermostatMassTranslation[i]);
  }
}

void Thermostat::NoseHooverNVT()
{
  /*
int i,j,k,l;
  int A,M,nc,nyosh,Type;
  REAL AA;
  REAL ScaleTranslationAdsorbates;
  REAL UKineticTranslationAdsorbates;

  M=therm_baro_stats.ThermostatChainLength;
  nc=therm_baro_stats.NumberOfRespaSteps;
  nyosh=therm_baro_stats.NumberOfYoshidaSuzukiSteps;

  ScaleTranslationAdsorbates=1.0;

  UKineticTranslationAdsorbates=2.0*GetTranslationKineticEnergyAdsorbates();

  ThermostatForceTranslationAdsorbates[CurrentSystem][0]=
        (UKineticTranslationAdsorbates-ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][0])/
                ThermostatMassTranslationAdsorbates[CurrentSystem][0];

  for(i=0;i<nc;i++)
  {
    for(j=0;j<nyosh;j++)
    {
      ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-1]+=ThermostatForceTranslationAdsorbates[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);

      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-k-1]);
        ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-k-2]=ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-k-2]*SQR(AA)+
                                       ThermostatForceTranslationAdsorbates[CurrentSystem][M-k-2]*AA*w[j]*DeltaT/(4.0*nc);
      }

      AA=exp(-(w[j]*DeltaT/(2.0*nc))*ThermostatVelocityTranslationAdsorbates[CurrentSystem][0]);
      ScaleTranslationAdsorbates*=AA;
      ThermostatForceTranslationAdsorbates[CurrentSystem][0]=(SQR(ScaleTranslationAdsorbates)*UKineticTranslationAdsorbates-
              ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][0])/
                                          ThermostatMassTranslationAdsorbates[CurrentSystem][0];

      for(k=0;k<M;k++)
        ThermostatPositionTranslationAdsorbates[CurrentSystem][k]+=ThermostatVelocityTranslationAdsorbates[CurrentSystem][k]*w[j]*DeltaT/(2.0*nc);

      for(k=0;k<M-1;k++)
      {
        AA=exp(-(w[j]*DeltaT/(8.0*nc))*ThermostatVelocityTranslationAdsorbates[CurrentSystem][k+1]);
        ThermostatVelocityTranslationAdsorbates[CurrentSystem][k]=ThermostatVelocityTranslationAdsorbates[CurrentSystem][k]*SQR(AA)+
                                 ThermostatForceTranslationAdsorbates[CurrentSystem][k]*AA*(w[j]*DeltaT/(4.0*nc));

        ThermostatForceTranslationAdsorbates[CurrentSystem][k+1]=
                (ThermostatMassTranslationAdsorbates[CurrentSystem][k]*SQR(ThermostatVelocityTranslationAdsorbates[CurrentSystem][k])-
                ThermostatDegreesOfFreedomTranslationAdsorbates[CurrentSystem][k+1])/ThermostatMassTranslationAdsorbates[CurrentSystem][k+1];
      }
      ThermostatVelocityTranslationAdsorbates[CurrentSystem][M-1]+=ThermostatForceTranslationAdsorbates[CurrentSystem][M-1]*w[j]*DeltaT/(4.0*nc);
    }
  }
  */
}

