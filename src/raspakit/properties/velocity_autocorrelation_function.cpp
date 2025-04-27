module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <print>
#include <source_location>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#endif

module property_vacf;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <string>;
import <iostream>;
import <fstream>;
import <sstream>;
import <tuple>;
import <vector>;
import <algorithm>;
import <complex>;
import <format>;
import <numbers>;
import <span>;
import <array>;
import <cmath>;
import <exception>;
import <source_location>;
import <print>;
import <filesystem>;
#endif

import archive;
import double3;
import double4;
import atom;
import simulationbox;
import forcefield;
import component;
import averages;

void PropertyVelocityAutoCorrelationFunction::addSample(size_t currentCycle, const std::vector<Component> &components,
                                                        const std::vector<size_t> &numberOfMoleculesPerComponent,
                                                        std::vector<Molecule> &moleculePositions)
{
  if (currentCycle % sampleEvery != 0uz) return;

  for (size_t currentBuffer = 0; currentBuffer < numberOfBuffersVACF; ++currentBuffer)
  {
    if (countVACF[currentBuffer] == 0)
    {
      size_t molecule_index{0};
      for (size_t i = 0; i != components.size(); ++i)
      {
        originOnsagerVACF[currentBuffer][i] = double3(0.0, 0.0, 0.0);
        for (size_t m = 0; m != numberOfMoleculesPerComponent[i]; ++m)
        {
          double3 value = moleculePositions[molecule_index].velocity;
          originVACF[currentBuffer][molecule_index] = value;
          originOnsagerVACF[currentBuffer][i] += value;
          ++molecule_index;
        }
      }
    }

    if (countVACF[currentBuffer] >= 0)
    {
      size_t index = static_cast<size_t>(countVACF[currentBuffer]);

      for (size_t i = 0; i < numberOfComponents; ++i)
      {
        acfVACF[currentBuffer][i][index] = double4(0.0, 0.0, 0.0, 0.0);
        sumVel[i] = double3(0.0, 0.0, 0.0);
      }

      size_t molecule_index{0};
      for (size_t i = 0; i != components.size(); ++i)
      {
        for (size_t m = 0; m != numberOfMoleculesPerComponent[i]; ++m)
        {
          double3 value = moleculePositions[molecule_index].velocity;
          acfVACF[currentBuffer][i][index].x += value.x * originVACF[currentBuffer][molecule_index].x;
          acfVACF[currentBuffer][i][index].y += value.y * originVACF[currentBuffer][molecule_index].y;
          acfVACF[currentBuffer][i][index].z += value.z * originVACF[currentBuffer][molecule_index].z;
          acfVACF[currentBuffer][i][index].w += value.x * originVACF[currentBuffer][molecule_index].x +
                                                value.y * originVACF[currentBuffer][molecule_index].y +
                                                value.z * originVACF[currentBuffer][molecule_index].z;

          sumVel[i] += value;
          ++molecule_index;
        }
      }

      for (size_t i = 0; i < numberOfComponents; ++i)
      {
        for (size_t j = 0; j < numberOfComponents; ++j)
        {
          acfOnsagerVACF[currentBuffer][i][j][index].x = sumVel[i].x * originOnsagerVACF[currentBuffer][j].x;
          acfOnsagerVACF[currentBuffer][i][j][index].y = sumVel[i].y * originOnsagerVACF[currentBuffer][j].y;
          acfOnsagerVACF[currentBuffer][i][j][index].z = sumVel[i].z * originOnsagerVACF[currentBuffer][j].z;
          acfOnsagerVACF[currentBuffer][i][j][index].w = sumVel[i].x * originOnsagerVACF[currentBuffer][j].x +
                                                         sumVel[i].y * originOnsagerVACF[currentBuffer][j].y +
                                                         sumVel[i].z * originOnsagerVACF[currentBuffer][j].z;
        }
      }
    }
    countVACF[currentBuffer]++;
  }

  // accumulate the vacf
  for (size_t currentBuffer = 0; currentBuffer < numberOfBuffersVACF; ++currentBuffer)
  {
    if (countVACF[currentBuffer] == static_cast<std::make_signed_t<std::size_t>>(bufferLengthVACF))
    {
      for (size_t k = 0; k < numberOfComponents; ++k)
      {
        for (size_t i = 0; i < bufferLengthVACF; ++i)
        {
          accumulatedAcfVACF[k][i].x += acfVACF[currentBuffer][k][i].x;
          accumulatedAcfVACF[k][i].y += acfVACF[currentBuffer][k][i].y;
          accumulatedAcfVACF[k][i].z += acfVACF[currentBuffer][k][i].z;
          accumulatedAcfVACF[k][i].w += acfVACF[currentBuffer][k][i].w;
        }
      }
      for (size_t k = 0; k < numberOfComponents; ++k)
      {
        for (size_t l = 0; l < numberOfComponents; ++l)
        {
          for (size_t i = 0; i < bufferLengthVACF; ++i)
          {
            accumulatedAcfOnsagerVACF[k][l][i].x += acfOnsagerVACF[currentBuffer][k][l][i].x;
            accumulatedAcfOnsagerVACF[k][l][i].y += acfOnsagerVACF[currentBuffer][k][l][i].y;
            accumulatedAcfOnsagerVACF[k][l][i].z += acfOnsagerVACF[currentBuffer][k][l][i].z;
            accumulatedAcfOnsagerVACF[k][l][i].w += acfOnsagerVACF[currentBuffer][k][l][i].w;
          }
        }
      }
      countVACF[currentBuffer] = 0;
      ++countAccumulatedVACF;
    }
  }
}

void PropertyVelocityAutoCorrelationFunction::writeOutput(size_t systemId, const std::vector<Component> &components,
                                                          const std::vector<size_t> &numberOfMoleculesPerComponent,
                                                          double deltaT, size_t currentCycle)
{
  if (currentCycle % writeEvery != 0uz) return;

  std::filesystem::create_directory("vacf");

  for (size_t i = 0; i < components.size(); ++i)
  {
    std::ofstream stream_vacf_self_output(std::format("vacf/vacf_self_{}.s{}.txt", components[i].name, systemId));

    stream_vacf_self_output << std::format("# vacf, number of counts: {}\n", countAccumulatedVACF);
    stream_vacf_self_output << "# column 1: time [ps]\n";
    stream_vacf_self_output << "# column 2: vacf xyz [A^2]\n";
    stream_vacf_self_output << "# column 3: vacf x [A^2]\n";
    stream_vacf_self_output << "# column 4: vacf y [A^2]\n";
    stream_vacf_self_output << "# column 5: vacf z [A^2]\n";
    stream_vacf_self_output << "# column 6: number of samples [-]\n";

    double fac = 1.0 / static_cast<double>(numberOfMoleculesPerComponent[i] * countAccumulatedVACF);

    for (size_t k = 0; k < bufferLengthVACF; ++k)
    {
      stream_vacf_self_output << std::format(
          "{} {} {} {} {} (count: {})\n", static_cast<double>(k * sampleEvery) * deltaT,
          fac * accumulatedAcfVACF[i][k].w, fac * accumulatedAcfVACF[i][k].x, fac * accumulatedAcfVACF[i][k].y,
          fac * accumulatedAcfVACF[i][k].z, countAccumulatedVACF);
    }
  }

  for (size_t i = 0; i < components.size(); ++i)
  {
    double fac = 1.0 / static_cast<double>(numberOfMoleculesPerComponent[i] * countAccumulatedVACF);
    for (size_t j = 0; j < components.size(); ++j)
    {
      std::ofstream stream_vacf_onsager_output(
          std::format("vacf/vacf_onsager_{}_{}.s{}.txt", components[i].name, components[j].name, systemId));

      stream_vacf_onsager_output << std::format("# vacf, number of counts: {}\n", countAccumulatedVACF);
      stream_vacf_onsager_output << "# column 1: time [ps]\n";
      stream_vacf_onsager_output << "# column 2: vacf xyz [A^2]\n";
      stream_vacf_onsager_output << "# column 3: vacf x [A^2]\n";
      stream_vacf_onsager_output << "# column 4: vacf y [A^2]\n";
      stream_vacf_onsager_output << "# column 5: vacf z [A^2]\n";
      stream_vacf_onsager_output << "# column 6: number of samples [-]\n";
      for (size_t k = 0; k < bufferLengthVACF; ++k)
      {
        stream_vacf_onsager_output << std::format(
            "{} {} {} {} {} (count: {})\n", static_cast<double>(k * sampleEvery) * deltaT,
            fac * accumulatedAcfOnsagerVACF[i][j][k].w, fac * accumulatedAcfOnsagerVACF[i][j][k].x,
            fac * accumulatedAcfOnsagerVACF[i][j][k].y, fac * accumulatedAcfOnsagerVACF[i][j][k].z,
            countAccumulatedVACF);
      }
    }
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyVelocityAutoCorrelationFunction &vacf)
{
  archive << vacf.versionNumber;

  archive << vacf.sampleEvery;
  archive << vacf.writeEvery;
  archive << vacf.numberOfComponents;
  archive << vacf.numberOfParticles;

  archive << vacf.numberOfBuffersVACF;
  archive << vacf.bufferLengthVACF;

  archive << vacf.originVACF;
  archive << vacf.originOnsagerVACF;
  archive << vacf.acfVACF;
  archive << vacf.acfOnsagerVACF;

  archive << vacf.accumulatedAcfVACF;
  archive << vacf.accumulatedAcfOnsagerVACF;

  archive << vacf.countVACF;
  archive << vacf.countAccumulatedVACF;
  archive << vacf.sumVel;

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyVelocityAutoCorrelationFunction &vacf)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > vacf.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'PropertyVelocityAutoCorrelationFunction' at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> vacf.sampleEvery;
  archive >> vacf.writeEvery;
  archive >> vacf.numberOfComponents;
  archive >> vacf.numberOfParticles;

  archive >> vacf.numberOfBuffersVACF;
  archive >> vacf.bufferLengthVACF;

  archive >> vacf.originVACF;
  archive >> vacf.originOnsagerVACF;
  archive >> vacf.acfVACF;
  archive >> vacf.acfOnsagerVACF;

  archive >> vacf.accumulatedAcfVACF;
  archive >> vacf.accumulatedAcfOnsagerVACF;

  archive >> vacf.countVACF;
  archive >> vacf.countAccumulatedVACF;
  archive >> vacf.sumVel;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyVelocityAutoCorrelationFunction: Error in binary restart\n"));
  }
#endif

  return archive;
}
