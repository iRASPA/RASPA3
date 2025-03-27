module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#endif

export module property_vacf;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <array>;
import <optional>;
import <cmath>;
import <string>;
import <span>;
import <tuple>;
#endif

import archive;
import double3;
import double4;
import atom;
import molecule;
import molecule;
import simulationbox;
import forcefield;
import component;

// Computes Velocity Auto-Correlation Function (VASF)

export struct PropertyVelocityAutoCorrelationFunction
{
  PropertyVelocityAutoCorrelationFunction() {};

  PropertyVelocityAutoCorrelationFunction(size_t numberOfComponents, size_t numberOfParticles,
                                          size_t numberOfBuffersVACF, size_t bufferLengthVACF, size_t sampleEvery,
                                          size_t writeEvery)
      : numberOfComponents(numberOfComponents),
        numberOfParticles(numberOfParticles),
        numberOfBuffersVACF(numberOfBuffersVACF),
        bufferLengthVACF(bufferLengthVACF),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        originVACF(numberOfBuffersVACF, std::vector<double3>(numberOfParticles, double3(0.0, 0.0, 0.0))),
        originOnsagerVACF(numberOfBuffersVACF, std::vector<double3>(numberOfComponents, double3(0.0, 0.0, 0.0))),
        acfVACF(numberOfBuffersVACF,
                std::vector<std::vector<double4>>(numberOfComponents,
                                                  std::vector<double4>(bufferLengthVACF, double4(0.0, 0.0, 0.0, 0.0)))),
        acfOnsagerVACF(
            numberOfBuffersVACF,
            std::vector<std::vector<std::vector<double4>>>(
                numberOfComponents,
                std::vector<std::vector<double4>>(
                    numberOfComponents, std::vector<double4>(bufferLengthVACF, double4(0.0, 0.0, 0.0, 0.0))))),
        accumulatedAcfVACF(numberOfComponents, std::vector<double4>(bufferLengthVACF, double4(0.0, 0.0, 0.0, 0.0))),
        accumulatedAcfOnsagerVACF(
            numberOfBuffersVACF,
            std::vector<std::vector<double4>>(numberOfComponents,
                                              std::vector<double4>(bufferLengthVACF, double4(0.0, 0.0, 0.0, 0.0)))),
        countVACF(numberOfBuffersVACF),
        countAccumulatedVACF(0),
        sumVel(numberOfComponents)
  {
    // trick to space the origins evenly (see for example Rapaport 2004)
    for (size_t currentBuffer = 0; currentBuffer < numberOfBuffersVACF; ++currentBuffer)
    {
      countVACF[currentBuffer] = -static_cast<int64_t>(currentBuffer * bufferLengthVACF / numberOfBuffersVACF);
    }
  }

  uint64_t versionNumber{1};

  size_t numberOfComponents;
  size_t numberOfParticles;

  size_t numberOfBuffersVACF;
  size_t bufferLengthVACF;

  size_t sampleEvery;
  size_t writeEvery;

  std::vector<std::vector<double3>> originVACF;
  std::vector<std::vector<double3>> originOnsagerVACF;
  std::vector<std::vector<std::vector<double4>>> acfVACF;
  std::vector<std::vector<std::vector<std::vector<double4>>>> acfOnsagerVACF;

  std::vector<std::vector<double4>> accumulatedAcfVACF;
  std::vector<std::vector<std::vector<double4>>> accumulatedAcfOnsagerVACF;

  std::vector<int64_t> countVACF;
  size_t countAccumulatedVACF;
  std::vector<double3> sumVel;

  void addSample(size_t currentCycle, const std::vector<Component> &components,
                 const std::vector<size_t> &numberOfMoleculesPerComponent, std::vector<Molecule> &molecules);
  void writeOutput(size_t systemId, const std::vector<Component> &components,
                   const std::vector<size_t> &numberOfMoleculesPerComponent, double deltaT, size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyVelocityAutoCorrelationFunction &msd);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive,
                                            PropertyVelocityAutoCorrelationFunction &msd);
};
