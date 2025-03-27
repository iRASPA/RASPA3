module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#endif

export module property_msd;

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

// Computes Mean Squared Displacement (MSD) using order-N algorithm

export struct PropertyMeanSquaredDisplacement
{
  PropertyMeanSquaredDisplacement() {};

  PropertyMeanSquaredDisplacement(size_t numberOfComponents, size_t numberOfParticles, size_t sampleEvery,
                                  size_t writeEvery, size_t numberOfBlockElementsMSD)
      : sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        numberOfComponents(numberOfComponents),
        numberOfParticles(numberOfParticles),
        numberOfBlockElementsMSD(numberOfBlockElementsMSD),
        maxNumberOfBlocksMSD(1),
        blockLengthMSD(maxNumberOfBlocksMSD),
        msdSelfCount(maxNumberOfBlocksMSD, std::vector<std::vector<size_t>>(
                                               numberOfComponents, std::vector<size_t>(numberOfBlockElementsMSD, 0uz))),
        blockDataMSDSelf(maxNumberOfBlocksMSD,
                         std::vector<std::vector<double3>>(numberOfParticles,
                                                           std::vector<double3>(numberOfBlockElementsMSD, double3()))),
        msdSelf(maxNumberOfBlocksMSD,
                std::vector<std::vector<double4>>(numberOfComponents,
                                                  std::vector<double4>(numberOfBlockElementsMSD, double4()))),
        msdOnsagerCount(
            maxNumberOfBlocksMSD,
            std::vector<std::vector<size_t>>(numberOfComponents, std::vector<size_t>(numberOfBlockElementsMSD, 0uz))),
        blockDataMSDOnsager(maxNumberOfBlocksMSD,
                            std::vector<std::vector<double3>>(
                                numberOfComponents, std::vector<double3>(numberOfBlockElementsMSD, double3()))),
        msdOnsager(
            maxNumberOfBlocksMSD,
            std::vector<std::vector<std::vector<double4>>>(
                numberOfComponents, std::vector<std::vector<double4>>(
                                        numberOfComponents, std::vector<double4>(numberOfBlockElementsMSD, double4()))))
  {
  }

  uint64_t versionNumber{1};

  size_t numberOfBlocks;
  size_t sampleEvery;
  size_t writeEvery;
  size_t numberOfComponents;
  size_t numberOfParticles;
  size_t numberOfBlockElementsMSD;
  size_t maxNumberOfBlocksMSD;
  size_t countMSD{0uz};
  size_t numberOfBlocksMSD;
  std::vector<size_t> blockLengthMSD;

  std::vector<std::vector<std::vector<size_t>>> msdSelfCount;
  std::vector<std::vector<std::vector<double3>>> blockDataMSDSelf;
  std::vector<std::vector<std::vector<double4>>> msdSelf;

  std::vector<std::vector<std::vector<size_t>>> msdOnsagerCount;
  std::vector<std::vector<std::vector<double3>>> blockDataMSDOnsager;
  std::vector<std::vector<std::vector<std::vector<double4>>>> msdOnsager;

  void addSample(size_t currentCycle, const std::vector<Component> &components,
                 const std::vector<size_t> &numberOfMoleculesPerComponent, std::vector<Molecule> &molecules);
  void writeOutput(size_t systemId, const std::vector<Component> &components,
                   const std::vector<size_t> &numberOfMoleculesPerComponent, double deltaT, size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyMeanSquaredDisplacement &msd);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyMeanSquaredDisplacement &msd);
};
