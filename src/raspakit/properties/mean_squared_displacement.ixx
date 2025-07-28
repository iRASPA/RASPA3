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
import std;
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

  PropertyMeanSquaredDisplacement(std::size_t numberOfComponents, std::size_t numberOfParticles,
                                  std::size_t sampleEvery, std::size_t writeEvery, std::size_t numberOfBlockElementsMSD)
      : sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        numberOfComponents(numberOfComponents),
        numberOfParticles(numberOfParticles),
        numberOfBlockElementsMSD(numberOfBlockElementsMSD),
        maxNumberOfBlocksMSD(1),
        blockLengthMSD(maxNumberOfBlocksMSD),
        msdSelfCount(maxNumberOfBlocksMSD,
                     std::vector<std::vector<std::size_t>>(numberOfComponents,
                                                           std::vector<std::size_t>(numberOfBlockElementsMSD, 0uz))),
        blockDataMSDSelf(maxNumberOfBlocksMSD,
                         std::vector<std::vector<double3>>(numberOfParticles,
                                                           std::vector<double3>(numberOfBlockElementsMSD, double3()))),
        msdSelf(maxNumberOfBlocksMSD,
                std::vector<std::vector<double4>>(numberOfComponents,
                                                  std::vector<double4>(numberOfBlockElementsMSD, double4()))),
        msdOnsagerCount(maxNumberOfBlocksMSD,
                        std::vector<std::vector<std::size_t>>(numberOfComponents,
                                                              std::vector<std::size_t>(numberOfBlockElementsMSD, 0uz))),
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

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks;
  std::size_t sampleEvery;
  std::size_t writeEvery;
  std::size_t numberOfComponents;
  std::size_t numberOfParticles;
  std::size_t numberOfBlockElementsMSD;
  std::size_t maxNumberOfBlocksMSD;
  std::size_t countMSD{0uz};
  std::size_t numberOfBlocksMSD;
  std::vector<std::size_t> blockLengthMSD;

  std::vector<std::vector<std::vector<std::size_t>>> msdSelfCount;
  std::vector<std::vector<std::vector<double3>>> blockDataMSDSelf;
  std::vector<std::vector<std::vector<double4>>> msdSelf;

  std::vector<std::vector<std::vector<std::size_t>>> msdOnsagerCount;
  std::vector<std::vector<std::vector<double3>>> blockDataMSDOnsager;
  std::vector<std::vector<std::vector<std::vector<double4>>>> msdOnsager;

  void addSample(std::size_t currentCycle, const std::vector<Component> &components,
                 const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::vector<Molecule> &molecules);
  void writeOutput(std::size_t systemId, const std::vector<Component> &components,
                   const std::vector<std::size_t> &numberOfMoleculesPerComponent, double deltaT,
                   std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyMeanSquaredDisplacement &msd);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyMeanSquaredDisplacement &msd);
};
