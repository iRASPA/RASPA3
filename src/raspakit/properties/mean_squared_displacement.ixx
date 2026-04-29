module;

export module property_msd;

import std;

import archive;
import double3;
import double4;
import atom;
import molecule;
import molecule;
import simulationbox;
import forcefield;
import component;
import mean_squared_displacement_data;

// Computes Mean Squared Displacement (MSD) using order-N algorithm

export struct PropertyMeanSquaredDisplacement
{
  PropertyMeanSquaredDisplacement() {};

  PropertyMeanSquaredDisplacement(std::size_t sampleEvery, std::optional<std::size_t> writeEvery):
        sampleEvery(sampleEvery),
        writeEvery(writeEvery)
  {
  }

  PropertyMeanSquaredDisplacement(const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::size_t numberOfParticles,
                                  double timeStep, std::size_t numberOfBlockElementsMSD, 
                                  std::size_t sampleEvery, std::optional<std::size_t> writeEvery):
        numberOfMoleculesPerComponent(numberOfMoleculesPerComponent),
        numberOfComponents(numberOfMoleculesPerComponent.size()),
        numberOfParticles(numberOfParticles),
        timeStep(timeStep),
        numberOfBlockElementsMSD(numberOfBlockElementsMSD),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
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
  std::vector<std::size_t> numberOfMoleculesPerComponent;
  std::size_t numberOfComponents;
  std::size_t numberOfParticles;
  double timeStep;
  std::size_t numberOfBlockElementsMSD;
  std::size_t sampleEvery;
  std::optional<std::size_t> writeEvery;
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

  void addSample(std::size_t currentCycle, std::vector<Molecule> &molecules);

  std::vector<std::vector<MeanSquaredDisplacementData>> result();

  void writeOutput(std::size_t systemId, const std::vector<Component> &components,
                   std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyMeanSquaredDisplacement &msd);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyMeanSquaredDisplacement &msd);
};
