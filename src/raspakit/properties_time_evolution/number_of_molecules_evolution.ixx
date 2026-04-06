module;

export module property_number_of_molecules_evolution;

import std;

import double4;

import archive;
import averages;
import simulationbox;
import component;

export struct PropertyNumberOfMoleculesEvolution
{
  PropertyNumberOfMoleculesEvolution() {};

  PropertyNumberOfMoleculesEvolution(std::size_t numberOfCycles, std::size_t numberOfComponents, 
                                     std::size_t sampleEvery, std::optional<std::size_t> writeEvery)
      : sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        totalSamples(numberOfCycles / sampleEvery),
        result(numberOfComponents, std::vector<std::size_t>(totalSamples))
  {
  }

  std::uint64_t versionNumber{1};

  std::size_t sampleEvery;
  std::optional<std::size_t> writeEvery;
  std::size_t totalSamples;
  std::vector<std::vector<std::size_t>> result;

  void addSample(std::size_t absoluteCurrentCycle, const std::vector<std::size_t> &numberOfIntegerMoleculesPerComponent);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyNumberOfMoleculesEvolution &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyNumberOfMoleculesEvolution &hist);
};
