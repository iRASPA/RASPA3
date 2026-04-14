module;

export module property_volume_evolution;

import std;

import double4;

import archive;
import averages;
import simulationbox;
import component;

export struct PropertyVolumeEvolution
{
  PropertyVolumeEvolution() {};

  PropertyVolumeEvolution(std::size_t numberOfCycles, std::size_t sampleEvery, std::optional<std::size_t> writeEvery)
      : sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        totalSamples(numberOfCycles / sampleEvery),
        result(totalSamples)
  {
  }

  std::uint64_t versionNumber{1};

  std::size_t sampleEvery;
  std::optional<std::size_t> writeEvery;
  std::size_t totalSamples;
  std::vector<double> result;

  void addSample(std::size_t absoluteCurrentCycle, double currentVolume);

  void writeOutput(std::size_t systemId, std::size_t absoluteCurrentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyVolumeEvolution &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyVolumeEvolution &hist);
};
