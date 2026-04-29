module;

export module property_conserved_energy_evolution;

import std;

import double4;

import archive;
import running_energy;

export struct PropertyConservedEnergyEvolution
{
  PropertyConservedEnergyEvolution() {};

  PropertyConservedEnergyEvolution(std::size_t sampleEvery, std::optional<std::size_t> writeEvery):
        sampleEvery(sampleEvery),
        writeEvery(writeEvery)
  {
  }


  PropertyConservedEnergyEvolution(std::size_t numberOfCycles, std::size_t sampleEvery, std::optional<std::size_t> writeEvery)
      : sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        totalSamples(numberOfCycles / sampleEvery),
        data(totalSamples)
  {
  }

  std::uint64_t versionNumber{1};

  std::size_t sampleEvery;
  std::optional<std::size_t> writeEvery;
  std::size_t totalSamples;
  std::vector<RunningEnergy> data;

  const std::vector<RunningEnergy> result() const {return data;}

  void addSample(std::size_t absoluteCurrentCycle, const RunningEnergy &data);

  void writeOutput(std::size_t systemId, std::size_t absoluteCurrentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyConservedEnergyEvolution &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyConservedEnergyEvolution &hist);
};
