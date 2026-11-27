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
    // an interval of zero is a division by zero at the first sample (and in the number of samples)
    if (sampleEvery == 0uz) throw std::runtime_error("PropertyConservedEnergyEvolution: the sampling interval must be at least 1\n");
    if (writeEvery.has_value() && writeEvery.value() == 0uz)
    {
      throw std::runtime_error("PropertyConservedEnergyEvolution: the writing interval must be at least 1\n");
    }
  }

  PropertyConservedEnergyEvolution(std::size_t numberOfProductionCycles, std::size_t sampleEvery,
                                   std::optional<std::size_t> writeEvery)
      : sampleEvery(sampleEvery),
        writeEvery(writeEvery)
  {
    // an interval of zero is a division by zero at the first sample (and in the number of samples)
    if (sampleEvery == 0uz) throw std::runtime_error("PropertyConservedEnergyEvolution: the sampling interval must be at least 1\n");
    if (writeEvery.has_value() && writeEvery.value() == 0uz)
    {
      throw std::runtime_error("PropertyConservedEnergyEvolution: the writing interval must be at least 1\n");
    }
    totalSamples = (numberOfProductionCycles + sampleEvery - 1) / sampleEvery;
    data.resize(totalSamples);
  }

  std::uint64_t versionNumber{1};

  std::size_t sampleEvery;
  std::optional<std::size_t> writeEvery;
  std::size_t totalSamples{0};  // (it was left undefined by the constructor without a number of cycles)
  std::vector<RunningEnergy> data;

  const std::vector<RunningEnergy> result() const {return data;}

  void addSample(std::size_t absoluteCurrentCycle, const RunningEnergy &data);

  /// Writes the samples collected up to 'absoluteCurrentCycle' when it is a multiple of 'writeEvery'; the last write
  /// ('lastWrite', with the number of cycles completed) skips that check and writes every collected sample.
  void writeOutput(std::size_t systemId, std::size_t absoluteCurrentCycle, bool lastWrite = false);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyConservedEnergyEvolution &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyConservedEnergyEvolution &hist);
};
