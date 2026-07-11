module;

export module property_energy_histogram;

import std;

import archive;
import averages;
import simulationbox;
import average_energy_type;
import units;
export import property_block_average;

export struct PropertyEnergyHistogram
{
  PropertyEnergyHistogram() {};

  PropertyEnergyHistogram(std::size_t numberOfBlocks, std::size_t numberOfBins, std::pair<double, double> valueRange,
                          std::size_t sampleEvery, std::optional<std::size_t> writeEvery)
      : numberOfBins(numberOfBins),
        valueRange(Units::KelvinToEnergy * valueRange.first, Units::KelvinToEnergy * valueRange.second),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        histogram(numberOfBlocks, 1, numberOfBins)
  {
  }

  std::uint64_t versionNumber{2};

  std::size_t numberOfBins;
  std::pair<double, double> valueRange;
  std::size_t sampleEvery;
  std::optional<std::size_t> writeEvery;
  BlockHistogram<AverageEnergyType> histogram;

  void addSample(std::size_t blockIndex, std::size_t currentCycle, AverageEnergyType energy, const double &weight);

  std::vector<AverageEnergyType> averagedProbabilityHistogram(std::size_t blockIndex) const
  {
    return histogram.averaged(blockIndex, 0);
  }
  std::vector<AverageEnergyType> averagedProbabilityHistogram() const { return histogram.averaged(0); }

  std::tuple<std::vector<double>, std::vector<AverageEnergyType>, std::vector<AverageEnergyType>> result() const;
  std::vector<double> bins() const;

  void writeOutput(std::size_t systemId, std::size_t currentCycle);

  std::string printSettings() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEnergyHistogram &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEnergyHistogram &hist);
};
