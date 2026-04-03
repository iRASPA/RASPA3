module;

export module property_energy_histogram;

import std;

import archive;
import averages;
import simulationbox;
import average_energy_type;

inline std::pair<double, double> pair_sum(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}


export struct PropertyEnergyHistogram
{

  PropertyEnergyHistogram() {};

  PropertyEnergyHistogram(std::size_t numberOfBlocks, std::size_t numberOfBins, std::pair<double, double> range,
                          std::size_t sampleEvery, std::size_t writeEvery)
      : numberOfBlocks(numberOfBlocks),
        numberOfBins(numberOfBins),
        range(range),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        bookKeepingEnergyHistogram(
            std::vector<std::vector<AverageEnergyType>>(numberOfBlocks, std::vector<AverageEnergyType>(numberOfBins))),
        numberOfCounts(numberOfBlocks)
  {
  }

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks;
  std::size_t numberOfBins;
  std::pair<double, double> range;
  std::size_t sampleEvery;
  std::size_t writeEvery;
  std::vector<std::vector<AverageEnergyType>> bookKeepingEnergyHistogram;
  std::vector<double> numberOfCounts;
  double totalNumberOfCounts{0.0};

  void addSample(std::size_t blockIndex, std::size_t currentCycle, AverageEnergyType energy, const double &weight);

  std::vector<AverageEnergyType> averagedProbabilityHistogram(std::size_t blockIndex) const;
  std::vector<AverageEnergyType> averagedProbabilityHistogram() const;

  std::tuple<std::vector<double>, std::vector<AverageEnergyType>, std::vector<AverageEnergyType>> result() const;
  std::vector<double> bins() const;

  void writeOutput(std::size_t systemId, std::size_t currentCycle);

  std::string printSettings() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEnergyHistogram &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEnergyHistogram &hist);
};
