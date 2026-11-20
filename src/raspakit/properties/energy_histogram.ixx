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

  /// The range [first, second) in Kelvin: a reversed or empty range is refused (the bins were laid out from 'first'
  /// over |second - first|, so a reversed range binned the interval above 'first' instead of the one below it).
  static std::pair<double, double> checkedValueRange(std::pair<double, double> range)
  {
    if (!(range.second > range.first))
      throw std::invalid_argument(std::format(
          "PropertyEnergyHistogram: the upper limit ({}) must lie above the lower limit ({})", range.second,
          range.first));
    return range;
  }

  static std::size_t checkedNumberOfBins(std::size_t numberOfBins)
  {
    if (numberOfBins == 0uz) throw std::invalid_argument("PropertyEnergyHistogram: the number of bins must be positive");
    return numberOfBins;
  }

  PropertyEnergyHistogram(std::size_t numberOfBlocks, std::size_t numberOfBins, std::pair<double, double> valueRange,
                          std::size_t sampleEvery, std::optional<std::size_t> writeEvery)
      : numberOfBins(checkedNumberOfBins(numberOfBins)),
        valueRange(Units::KelvinToEnergy * checkedValueRange(valueRange).first,
                   Units::KelvinToEnergy * valueRange.second),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        histogram(numberOfBlocks, 1, numberOfBins)
  {
  }

  std::uint64_t versionNumber{1};

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

  /// The lower edges of the bins in the unit the range was given in [K], with the averaged histogram and its error.
  std::tuple<std::vector<double>, std::vector<AverageEnergyType>, std::vector<AverageEnergyType>> result() const;
  /// The lower edges of the bins in the unit the range was given in [K].
  std::vector<double> bins() const;

  void writeOutput(std::size_t systemId, std::size_t currentCycle);

  std::string printSettings() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEnergyHistogram &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEnergyHistogram &hist);
};
