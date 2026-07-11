module;

export module property_number_of_molecules_histogram;

import std;

import double4;

import archive;
import averages;
import simulationbox;
import component;
export import property_block_average;

export struct PropertyNumberOfMoleculesHistogram
{
  PropertyNumberOfMoleculesHistogram() {};


    PropertyNumberOfMoleculesHistogram(std::size_t numberOfBlocks, std::pair<std::size_t, std::size_t> range,
                                       std::size_t sampleEvery, std::optional<std::size_t> writeEvery):
        numberOfBlocks(numberOfBlocks),
        numberOfBins(range.second - range.first),
        range(range),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery)
  {
  }

  PropertyNumberOfMoleculesHistogram(std::size_t numberOfBlocks, std::size_t numberOfComponents,
                                     std::pair<std::size_t, std::size_t> range,
                                     std::size_t sampleEvery, std::optional<std::size_t> writeEvery)
      : numberOfBlocks(numberOfBlocks),
        numberOfBins(range.second - range.first),
        range(range),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        histogram(numberOfBlocks, numberOfComponents, range.second - range.first)
  {
  }

  std::uint64_t versionNumber{2};

  std::size_t numberOfBlocks;
  std::size_t numberOfBins;
  std::pair<std::size_t, std::size_t> range;
  std::size_t sampleEvery;
  std::optional<std::size_t> writeEvery;
  BlockHistogram<double> histogram;

  void addSample(std::size_t blockIndex, std::size_t currentCycle,
                 std::vector<std::size_t> numberOfIntegerMoleculesPerComponent, const double &weight);

  std::vector<double> averagedProbabilityHistogram(std::size_t blockIndex, std::size_t component_id) const
  {
    return histogram.averaged(blockIndex, component_id);
  }
  std::vector<double> averagedProbabilityHistogram(std::size_t component_id) const
  {
    return histogram.averaged(component_id);
  }

  std::vector<std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>> result() const;

  void writeOutput(std::size_t systemId, std::vector<Component> &components, std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyNumberOfMoleculesHistogram &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyNumberOfMoleculesHistogram &hist);
};
