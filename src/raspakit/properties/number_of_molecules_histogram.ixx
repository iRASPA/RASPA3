module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numbers>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
#endif

export module property_number_of_molecules_histogram;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double4;

import archive;
import averages;
import simulationbox;
import component;

inline std::pair<double, double> pair_sum(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyNumberOfMoleculesHistogram
{
  PropertyNumberOfMoleculesHistogram() {};

  PropertyNumberOfMoleculesHistogram(std::size_t numberOfBlocks, std::pair<std::size_t, std::size_t> range,
                                     std::size_t size, std::size_t sampleEvery, std::size_t writeEvery)
      : numberOfBlocks(numberOfBlocks),
        numberOfBins(range.second - range.first),
        range(range),
        size(size),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        bookKeepingEnergyHistogram(std::vector<std::vector<std::vector<double>>>(
            numberOfBlocks, std::vector<std::vector<double>>(numberOfBins, std::vector<double>(size)))),
        numberOfCounts(numberOfBlocks)
  {
  }

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks;
  std::size_t numberOfBins;
  std::pair<std::size_t, std::size_t> range;
  std::size_t size;
  std::size_t sampleEvery;
  std::size_t writeEvery;
  std::vector<std::vector<std::vector<double>>> bookKeepingEnergyHistogram;
  std::vector<double> numberOfCounts;
  double totalNumberOfCounts{0.0};

  void addSample(std::size_t blockIndex, std::size_t currentCycle,
                 std::vector<std::size_t> numberOfIntegerMoleculesPerComponent, const double &weight);

  std::vector<std::vector<double>> averagedProbabilityHistogram(std::size_t blockIndex) const;
  std::vector<std::vector<double>> averagedProbabilityHistogram() const;
  std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> averageProbabilityHistogram() const;

  void writeOutput(std::size_t systemId, std::vector<Component> &components, std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyNumberOfMoleculesHistogram &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyNumberOfMoleculesHistogram &hist);
};
