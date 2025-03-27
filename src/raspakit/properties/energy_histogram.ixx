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

export module property_energy_histogram;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <array>;
import <optional>;
import <cmath>;
import <string>;
import <algorithm>;
import <numeric>;
import <numbers>;
import <tuple>;
import <iostream>;
import <fstream>;
#endif

import double4;

import archive;
import averages;
import simulationbox;

inline std::pair<double, double> pair_sum(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyEnergyHistogram
{
  PropertyEnergyHistogram() {};

  PropertyEnergyHistogram(size_t numberOfBlocks, size_t numberOfBins, std::pair<double, double> range,
                          size_t sampleEvery, size_t writeEvery)
      : numberOfBlocks(numberOfBlocks),
        numberOfBins(numberOfBins),
        range(range),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        bookKeepingEnergyHistogram(
            std::vector<std::vector<double4>>(numberOfBlocks, std::vector<double4>(numberOfBins))),
        numberOfCounts(numberOfBlocks)
  {
  }

  uint64_t versionNumber{1};

  size_t numberOfBlocks;
  size_t numberOfBins;
  std::pair<double, double> range;
  size_t sampleEvery;
  size_t writeEvery;
  std::vector<std::vector<double4>> bookKeepingEnergyHistogram;
  std::vector<double> numberOfCounts;
  double totalNumberOfCounts{0.0};

  void addSample(size_t blockIndex, size_t currentCycle, double4 energy, const double &weight);

  std::vector<double4> averagedProbabilityHistogram(size_t blockIndex) const;
  std::vector<double4> averagedProbabilityHistogram() const;
  std::pair<std::vector<double4>, std::vector<double4>> averageProbabilityHistogram() const;

  void writeOutput(size_t systemId, size_t currentCycle);

  std::string printSettings() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEnergyHistogram &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEnergyHistogram &hist);
};
