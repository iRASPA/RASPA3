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
import component;

inline std::pair<double, double> pair_sum(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyNumberOfMoleculesHistogram
{
  PropertyNumberOfMoleculesHistogram() {};

  PropertyNumberOfMoleculesHistogram(size_t numberOfBlocks, std::pair<size_t, size_t> range, size_t size,
                                     size_t sampleEvery, size_t writeEvery)
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

  uint64_t versionNumber{1};

  size_t numberOfBlocks;
  size_t numberOfBins;
  std::pair<size_t, size_t> range;
  size_t size;
  size_t sampleEvery;
  size_t writeEvery;
  std::vector<std::vector<std::vector<double>>> bookKeepingEnergyHistogram;
  std::vector<double> numberOfCounts;
  double totalNumberOfCounts{0.0};

  void addSample(size_t blockIndex, size_t currentCycle, std::vector<size_t> numberOfIntegerMoleculesPerComponent,
                 const double &weight);

  std::vector<std::vector<double>> averagedProbabilityHistogram(size_t blockIndex) const;
  std::vector<std::vector<double>> averagedProbabilityHistogram() const;
  std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> averageProbabilityHistogram() const;

  void writeOutput(size_t systemId, std::vector<Component> &components, size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyNumberOfMoleculesHistogram &hist);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyNumberOfMoleculesHistogram &hist);
};
