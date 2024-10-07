module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
#endif

export module property_temperature;

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

import archive;
import averages;
import simulationbox;

inline std::pair<double, double> pair_sum(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyTemperature
{
  PropertyTemperature() {};

  PropertyTemperature(size_t numberOfBlocks)
      : numberOfBlocks(numberOfBlocks), bookKeepingTemperature(numberOfBlocks, std::make_pair(0.0, 0.0))
  {
  }

  uint64_t versionNumber{1};

  size_t numberOfBlocks;
  std::vector<std::pair<double, double>> bookKeepingTemperature;

  inline void addSample(size_t blockIndex, const double &temperature, const double &weight)
  {
    bookKeepingTemperature[blockIndex].first += weight * temperature;
    bookKeepingTemperature[blockIndex].second += weight;
  }

  //====================================================================================================================

  double averagedTemperature(size_t blockIndex) const
  {
    return bookKeepingTemperature[blockIndex].first / bookKeepingTemperature[blockIndex].second;
  }

  double averagedTemperature() const
  {
    std::pair<double, double> summedBlocks = std::accumulate(
        bookKeepingTemperature.begin(), bookKeepingTemperature.end(), std::make_pair(0.0, 0.0), pair_sum);
    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<double, double> averageTemperature() const
  {
    double average = averagedTemperature();

    double sumOfSquares{};
    size_t numberOfSamples = 0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingTemperature[blockIndex].second / bookKeepingTemperature[0].second > 0.5)
      {
        double value = averagedTemperature(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError{};
    if (numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyTemperature &temp);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyTemperature &temp);
};
