module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

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

export module property_loading;

#ifdef USE_STD_IMPORT
import std;
#endif

import archive;
import int3;
import averages;
import loadings;
import component;

inline std::pair<Loadings, double> pair_sum(const std::pair<Loadings, double> &lhs,
                                            const std::pair<Loadings, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyLoading
{
  PropertyLoading() {};

  PropertyLoading(std::size_t numberOfBlocks, std::size_t numberOfComponents)
      : numberOfBlocks(numberOfBlocks),
        numberOfComponents(numberOfComponents),
        bookKeepingLoadings(numberOfBlocks)
  {
    // workaround for g++
    for(std::size_t i = 0; i < numberOfBlocks; ++i)
    {
      Loadings loadings = Loadings(numberOfComponents);
      bookKeepingLoadings[i] = {loadings, 0.0};
    }
  }

  bool operator==(PropertyLoading const &) const = default;

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks;
  std::size_t numberOfComponents;
  std::vector<std::pair<Loadings, double>> bookKeepingLoadings;

  void resize(std::size_t newNumberOfComponents)
  {
    numberOfComponents = newNumberOfComponents;

    // workaround for g++
    bookKeepingLoadings = std::vector<std::pair<Loadings, double>>(numberOfBlocks);
    for(std::size_t i = 0; i < numberOfBlocks; ++i)
    {
      Loadings loadings = Loadings(numberOfComponents);
      bookKeepingLoadings[i] = {loadings, 0.0};
    }
  }

  inline void addSample(std::size_t blockIndex, const Loadings &loading, const double &weight)
  {
    bookKeepingLoadings[blockIndex].first += weight * loading;
    bookKeepingLoadings[blockIndex].second += weight;
  }

  //====================================================================================================================

  Loadings averagedLoading(std::size_t blockIndex) const
  {
    return bookKeepingLoadings[blockIndex].first / std::max(1.0, bookKeepingLoadings[blockIndex].second);
  }

  Loadings averagedLoading() const
  {
    // workaround for g++
    Loadings loadings = Loadings(numberOfComponents);

    std::pair<Loadings, double> summedBlocks = std::make_pair(loadings, 0.0);

    for(const std::pair<Loadings, double> &bookKeepingLoading : bookKeepingLoadings)
    {
      summedBlocks.first += bookKeepingLoading.first;
      summedBlocks.second += bookKeepingLoading.second;
    }

    return summedBlocks.first / std::max(1.0, summedBlocks.second);;
  }

  std::pair<Loadings, Loadings> averageLoading() const
  {
    Loadings average = averagedLoading();

    Loadings sumOfSquares(numberOfComponents);
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingLoadings[blockIndex].second / std::max(1.0, bookKeepingLoadings[0].second) > 0.5)
      {
        Loadings value = averagedLoading(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }
    Loadings confidenceIntervalError(numberOfComponents);
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      Loadings standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      Loadings standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  std::pair<double, double> averageLoadingNumberOfMolecules(std::size_t comp) const;

  std::string writeAveragesStatistics(std::vector<Component> components, std::optional<double> frameworkMass,
                                      std::optional<int3> numberOfUnitCells) const;

  /**
   * \brief Returns a string representation of the PropertyLoading.
   *
   * Generates a simple string indicating a test representation of the PropertyLoading.
   *
   * \return A string representing the PropertyLoading.
   */
  std::string repr() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyLoading &l);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyLoading &l);
};
