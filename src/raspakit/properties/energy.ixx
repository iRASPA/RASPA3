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

export module property_energy;

#ifdef USE_STD_IMPORT
import std;
#endif

import archive;
import averages;
import energy_status;
import energy_status_inter;
import energy_status_intra;
import framework;
import component;
import json;

inline std::pair<EnergyStatus, double> pair_sum(const std::pair<EnergyStatus, double> &lhs,
                                                const std::pair<EnergyStatus, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyEnergy
{
  PropertyEnergy() {};

  PropertyEnergy(std::size_t numberOfBlocks, std::size_t numberOfExternalFields, std::size_t numberOfFrameworks,
                 std::size_t numberOfComponents)
      : numberOfBlocks(numberOfBlocks),
        numberOfComponents(numberOfComponents),
        bookKeepingEnergyStatus(std::vector<std::pair<EnergyStatus, double>>(
            numberOfBlocks,
            std::make_pair(EnergyStatus(numberOfExternalFields, numberOfFrameworks, numberOfComponents), 0.0)))
  {
  }

  std::uint64_t versionNumber{1};
  std::size_t numberOfBlocks{5};
  std::size_t numberOfExternalFields{1};
  std::size_t numberOfFrameworks{1};
  std::size_t numberOfComponents{1};
  std::vector<std::pair<EnergyStatus, double>> bookKeepingEnergyStatus;

  void resize(std::size_t newNumberOfFrameworks, std::size_t newNumberOfComponents)
  {
    numberOfComponents = newNumberOfComponents;
    numberOfFrameworks = newNumberOfFrameworks;
    bookKeepingEnergyStatus = std::vector<std::pair<EnergyStatus, double>>(
        numberOfBlocks,
        std::make_pair(EnergyStatus(numberOfExternalFields, numberOfFrameworks, numberOfComponents), 0.0));
  }

  inline void addSample(std::size_t blockIndex, const EnergyStatus &energyStatus, const double &weight)
  {
    bookKeepingEnergyStatus[blockIndex].first += weight * energyStatus;
    bookKeepingEnergyStatus[blockIndex].second += weight;
  }

  //====================================================================================================================

  EnergyStatus averagedEnergy(std::size_t blockIndex) const
  {
    return bookKeepingEnergyStatus[blockIndex].first / std::max(1.0, bookKeepingEnergyStatus[blockIndex].second);
  }

  EnergyStatus averagedEnergy() const
  {
    std::pair<EnergyStatus, double> summedBlocks = std::accumulate(
        bookKeepingEnergyStatus.begin(), bookKeepingEnergyStatus.end(),
        std::make_pair(EnergyStatus(numberOfExternalFields, numberOfFrameworks, numberOfComponents), 0.0), pair_sum);
    return summedBlocks.first / std::max(1.0, summedBlocks.second);
  }

  std::pair<EnergyStatus, EnergyStatus> averageEnergy() const
  {
    EnergyStatus average = averagedEnergy();

    // Use bins that are at least 90% filled for the computation of the error
    EnergyStatus sumOfSquares(numberOfExternalFields, numberOfFrameworks, numberOfComponents);
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingEnergyStatus[blockIndex].second / std::max(1.0, bookKeepingEnergyStatus[0].second) > 0.5)
      {
        EnergyStatus value = averagedEnergy(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }
    EnergyStatus confidenceIntervalError(numberOfExternalFields, numberOfFrameworks, numberOfComponents);
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      EnergyStatus standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      EnergyStatus standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }
    return std::make_pair(average, confidenceIntervalError);
  }

  std::vector<EnergyStatus> blockEnergy() const
  {
    std::vector<EnergyStatus> blockEnergies(numberOfBlocks);
    std::transform(bookKeepingEnergyStatus.begin(), bookKeepingEnergyStatus.end(), blockEnergies.begin(),
                   [](std::pair<EnergyStatus, double> block) { return block.first / block.second; });
    return blockEnergies;
  }

  std::string writeAveragesStatistics(bool externalField, std::optional<Framework> &framework,
                                      std::vector<Component> &components) const;
  nlohmann::json jsonAveragesStatistics(bool externalField, std::optional<Framework> &framework,
                                        std::vector<Component> &components) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEnergy &e);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEnergy &e);
};
