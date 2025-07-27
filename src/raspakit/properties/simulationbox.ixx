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

export module property_simulationbox;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import averages;
import simulationbox;

inline std::pair<SimulationBox, double> pair_sum(const std::pair<SimulationBox, double> &lhs,
                                                 const std::pair<SimulationBox, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertySimulationBox
{
  PropertySimulationBox() {};

  PropertySimulationBox(std::size_t numberOfBlocks)
      : numberOfBlocks(numberOfBlocks), bookKeepingSimulationBox(numberOfBlocks, std::make_pair(SimulationBox(), 0.0))
  {
  }

  bool operator==(PropertySimulationBox const &) const = default;

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks;
  std::vector<std::pair<SimulationBox, double>> bookKeepingSimulationBox;

  inline void addSample(std::size_t blockIndex, const SimulationBox &box, const double &weight)
  {
    bookKeepingSimulationBox[blockIndex].first += weight * box;
    bookKeepingSimulationBox[blockIndex].second += weight;
  }

  //====================================================================================================================

  SimulationBox averagedSimulationBox(std::size_t blockIndex) const
  {
    return bookKeepingSimulationBox[blockIndex].first / bookKeepingSimulationBox[blockIndex].second;
  }

  SimulationBox averagedSimulationBox() const
  {
    std::pair<SimulationBox, double> summedBlocks =
        std::accumulate(bookKeepingSimulationBox.begin(), bookKeepingSimulationBox.end(),
                        std::make_pair(SimulationBox(), 0.0), pair_sum);
    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<SimulationBox, SimulationBox> averageSimulationBox() const
  {
    SimulationBox average = averagedSimulationBox();

    SimulationBox sumOfSquares{};
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingSimulationBox[blockIndex].second / bookKeepingSimulationBox[0].second > 0.5)
      {
        SimulationBox value = averagedSimulationBox(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    SimulationBox confidenceIntervalError{};
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      SimulationBox standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      SimulationBox standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertySimulationBox &box);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertySimulationBox &box);
};
