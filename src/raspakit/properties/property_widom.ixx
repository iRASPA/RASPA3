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

export module property_widom;

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

inline std::pair<double, double> pair_acc(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyWidom
{
  PropertyWidom();

  PropertyWidom(size_t numberOfBlocks)
      : numberOfBlocks(numberOfBlocks), bookKeepingWidom(numberOfBlocks), bookKeepingDensity(numberOfBlocks)
  {
  }

  bool operator==(PropertyWidom const &) const = default;

  uint64_t versionNumber{1};

  size_t numberOfBlocks;
  std::vector<std::pair<double, double>> bookKeepingWidom;
  std::vector<std::pair<double, double>> bookKeepingDensity;

  std::string writeAveragesStatistics(double beta, std::optional<double> imposedChemicalPotential,
                                      std::optional<double> imposedFugacity) const;

  inline void addWidomSample(size_t blockIndex, double WidomValue, double weight)
  {
    bookKeepingWidom[blockIndex].first += weight * WidomValue;
    bookKeepingWidom[blockIndex].second += weight;
  }

  inline void addDensitySample(size_t blockIndex, double density, double weight)
  {
    bookKeepingDensity[blockIndex].first += weight * density;
    bookKeepingDensity[blockIndex].second += weight;
  }

  //====================================================================================================================

  double averagedRosenbluthWeight(size_t blockIndex) const
  {
    return bookKeepingWidom[blockIndex].first / bookKeepingWidom[blockIndex].second;
  }

  double averagedRosenbluthWeight() const
  {
    std::pair<double, double> summedBlocks =
        std::accumulate(bookKeepingWidom.begin(), bookKeepingWidom.end(), std::make_pair(0.0, 0.0), pair_acc);
    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<double, double> averageRosenbluthWeight() const
  {
    double average = averagedRosenbluthWeight();

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingWidom[blockIndex].second / bookKeepingWidom[0].second > 0.5)
      {
        double value = averagedRosenbluthWeight(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
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

  //====================================================================================================================

  double averagedExcessChemicalPotential(size_t blockIndex, double beta) const
  {
    return -(1.0 / beta) * std::log((bookKeepingWidom[blockIndex].first / bookKeepingWidom[blockIndex].second));
  }

  double averagedExcessChemicalPotential(double beta) const
  {
    std::pair<double, double> summedBlocks =
        std::accumulate(bookKeepingWidom.begin(), bookKeepingWidom.end(), std::make_pair(0.0, 0.0), pair_acc);
    return -(1.0 / beta) * std::log((summedBlocks.first / summedBlocks.second));
  }

  std::pair<double, double> averageExcessChemicalPotential(double beta) const
  {
    double average = averagedExcessChemicalPotential(beta);

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingWidom[blockIndex].second / bookKeepingWidom[0].second > 0.5)
      {
        double value = averagedExcessChemicalPotential(blockIndex, beta) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
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

  //====================================================================================================================

  double averagedDensity(size_t blockIndex) const
  {
    return bookKeepingDensity[blockIndex].first / bookKeepingDensity[blockIndex].second;
  }

  double averagedDensity() const
  {
    std::pair<double, double> summedBlocks{0.0, 0.0};
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      summedBlocks.first += bookKeepingDensity[blockIndex].first;
      summedBlocks.second += bookKeepingDensity[blockIndex].second;
    }

    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<double, double> averageDensity() const
  {
    double average = averagedDensity();

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingDensity[blockIndex].second / bookKeepingDensity[0].second > 0.5)
      {
        double value = averagedDensity(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
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

  //====================================================================================================================

  double averagedIdealGasChemicalPotential(size_t blockIndex, double beta) const
  {
    return std::log(bookKeepingDensity[blockIndex].first / bookKeepingDensity[blockIndex].second) / beta;
  }

  double averagedIdealGasChemicalPotential(double beta) const
  {
    std::pair<double, double> summedBlocks{0.0, 0.0};
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      summedBlocks.first += bookKeepingDensity[blockIndex].first;
      summedBlocks.second += bookKeepingDensity[blockIndex].second;
    }

    return std::log(summedBlocks.first / summedBlocks.second) / beta;
  }

  std::pair<double, double> averageIdealGasChemicalPotential(double beta) const
  {
    double average = averagedIdealGasChemicalPotential(beta);

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingDensity[blockIndex].second / bookKeepingDensity[0].second > 0.5)
      {
        double value = averagedIdealGasChemicalPotential(blockIndex, beta) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
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

  //====================================================================================================================

  std::pair<double, double> averageTotalChemicalPotential(double beta) const
  {
    double average = averagedExcessChemicalPotential(beta) + averagedIdealGasChemicalPotential(beta);

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingWidom[blockIndex].second / bookKeepingWidom[0].second > 0.5)
      {
        double value =
            (averagedExcessChemicalPotential(blockIndex, beta) + averagedIdealGasChemicalPotential(blockIndex, beta)) -
            average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
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

  std::pair<double, double> averageFugacity(double beta) const
  {
    double average =
        std::exp(beta * (averagedExcessChemicalPotential(beta) + averagedIdealGasChemicalPotential(beta))) / beta;

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingWidom[blockIndex].second / bookKeepingWidom[0].second > 0.5)
      {
        double value = std::exp(beta * (averagedExcessChemicalPotential(blockIndex, beta) +
                                        averagedIdealGasChemicalPotential(blockIndex, beta))) /
                           beta -
                       average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
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

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyWidom &w);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyWidom &w);
};
