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

export module property_widom;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import int3;
import averages;

export struct PropertyWidom
{
  struct BookKeeping
  {
    double RosenbluthValue;
    double tailCorrection;
    double count;

    friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BookKeeping &w)
    {
      archive << w.RosenbluthValue;
      archive << w.tailCorrection;
      archive << w.count;

      return archive;
    }

    friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BookKeeping &w)
    {
      archive >> w.RosenbluthValue;
      archive >> w.tailCorrection;
      archive >> w.count;

      return archive;
    }
  };

  inline static BookKeeping pair_acc_widom(const BookKeeping &lhs, const BookKeeping &rhs)
  {
    return {lhs.RosenbluthValue + rhs.RosenbluthValue, lhs.tailCorrection + rhs.tailCorrection, lhs.count + rhs.count};
  }

  PropertyWidom();

  PropertyWidom(std::size_t numberOfBlocks)
      : numberOfBlocks(numberOfBlocks), bookKeepingWidom(numberOfBlocks), bookKeepingDensity(numberOfBlocks)
  {
  }

  bool operator==(PropertyWidom const &) const = default;

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks;
  std::vector<BookKeeping> bookKeepingWidom;
  std::vector<std::pair<double, double>> bookKeepingDensity;

  std::string writeAveragesRosenbluthWeightStatistics(double temperature, double volume,
                                                      std::optional<double> frameworkMass,
                                                      std::optional<int3> number_of_unit_cells) const;
  std::string writeAveragesChemicalPotentialStatistics(double beta, std::optional<double> imposedChemicalPotential,
                                                       std::optional<double> imposedFugacity) const;

  inline void addWidomSample(std::size_t blockIndex, double WidomValue, double tailCorrection, double weight)
  {
    bookKeepingWidom[blockIndex].RosenbluthValue += weight * WidomValue;
    bookKeepingWidom[blockIndex].tailCorrection += weight * tailCorrection;
    bookKeepingWidom[blockIndex].count += weight;
  }

  inline void addDensitySample(std::size_t blockIndex, double density, double weight)
  {
    bookKeepingDensity[blockIndex].first += weight * density;
    bookKeepingDensity[blockIndex].second += weight;
  }

  //====================================================================================================================

  double averagedChemicalPotentialTailCorrection(std::size_t blockIndex) const
  {
    return bookKeepingWidom[blockIndex].tailCorrection / bookKeepingWidom[blockIndex].count;
  }

  double averagedChemicalPotentialTailCorrection() const
  {
    std::pair<double, double> summedBlocks{0.0, 0.0};
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      summedBlocks.first += bookKeepingWidom[blockIndex].tailCorrection;
      summedBlocks.second += bookKeepingWidom[blockIndex].count;
    }

    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<double, double> averageChemicalPotentialTailCorrection() const
  {
    double average = averagedChemicalPotentialTailCorrection();

    double sumOfSquares = 0.0;
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingWidom[blockIndex].count / bookKeepingWidom[0].count > 0.5)
      {
        double value = averagedChemicalPotentialTailCorrection(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = std::sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedRosenbluthWeight(std::size_t blockIndex) const
  {
    return bookKeepingWidom[blockIndex].RosenbluthValue / bookKeepingWidom[blockIndex].count;
  }

  double averagedRosenbluthWeight() const
  {
    BookKeeping summedBlocks =
        std::accumulate(bookKeepingWidom.begin(), bookKeepingWidom.end(), BookKeeping{}, pair_acc_widom);
    return summedBlocks.RosenbluthValue / summedBlocks.count;
  }

  std::pair<double, double> averageRosenbluthWeight() const
  {
    double average = averagedRosenbluthWeight();

    double sumOfSquares = 0.0;
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingWidom[blockIndex].count / bookKeepingWidom[0].count > 0.5)
      {
        double value = averagedRosenbluthWeight(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = std::sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedExcessChemicalPotential(std::size_t blockIndex, double beta) const
  {
    return -(1.0 / beta) *
           std::log((bookKeepingWidom[blockIndex].RosenbluthValue / bookKeepingWidom[blockIndex].count));
  }

  double averagedExcessChemicalPotential(double beta) const
  {
    BookKeeping summedBlocks =
        std::accumulate(bookKeepingWidom.begin(), bookKeepingWidom.end(), BookKeeping{}, pair_acc_widom);
    return -(1.0 / beta) * std::log((summedBlocks.RosenbluthValue / summedBlocks.count));
  }

  std::pair<double, double> averageExcessChemicalPotential(double beta) const
  {
    double average = averagedExcessChemicalPotential(beta);

    double sumOfSquares = 0.0;
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingWidom[blockIndex].count / bookKeepingWidom[0].count > 0.5)
      {
        double value = averagedExcessChemicalPotential(blockIndex, beta) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = std::sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedDensity(std::size_t blockIndex) const
  {
    return bookKeepingDensity[blockIndex].first / bookKeepingDensity[blockIndex].second;
  }

  double averagedDensity() const
  {
    std::pair<double, double> summedBlocks{0.0, 0.0};
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
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
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
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
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = std::sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedIdealGasChemicalPotential(std::size_t blockIndex, double beta) const
  {
    return std::log(bookKeepingDensity[blockIndex].first / bookKeepingDensity[blockIndex].second) / beta;
  }

  double averagedIdealGasChemicalPotential(double beta) const
  {
    std::pair<double, double> summedBlocks{0.0, 0.0};
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
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
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
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
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = std::sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  std::pair<double, double> averageTotalChemicalPotential(double beta) const
  {
    double average = averagedExcessChemicalPotential(beta) + averagedIdealGasChemicalPotential(beta) +
                     averagedChemicalPotentialTailCorrection();

    double sumOfSquares = 0.0;
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingWidom[blockIndex].count / bookKeepingWidom[0].count > 0.5)
      {
        double value =
            (averagedExcessChemicalPotential(blockIndex, beta) + averagedIdealGasChemicalPotential(blockIndex, beta) +
             averagedChemicalPotentialTailCorrection(blockIndex)) -
            average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = std::sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  std::pair<double, double> averageFugacity(double beta) const
  {
    double average = std::exp(beta * (averagedExcessChemicalPotential(beta) + averagedIdealGasChemicalPotential(beta) +
                                      averagedChemicalPotentialTailCorrection())) /
                     beta;

    double sumOfSquares = 0.0;
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingWidom[blockIndex].count / bookKeepingWidom[0].count > 0.5)
      {
        double value = std::exp(beta * (averagedExcessChemicalPotential(blockIndex, beta) +
                                        averagedIdealGasChemicalPotential(
                                            blockIndex, beta + averagedChemicalPotentialTailCorrection(blockIndex)))) /
                           beta -
                       average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = std::sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyWidom &w);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyWidom &w);
};
