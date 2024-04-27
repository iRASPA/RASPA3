module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <array>
#include <optional>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>
#include <numbers>
#include <tuple>
#include <iostream>
#include <fstream>
#endif

export module property_pressure;

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
import double3x3;
import averages;
import units;

inline std::pair<double, double> pair_acc(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

inline std::pair<double3x3, double> 
pair_acc2(const std::pair<double3x3, double> &lhs, const std::pair<double3x3, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}


export struct PropertyPressure
{
  PropertyPressure() {};

  PropertyPressure(size_t numberOfBlocks) :
      numberOfBlocks(numberOfBlocks),
      bookKeepingExcessPressure(numberOfBlocks),
      bookKeepingIdealGasPressure(numberOfBlocks)
  {
  }

  bool operator==(PropertyPressure const&) const = default;

  uint64_t versionNumber{ 1 };
  size_t numberOfBlocks;
  std::vector<std::pair<double3x3, double>> bookKeepingExcessPressure;
  std::vector<std::pair<double, double>> bookKeepingIdealGasPressure;

  inline void addSample(size_t blockIndex, double idealGasPressureValue, double3x3 excessPressureValue, double weight)
  {
    bookKeepingIdealGasPressure[blockIndex].first += weight * idealGasPressureValue;
    bookKeepingIdealGasPressure[blockIndex].second += weight;

    bookKeepingExcessPressure[blockIndex].first += weight * excessPressureValue;
    bookKeepingExcessPressure[blockIndex].second += weight;
  }

  //====================================================================================================================

  double3x3 averagedExcessPressureTensor(size_t blockIndex) const
  {
    return bookKeepingExcessPressure[blockIndex].first / bookKeepingExcessPressure[blockIndex].second;
  }

  double3x3 averagedExcessPressureTensor() const
  {
    std::pair<double3x3,double> summedBlocks = 
      std::accumulate(bookKeepingExcessPressure.begin(), bookKeepingExcessPressure.end(), 
                      std::make_pair(double3x3(0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0), 0.0), pair_acc2);
    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<double3x3, double3x3> averageExcessPressureTensor() const
  {
    double3x3 average = averagedExcessPressureTensor();

    double3x3 sumOfSquares{};
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingExcessPressure[blockIndex].second / bookKeepingExcessPressure[0].second > 0.5)
      {
        double3x3 value = averagedExcessPressureTensor(blockIndex) - average;
        sumOfSquares += sqr(value);
        ++numberOfSamples;
      }
    }

    double3x3 confidenceIntervalError{};
    if(numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      double3x3 standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double3x3 standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double3x3 averagedIdealGasPressureTensor(size_t blockIndex) const
  {
    return (bookKeepingIdealGasPressure[blockIndex].first / 
            bookKeepingIdealGasPressure[blockIndex].second) * double3x3::identity();
  }

  double3x3 averagedIdealGasPressureTensor() const
  {
    std::pair<double3x3,double> summedBlocks{};
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      summedBlocks.first += bookKeepingIdealGasPressure[blockIndex].first * double3x3::identity();
      summedBlocks.second += bookKeepingIdealGasPressure[blockIndex].second;
    }

    return (summedBlocks.first / summedBlocks.second);
  }

  std::pair<double3x3, double3x3> averageIdealGasPressureTensor() const
  {
    double3x3 average = averagedIdealGasPressureTensor();

    double3x3 sumOfSquares{};
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingIdealGasPressure[blockIndex].second / bookKeepingIdealGasPressure[0].second > 0.5)
      {
        double3x3 value = averagedIdealGasPressureTensor(blockIndex) - average;
        sumOfSquares += sqr(value);
        ++numberOfSamples;
      }
    }

    double3x3 confidenceIntervalError{};
    if(numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      double3x3 standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double3x3 standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }


  //====================================================================================================================

  double3x3 averagedPressureTensor(size_t blockIndex) const
  {
    return (bookKeepingIdealGasPressure[blockIndex].first / bookKeepingIdealGasPressure[blockIndex].second) * double3x3::identity() + 
           (bookKeepingExcessPressure[blockIndex].first / bookKeepingExcessPressure[blockIndex].second);
  }

  std::pair<double3x3, double3x3> averagePressureTensor() const
  {
    double3x3 average = averagedExcessPressureTensor()  + averagedIdealGasPressureTensor();

    double3x3 sumOfSquares{};
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingExcessPressure[blockIndex].second / bookKeepingExcessPressure[0].second > 0.5)
      {
        double3x3 value = 
          (averagedExcessPressureTensor(blockIndex) + averagedIdealGasPressureTensor(blockIndex)) - average;
        sumOfSquares += sqr(value);
        ++numberOfSamples;
      }
    }

    double3x3 confidenceIntervalError{};
    if(numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      double3x3 standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double3x3 standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedExcessPressure(size_t blockIndex) const
  {
    return bookKeepingExcessPressure[blockIndex].first.trace() / (3.0 * bookKeepingExcessPressure[blockIndex].second);
  }

  double averagedExcessPressure() const
  {
    std::pair<double3x3,double> summedBlocks = 
      std::accumulate(bookKeepingExcessPressure.begin(), bookKeepingExcessPressure.end(), 
                      std::make_pair(double3x3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0.0), pair_acc2);
    return summedBlocks.first.trace() / (3.0 * summedBlocks.second);
  }

  std::pair<double, double> averageExcessPressure() const
  {
    double average = averagedExcessPressure();

    double sumOfSquares{ 0.0 };
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingExcessPressure[blockIndex].second / bookKeepingExcessPressure[0].second > 0.5)
      {
        double value = averagedExcessPressure(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError{0.0};
    if(numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedIdealGasPressure(size_t blockIndex) const
  {
    return (bookKeepingIdealGasPressure[blockIndex].first / bookKeepingIdealGasPressure[blockIndex].second);
  }

  double averagedIdealGasPressure() const
  {
    std::pair<double,double> summedBlocks{0.0, 0.0};
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
        summedBlocks.first += bookKeepingIdealGasPressure[blockIndex].first;
        summedBlocks.second += bookKeepingIdealGasPressure[blockIndex].second;
    }

    return (summedBlocks.first / summedBlocks.second);
  }

  std::pair<double, double> averageIdealGasPressure() const
  {
    double average = averagedIdealGasPressure();

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingIdealGasPressure[blockIndex].second / bookKeepingIdealGasPressure[0].second > 0.5)
      {
        double value = averagedIdealGasPressure(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
    if(numberOfSamples >= 3)
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

  double averagedPressure(size_t blockIndex) const
  {
    return (bookKeepingIdealGasPressure[blockIndex].first / bookKeepingIdealGasPressure[blockIndex].second) + 
           (bookKeepingExcessPressure[blockIndex].first.trace() / (3.0 * bookKeepingExcessPressure[blockIndex].second));
  }

  std::pair<double, double> averagePressure() const
  {
    double average = averagedExcessPressure()  + averagedIdealGasPressure();

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingExcessPressure[blockIndex].second / bookKeepingExcessPressure[0].second > 0.5)
      {
        double value = (averagedExcessPressure(blockIndex) + averagedIdealGasPressure(blockIndex)) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError = 0.0;
    if(numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  std::string writeAveragesStatistics() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPressure &e);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPressure &e);
};
