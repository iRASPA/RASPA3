export module property_widom;

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

import averages;

inline std::pair<double, double> pair_acc(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyWidom
{
  PropertyWidom(size_t numberOfBlocks) :
      numberOfBlocks(numberOfBlocks),
      bookKeepingWidom(numberOfBlocks),
      bookKeepingDensity(numberOfBlocks)
  {
  }

  size_t numberOfBlocks;
  std::vector<std::pair<double, double>> bookKeepingWidom;
  std::vector<std::pair<double, double>> bookKeepingDensity;

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
    std::pair<double,double> summedBlocks = std::accumulate(bookKeepingWidom.begin(), bookKeepingWidom.end(), std::make_pair(0.0, 0.0), pair_acc);
    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<double, double> averageRosenbluthWeight() const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedRosenbluthWeight();

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = averagedRosenbluthWeight(blockIndex) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedExcessChemicalPotential(size_t blockIndex, double Beta) const
  {
    return -(1.0/Beta) * std::log((bookKeepingWidom[blockIndex].first / bookKeepingWidom[blockIndex].second));
  }

  double averagedExcessChemicalPotential(double Beta) const
  {
    std::pair<double,double> summedBlocks = std::accumulate (bookKeepingWidom.begin(), bookKeepingWidom.end(), std::make_pair(0.0, 0.0), pair_acc);
    return -(1.0/Beta) * std::log((summedBlocks.first / summedBlocks.second));
  }

  std::pair<double, double> averageExcessChemicalPotential(double Beta) const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedExcessChemicalPotential(Beta);

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = averagedExcessChemicalPotential(blockIndex, Beta) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedDensity(size_t blockIndex) const
  {
    return bookKeepingDensity[blockIndex].first / bookKeepingDensity[blockIndex].second;
  }

  double averagedDensity() const
  {
    std::pair<double,double> summedBlocks{0.0, 0.0};
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
        summedBlocks.first += bookKeepingDensity[blockIndex].first;
        summedBlocks.second += bookKeepingDensity[blockIndex].second;
    }

    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<double, double> averageDensity() const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedDensity();

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = averagedDensity(blockIndex) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedIdealGasChemicalPotential(size_t blockIndex, double Beta) const
  {
    return std::log(bookKeepingDensity[blockIndex].first / bookKeepingDensity[blockIndex].second) / Beta;
  }

  double averagedIdealGasChemicalPotential(double Beta) const
  {
    std::pair<double,double> summedBlocks{0.0, 0.0};
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
        summedBlocks.first += bookKeepingDensity[blockIndex].first;
        summedBlocks.second += bookKeepingDensity[blockIndex].second;
    }

    return std::log(summedBlocks.first / summedBlocks.second) / Beta;
  }

  std::pair<double, double> averageIdealGasChemicalPotential(double Beta) const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedIdealGasChemicalPotential(Beta);

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = averagedIdealGasChemicalPotential(blockIndex, Beta) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }


  //====================================================================================================================

  std::pair<double, double> averageTotalChemicalPotential(double Beta) const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedExcessChemicalPotential(Beta)  + averagedIdealGasChemicalPotential(Beta);

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = (averagedExcessChemicalPotential(blockIndex, Beta) + averagedIdealGasChemicalPotential(blockIndex, Beta)) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }

  std::pair<double, double> averageFugacity(double Beta) const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = std::exp(Beta*(averagedExcessChemicalPotential(Beta)  + averagedIdealGasChemicalPotential(Beta))) / Beta;

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = std::exp(Beta * (averagedExcessChemicalPotential(blockIndex, Beta) + averagedIdealGasChemicalPotential(blockIndex, Beta))) / Beta - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }


};
