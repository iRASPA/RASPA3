export module property_dudlambda;

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

inline std::pair<double, double> pair_sum(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyDUDlambda
{
  PropertyDUDlambda(size_t numberOfBlocks, size_t numberOfBins) :
      numberOfBlocks(numberOfBlocks),
      numberOfBins(numberOfBins),
      bookKeepingDUdlambda(std::vector<std::vector<std::pair<double, double>>>(numberOfBlocks,
                           std::vector<std::pair<double, double>>(numberOfBins))),
      bookKeepingDensity(numberOfBlocks)
  {
  }

  size_t numberOfBlocks;
  size_t numberOfBins;
  std::vector<std::vector<std::pair<double, double>>> bookKeepingDUdlambda;
  std::vector<std::pair<double, double>> bookKeepingDensity;

  inline void addSample(size_t blockIndex, size_t binIndex, double dUdlambda)
  {
    bookKeepingDUdlambda[blockIndex][binIndex].first += binIndex < numberOfBins/2 ? 2.0 * dUdlambda : 0.0;
    bookKeepingDUdlambda[blockIndex][binIndex].second += 1.0;
  }

   inline void addDensitySample(size_t blockIndex, double density, double weight)
  {
    bookKeepingDensity[blockIndex].first += weight * density;
    bookKeepingDensity[blockIndex].second += weight;
  }

  //====================================================================================================================

  std::vector<double> averagedDUdlambda(size_t blockIndex) const
  {
    std::vector<double> averagedData(numberOfBins);
    std::transform(bookKeepingDUdlambda[blockIndex].cbegin(), bookKeepingDUdlambda[blockIndex].cend(), averagedData.begin(),
                   [&](const std::pair<double, double> &sample){return sample.first / std::max(1.0, sample.second);});
    return averagedData;
  }

  std::vector<double> averagedDUdlambda() const
  {
    std::vector<std::pair<double,double>> summedBlocks(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform (summedBlocks.begin(), summedBlocks.end(), bookKeepingDUdlambda[blockIndex].begin(), summedBlocks.begin(), pair_sum);
    }

    std::vector<double> averagedData(numberOfBins);
    std::transform(summedBlocks.begin(), summedBlocks.end(), averagedData.begin(),
                   [&](std::pair<double, double> &sample){return sample.first / std::max(1.0, sample.second);});
    return averagedData;
  }

  std::pair<std::vector<double>, std::vector<double>> averageDuDlambda() const
  {
    size_t degreesOfFreedom = numberOfBlocks - 1;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    std::vector<double> average = averagedDUdlambda();

    std::vector<double> sumOfSquares(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::vector<double> blockAverage = averagedDUdlambda(blockIndex);
      for(size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
      {
        double value = blockAverage[binIndex] - average[binIndex];
        sumOfSquares[binIndex] += value * value;
      }
    }
    std::vector<double> standardDeviation(numberOfBins);
    std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(),
                   [&](const double &sumofsquares){return std::sqrt(sumofsquares / static_cast<double>(degreesOfFreedom));});

    std::vector<double> standardError(numberOfBins);
    std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
                   [&](const double &sigma){return sigma / sqrt(static_cast<double>(numberOfBlocks));});

    std::vector<double> confidenceIntervalError(numberOfBins);
    std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
                   [&](const double &error){return intermediateStandardNormalDeviate * error;});

    return std::make_pair(average, confidenceIntervalError);
  }


  //====================================================================================================================

  double averagedExcessChemicalPotential(size_t blockIndex) const
  {
    std::vector<double> averagedData(numberOfBins);
    std::transform(bookKeepingDUdlambda[blockIndex].begin(), bookKeepingDUdlambda[blockIndex].end(), averagedData.begin(),
                   [&](const std::pair<double, double> &sample){return sample.first / std::max(1.0, sample.second);});

    //return std::accumulate(averagedData.begin(), averagedData.end(), 0.0) / static_cast<double>(numberOfBins - 1);
    double sum = 0.0;
    for(size_t i = 1; i != averagedData.size(); ++i)
    {
        sum += 0.5 * (averagedData[i] + averagedData[i-1]);
    }
    return sum / static_cast<double>(numberOfBins - 1);
  }

  double averagedExcessChemicalPotential() const
  {
    std::vector<std::pair<double,double>> summedBlocks(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform (summedBlocks.begin(), summedBlocks.end(), bookKeepingDUdlambda[blockIndex].begin(), summedBlocks.begin(), pair_sum);
    }

    std::vector<double> averagedData(numberOfBins);
    std::transform(summedBlocks.begin(), summedBlocks.end(), averagedData.begin(),
                   [&](std::pair<double, double> &sample){return sample.first / std::max(1.0, sample.second);});

    double sum = 0.0;
    for(size_t i = 1; i != averagedData.size(); ++i)
    {
        sum += 0.5 * (averagedData[i] + averagedData[i-1]);
    }
    return sum / static_cast<double>(numberOfBins - 1);
  }

  std::pair<double, double> averageExcessChemicalPotential() const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedExcessChemicalPotential();

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = averagedExcessChemicalPotential(blockIndex) - average;
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
    double average = averagedExcessChemicalPotential()  + averagedIdealGasChemicalPotential(Beta);

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingDensity[blockIndex].second / bookKeepingDensity[0].second > 0.5)
      {
        double value = (averagedExcessChemicalPotential(blockIndex) + averagedIdealGasChemicalPotential(blockIndex, Beta)) - average;
        sumOfSquares += value * value;
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

  std::pair<double, double> averageFugacity(double Beta) const
  {
    double average = std::exp(Beta*(averagedExcessChemicalPotential()  + averagedIdealGasChemicalPotential(Beta))) / Beta;

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingDensity[blockIndex].second / bookKeepingDensity[0].second > 0.5)
      {
        double value = std::exp(Beta * (averagedExcessChemicalPotential(blockIndex) + averagedIdealGasChemicalPotential(blockIndex, Beta))) / Beta - average;
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
};
