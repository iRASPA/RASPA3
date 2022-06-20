export module property_lambda_probability_histogram;

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

export struct PropertyLambdaProbabilityHistogram
{
  PropertyLambdaProbabilityHistogram(size_t numberOfBlocks, size_t numberOfBins) :
      numberOfBlocks(numberOfBlocks),
      numberOfBins(numberOfBins),
      sumProperty(std::vector<std::vector<double>>(numberOfBlocks, std::vector<double>(numberOfBins))),
      sumWeightedDensity(std::vector<std::pair<double, double>>(numberOfBlocks))
  {
  }

  size_t numberOfBlocks;
  size_t numberOfBins;
  std::vector<std::vector<double>> sumProperty;
  std::vector<std::pair<double, double>> sumWeightedDensity;

  inline void addSample(size_t blockIndex, size_t binIndex, double density, double weight)
  {
    sumProperty[blockIndex][binIndex] += 1.0;
    sumWeightedDensity[blockIndex].first += weight * density;
    sumWeightedDensity[blockIndex].second += weight;
  }

  //====================================================================================================================

  double averagedDensity(size_t blockIndex) const
  {
    return sumWeightedDensity[blockIndex].first / sumWeightedDensity[blockIndex].second;
  }

  double averagedDensity() const
  {
    std::pair<double,double> summedBlocks{0.0, 0.0};
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
        summedBlocks.first += sumWeightedDensity[blockIndex].first;
        summedBlocks.second += sumWeightedDensity[blockIndex].second;
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
    return -std::log(1.0 / (sumWeightedDensity[blockIndex].first / sumWeightedDensity[blockIndex].second)) / Beta;
  }

  double averagedIdealGasChemicalPotential(double Beta) const
  {
    std::pair<double,double> summedBlocks{0.0, 0.0};
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
        summedBlocks.first += sumWeightedDensity[blockIndex].first;
        summedBlocks.second += sumWeightedDensity[blockIndex].second;
    }

    return -std::log(1.0 / (summedBlocks.first / summedBlocks.second)) / Beta;
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

  std::vector<double> averagedProbabilityHistogram(size_t blockIndex) const
  {
    std::vector<double> averagedData(numberOfBins);
    std::transform(sumProperty[blockIndex].begin(), sumProperty[blockIndex].end(), averagedData.begin(),
                   [&](const double &sample){return sample;});
    return averagedData;
  }

  std::vector<double> averagedProbabilityHistogram() const
  {
    std::vector<double> summedBlocks(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform(summedBlocks.begin(), summedBlocks.end(), sumProperty[blockIndex].begin(), summedBlocks.begin(), 
                     [](const double & a, const double & b){ return a + b; });
    }

    return summedBlocks;
  }

  std::pair<std::vector<double>, std::vector<double>> averageProbabilityHistogram() const
  {
    size_t degreesOfFreedom = numberOfBlocks - 1;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    std::vector<double> average = averagedProbabilityHistogram();

    std::vector<double> sumOfSquares(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::vector<double> blockAverage = averagedProbabilityHistogram(blockIndex);
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
/*
  std::vector<double> averagedNormalizedProbabilityHistogram(size_t blockIndex) const
  {
    std::vector<double> averagedData(numberOfBins);
    double count = std::accumulate(sumProperty[blockIndex].cbegin(), sumProperty[blockIndex].cend(), 0.0,
                            [&](double acc, const std::pair<double, double> &sample){return acc + sample.second;});
    std::transform(sumProperty[blockIndex].cbegin(), sumProperty[blockIndex].cend(), averagedData.begin(),
                   [&](const std::pair<double, double> &sample){return sample.first / std::max(1.0, count);});
    return averagedData;
  }

  std::vector<double> averagedNormalizedProbabilityHistogram() const
  {
    std::vector<std::pair<double,double>> summedBlocks(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform (summedBlocks.begin(), summedBlocks.end(), sumProperty[blockIndex].begin(), summedBlocks.begin(), pair_sum<double,double>());
    }
    double count = std::accumulate(summedBlocks.cbegin(), summedBlocks.cend(), 0.0,
                            [&](double acc, const std::pair<double, double> &sample){return acc + sample.second;});

    std::vector<double> averagedData(numberOfBins);
    std::transform(summedBlocks.begin(), summedBlocks.end(), averagedData.begin(),
                   [&](std::pair<double, double> &sample){return sample.first / std::max(1.0, count);});
    return averagedData;
  }

  std::pair<std::vector<double>, std::vector<double>> averageNormalizedProbabilityHistogram() const
  {
    size_t degreesOfFreedom = numberOfBlocks - 1;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    std::vector<double> average = averagedProbabilityHistogram();

    std::vector<double> sumOfSquares(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::vector<double> blockAverage = averagedProbabilityHistogram(blockIndex);
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
*/
  //===============

  std::vector<double> averagedLandauFreeEnergyHistogram(size_t blockIndex, double Beta) const
  {
    std::vector<double> averagedData(numberOfBins);
    std::transform(sumProperty[blockIndex].cbegin(), sumProperty[blockIndex].cend(), averagedData.begin(),
                   [&](const double &sample){return -std::log(sample) / Beta;});
    return averagedData;
  }

  std::vector<double> averagedLandauFreeEnergyHistogram(double Beta) const
  {
    // sum all blocks into one
    std::vector<double> summedBlocks(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform(summedBlocks.begin(), summedBlocks.end(), sumProperty[blockIndex].begin(), summedBlocks.begin(), 
                     [](const double & a, const double & b){ return a + b; });
    }

    std::vector<double> averagedData(numberOfBins);
    std::transform(summedBlocks.cbegin(), summedBlocks.cend(), averagedData.begin(),
                   [&](const double &sample){return -std::log(sample) / Beta;});
    return averagedData;
  }

  std::pair<std::vector<double>, std::vector<double>> averageLandauFreeEnergyHistogram(double Beta) const
  {
    size_t degreesOfFreedom = numberOfBlocks - 1;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    std::vector<double> average = averagedLandauFreeEnergyHistogram(Beta);

    std::vector<double> sumOfSquares(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::vector<double> blockAverage = averagedLandauFreeEnergyHistogram(blockIndex, Beta);
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

  double averagedExcessChemicalPotential(size_t blockIndex, double Beta) const
  {
    size_t lastBin = numberOfBins - 1;
    double averagedProbabilityFirstValue = sumProperty[blockIndex][0];
    double averagedProbabilityLastValue = sumProperty[blockIndex][lastBin];
    return -std::log(averagedProbabilityLastValue / averagedProbabilityFirstValue) / Beta;
  }

  double averagedExcessChemicalPotential(double Beta) const
  {
    std::vector<double> summedBlocks(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform (summedBlocks.begin(), summedBlocks.end(), sumProperty[blockIndex].begin(), summedBlocks.begin(),
                     [](const double & a, const double & b){ return a + b; });
    }

    size_t lastBin = numberOfBins - 1;
    double averagedProbabilityFirstValue = summedBlocks[0];
    double averagedProbabilityLastValue = summedBlocks[lastBin];
    return -std::log(averagedProbabilityLastValue / averagedProbabilityFirstValue) / Beta;
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

  std::pair<double, double> averageTotalChemicalPotential(double Beta, double bias) const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedExcessChemicalPotential(Beta)  + averagedIdealGasChemicalPotential(Beta) + bias;

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = (averagedExcessChemicalPotential(blockIndex, Beta) + averagedIdealGasChemicalPotential(blockIndex, Beta) + bias) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }

  std::pair<double, double> averageFugacity(double Beta, double bias) const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = std::exp(Beta*(averagedExcessChemicalPotential(Beta)  + averagedIdealGasChemicalPotential(Beta) + bias)) / Beta;

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = std::exp(Beta * (averagedExcessChemicalPotential(blockIndex, Beta) + averagedIdealGasChemicalPotential(blockIndex, Beta) + bias)) / Beta - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }
};
