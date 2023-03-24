export module property_loading;

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
import loadings;
import component;

inline std::pair<Loadings, double> pair_sum(const std::pair<Loadings, double> &lhs, const std::pair<Loadings, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyLoading
{
  PropertyLoading(size_t numberOfBlocks, size_t numberOfComponents) :
      numberOfBlocks(numberOfBlocks),
      numberOfComponents(numberOfComponents),
      bookKeepingLoadings(std::vector<std::pair<Loadings, double>>(numberOfBlocks, std::make_pair(Loadings(numberOfComponents), 0.0)))
  {
  }

  size_t numberOfBlocks;
  size_t numberOfComponents;
  std::vector<std::pair<Loadings, double>> bookKeepingLoadings;

  void resize(size_t newNumberOfComponents)
  {
    numberOfComponents = newNumberOfComponents;
    bookKeepingLoadings = std::vector<std::pair<Loadings, double>>(numberOfBlocks, std::make_pair(Loadings(numberOfComponents), 0.0));
  }

  inline void addSample(size_t blockIndex, const Loadings &loading, const double &weight)
  {
    bookKeepingLoadings[blockIndex].first += weight * loading;
    bookKeepingLoadings[blockIndex].second += weight;
  }

  //====================================================================================================================

  Loadings averagedLoading(size_t blockIndex) const
  {
    return bookKeepingLoadings[blockIndex].first / bookKeepingLoadings[blockIndex].second;
  }

  Loadings averagedLoading() const
  {
    std::pair<Loadings,double> summedBlocks = std::accumulate (bookKeepingLoadings.begin(), bookKeepingLoadings.end(), std::make_pair(Loadings(numberOfComponents), 0.0), pair_sum);
    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<Loadings, Loadings> averageLoading() const
  {
    Loadings average = averagedLoading();

    Loadings sumOfSquares(numberOfComponents);
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingLoadings[blockIndex].second / bookKeepingLoadings[0].second > 0.5)
      {
        Loadings value = averagedLoading(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }
    Loadings confidenceIntervalError(numberOfComponents);
    if(numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      Loadings standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      Loadings standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  std::string writeAveragesStatistics(std::vector<Component> components, std::optional<double> frameworkMass) const;
};
