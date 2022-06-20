export module property_enthalpy;

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
import enthalpy_of_adsorption;

inline std::pair<double, double> pair_sum(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyEnthalpy
{
  PropertyEnthalpy(size_t numberOfBlocks, size_t numberOfComponents) :
      numberOfBlocks(numberOfBlocks),
      numberOfComponents(numberOfComponents),
      bookKeepingEnthalpyOfAdsorptionTerms(std::vector<std::pair<EnthalpyOfAdsorptionTerms, double>>(numberOfBlocks, std::make_pair(EnthalpyOfAdsorptionTerms(numberOfComponents), 0.0)))
  {
  }

  size_t numberOfBlocks;
  size_t numberOfComponents;
  std::vector<std::pair<EnthalpyOfAdsorptionTerms, double>> bookKeepingEnthalpyOfAdsorptionTerms;

  void resize(size_t newNumberOfComponents)
  {
      numberOfComponents = newNumberOfComponents;
      bookKeepingEnthalpyOfAdsorptionTerms = std::vector<std::pair<EnthalpyOfAdsorptionTerms, double>>(numberOfBlocks, std::make_pair(EnthalpyOfAdsorptionTerms(numberOfComponents), 0.0));
  }

  inline void addSample(size_t blockIndex, const EnthalpyOfAdsorptionTerms &terms, const double &weight)
  {
    bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].first += weight * terms;
    bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].second += weight;
  }

  //====================================================================================================================

  EnthalpyOfAdsorption averagedEnthalpy(size_t blockIndex) const
  {
    return (bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].first / bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].second).compositeProperty();
  }

  EnthalpyOfAdsorption averagedEnthalpy() const
  {
    EnthalpyOfAdsorption average(numberOfComponents);
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].second / bookKeepingEnthalpyOfAdsorptionTerms[0].second > 0.5)
      {
        average += averagedEnthalpy(blockIndex);
        ++numberOfSamples;
      }
    }
    return (1.0 / static_cast<double>(numberOfSamples)) * average;
  }

  std::pair<EnthalpyOfAdsorption, EnthalpyOfAdsorption> averageEnthalpy() const
  {
    EnthalpyOfAdsorption average = averagedEnthalpy();

    EnthalpyOfAdsorption sumOfSquares(numberOfComponents);
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].second / bookKeepingEnthalpyOfAdsorptionTerms[0].second > 0.5)
      {
        EnthalpyOfAdsorption value = averagedEnthalpy(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }
    EnthalpyOfAdsorption confidenceIntervalError(numberOfComponents);
    if(numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      EnthalpyOfAdsorption standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      EnthalpyOfAdsorption standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }
};
