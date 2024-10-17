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

export module property_enthalpy;

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
import loadings;
import enthalpy_of_adsorption;
import component;
import json;
import units;

inline std::pair<double, double> pair_sum(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyEnthalpy
{
  PropertyEnthalpy() {};
  PropertyEnthalpy(size_t numberOfBlocks, size_t numberOfComponents)
      : numberOfBlocks(numberOfBlocks),
        numberOfComponents(numberOfComponents),
        bookKeepingEnthalpyOfAdsorptionTerms(std::vector<std::pair<EnthalpyOfAdsorptionTerms, double>>(
            numberOfBlocks, std::make_pair(EnthalpyOfAdsorptionTerms(numberOfComponents), 0.0)))
  {
  }

  bool operator==(PropertyEnthalpy const &) const = default;

  uint64_t versionNumber{1};
  size_t numberOfBlocks;
  size_t numberOfComponents;
  std::vector<std::pair<EnthalpyOfAdsorptionTerms, double>> bookKeepingEnthalpyOfAdsorptionTerms;

  void resize(size_t newNumberOfComponents)
  {
    numberOfComponents = newNumberOfComponents;
    bookKeepingEnthalpyOfAdsorptionTerms = std::vector<std::pair<EnthalpyOfAdsorptionTerms, double>>(
        numberOfBlocks, std::make_pair(EnthalpyOfAdsorptionTerms(numberOfComponents), 0.0));
  }

  inline void addSample(size_t blockIndex, const EnthalpyOfAdsorptionTerms &terms, const double &weight)
  {
    bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].first += weight * terms;
    bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].second += weight;
  }

  //====================================================================================================================

  EnthalpyOfAdsorption averagedEnthalpy(size_t blockIndex) const
  {
    return (bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].first /
            bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].second)
        .compositeProperty();
  }

  EnthalpyOfAdsorption averagedEnthalpy() const
  {
    EnthalpyOfAdsorption average(numberOfComponents);
    size_t numberOfSamples = 0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].second / bookKeepingEnthalpyOfAdsorptionTerms[0].second >
          0.5)
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
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingEnthalpyOfAdsorptionTerms[blockIndex].second / bookKeepingEnthalpyOfAdsorptionTerms[0].second >
          0.5)
      {
        EnthalpyOfAdsorption value = averagedEnthalpy(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }
    EnthalpyOfAdsorption confidenceIntervalError(numberOfComponents);
    if (numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      EnthalpyOfAdsorption standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      EnthalpyOfAdsorption standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  std::vector<double> blockEnthalpies(size_t &k, double &idealGasTerm) const
  {
    std::vector<double> enthalpy(numberOfBlocks);

    std::transform(
        bookKeepingEnthalpyOfAdsorptionTerms.begin(), bookKeepingEnthalpyOfAdsorptionTerms.end(), enthalpy.begin(),
        [&](const std::pair<EnthalpyOfAdsorptionTerms, double> &block) {
          return Units::EnergyToKelvin * ((block.first / block.second).compositeProperty().values[k] - idealGasTerm);
        });
    return enthalpy;
  }

  std::string writeAveragesStatistics(std::vector<size_t> &swappableComponents,
                                      std::vector<Component> &components) const;
  nlohmann::json jsonAveragesStatistics(std::vector<size_t> &swappableComponents,
                                        std::vector<Component> &components) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyEnthalpy &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEnthalpy &p);
};
