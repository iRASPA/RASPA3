module;

export module property_partial_molar_properties;

import std;

import archive;
import averages;
import partial_molar_properties_data;
import component;
import json;
import units;

/**
 * \brief Block-averaged accumulator for partial molar internal energy and partial molar volume.
 *
 * Mirrors PropertyEnthalpy: raw fluctuation terms are accumulated per block, averaged, and
 * combined through the matrix fluctuation formula (PartialMolarPropertiesTerms::compositeProperty).
 * Confidence intervals are estimated from the spread of the block estimates.
 */
export struct PropertyPartialMolarProperties
{
  PropertyPartialMolarProperties() {};
  PropertyPartialMolarProperties(std::size_t numberOfBlocks, std::size_t numberOfComponents)
      : numberOfBlocks(numberOfBlocks),
        numberOfComponents(numberOfComponents),
        bookKeepingPartialMolarPropertiesTerms(std::vector<std::pair<PartialMolarPropertiesTerms, double>>(
            numberOfBlocks, std::make_pair(PartialMolarPropertiesTerms(numberOfComponents), 0.0)))
  {
  }

  bool operator==(PropertyPartialMolarProperties const &) const = default;

  std::uint64_t versionNumber{1};
  std::size_t numberOfBlocks;
  std::size_t numberOfComponents;
  std::vector<std::pair<PartialMolarPropertiesTerms, double>> bookKeepingPartialMolarPropertiesTerms;

  void resize(std::size_t newNumberOfComponents)
  {
    numberOfComponents = newNumberOfComponents;
    bookKeepingPartialMolarPropertiesTerms = std::vector<std::pair<PartialMolarPropertiesTerms, double>>(
        numberOfBlocks, std::make_pair(PartialMolarPropertiesTerms(numberOfComponents), 0.0));
  }

  inline void addSample(std::size_t blockIndex, const PartialMolarPropertiesTerms &terms, const double &weight)
  {
    bookKeepingPartialMolarPropertiesTerms[blockIndex].first += weight * terms;
    bookKeepingPartialMolarPropertiesTerms[blockIndex].second += weight;
  }

  //====================================================================================================================

  PartialMolarPropertiesData averagedProperties(std::size_t blockIndex) const
  {
    return (bookKeepingPartialMolarPropertiesTerms[blockIndex].first /
            std::max(1.0, bookKeepingPartialMolarPropertiesTerms[blockIndex].second))
        .compositeProperty();
  }

  PartialMolarPropertiesData averagedProperties() const
  {
    PartialMolarPropertiesData average(numberOfComponents);
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingPartialMolarPropertiesTerms[blockIndex].second /
              std::max(1.0, bookKeepingPartialMolarPropertiesTerms[0].second) >
          0.5)
      {
        average += averagedProperties(blockIndex);
        ++numberOfSamples;
      }
    }
    return (1.0 / static_cast<double>(std::max(1uz, numberOfSamples))) * average;
  }

  std::pair<PartialMolarPropertiesData, PartialMolarPropertiesData> averageProperties() const
  {
    PartialMolarPropertiesData average = averagedProperties();

    PartialMolarPropertiesData sumOfSquares(numberOfComponents);
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingPartialMolarPropertiesTerms[blockIndex].second /
              std::max(1.0, bookKeepingPartialMolarPropertiesTerms[0].second) >
          0.5)
      {
        PartialMolarPropertiesData value = averagedProperties(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }
    PartialMolarPropertiesData confidenceIntervalError(numberOfComponents);
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      PartialMolarPropertiesData standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      PartialMolarPropertiesData standardError =
          (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  std::string writeAveragesStatistics(std::vector<std::size_t> &swappableComponents,
                                      std::vector<Component> &components) const;
  nlohmann::json jsonAveragesStatistics(std::vector<std::size_t> &swappableComponents,
                                        std::vector<Component> &components) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPartialMolarProperties &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPartialMolarProperties &p);
};
