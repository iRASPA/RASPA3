module;

export module property_widom;

import std;

import archive;
import int3;
import averages;
import widom_data;


inline std::pair<double, double> pair_acc_widom(const std::pair<double, double> &lhs,
                                                const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

inline std::pair<WidomData, double> pair_acc_widom2(const std::pair<WidomData, double> &lhs,
                                                    const std::pair<WidomData, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyWidom
{
  PropertyWidom();

  PropertyWidom(std::size_t numberOfBlocks)
      : numberOfBlocks(numberOfBlocks), 
        bookKeepingRosenbluthWeight(numberOfBlocks),
        bookKeepingFugacity(numberOfBlocks)
  {
  }

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks;
  std::vector<std::pair<double, double>> bookKeepingRosenbluthWeight;
  std::vector<std::pair<WidomData, double>> bookKeepingFugacity;

  std::string writeAveragesRosenbluthWeightStatistics(double temperature, double volume,
                                                      std::optional<double> frameworkMass,
                                                      std::optional<int3> number_of_unit_cells) const;
  std::string writeAveragesChemicalPotentialStatistics(double beta, std::optional<double> imposedChemicalPotential,
                                                       std::optional<double> imposedFugacity) const;

  inline void addWidomSample(std::size_t blockIndex, double RosenbluthValue, std::size_t N, double V, double weight)
  {
    bookKeepingRosenbluthWeight[blockIndex].first += weight * RosenbluthValue;
    bookKeepingRosenbluthWeight[blockIndex].second += weight;

    bookKeepingFugacity[blockIndex].first.total += 0.0;
    bookKeepingFugacity[blockIndex].first.excess += weight * RosenbluthValue;
    bookKeepingFugacity[blockIndex].first.idealGas += weight * static_cast<double>(N) / V;

    bookKeepingFugacity[blockIndex].second += weight;
  }

  //====================================================================================================================

  double averagedRosenbluthWeight() const
  {
    std::pair<double, double> summedBlocks = std::accumulate(
        bookKeepingRosenbluthWeight.begin(), bookKeepingRosenbluthWeight.end(),
        std::make_pair(0.0, 0.0), pair_acc_widom);

    return summedBlocks.first / summedBlocks.second;
  }


  double averagedRosenbluthWeight(std::size_t blockIndex) const
  {
    return bookKeepingRosenbluthWeight[blockIndex].first / bookKeepingRosenbluthWeight[blockIndex].second;
  }

  std::pair<double, double> result() const
  {
    double average = averagedRosenbluthWeight();

    double sumOfSquares{};
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingRosenbluthWeight[blockIndex].second / bookKeepingRosenbluthWeight[0].second > 0.5)
      {
        double value = averagedRosenbluthWeight(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError{};
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
 
  WidomData averagedChemicalPotential(double beta) const
  {
    std::pair<WidomData, double> summedBlocks = std::accumulate(
        bookKeepingFugacity.begin(), bookKeepingFugacity.end(),
        std::make_pair(WidomData(), 0.0), pair_acc_widom2);

    return WidomData(-(1.0 / beta) * std::log(summedBlocks.first.excess / summedBlocks.second) +
                      (1.0 / beta) * std::log(summedBlocks.first.idealGas / summedBlocks.second),
                     -(1.0 / beta) * std::log(summedBlocks.first.excess / summedBlocks.second),
                      (1.0 / beta) * std::log(summedBlocks.first.idealGas / summedBlocks.second));
  }


  WidomData averagedChemicalPotential(std::size_t blockIndex, double beta) const
  {
    return WidomData(-(1.0 / beta) * std::log(bookKeepingFugacity[blockIndex].first.excess / bookKeepingFugacity[blockIndex].second) +
                      (1.0 / beta) * std::log(bookKeepingFugacity[blockIndex].first.idealGas / bookKeepingFugacity[blockIndex].second),
                     -(1.0 / beta) * std::log(bookKeepingFugacity[blockIndex].first.excess / bookKeepingFugacity[blockIndex].second),
                      (1.0 / beta) * std::log(bookKeepingFugacity[blockIndex].first.idealGas / bookKeepingFugacity[blockIndex].second));
  }

  std::pair<WidomData, WidomData> chemicalPotentialResult(double beta) const
  {
    WidomData average = averagedChemicalPotential(beta);

    WidomData sumOfSquares{};
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingFugacity[blockIndex].second / bookKeepingFugacity[0].second > 0.5)
      {
        WidomData value = averagedChemicalPotential(blockIndex, beta) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    WidomData confidenceIntervalError{};
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      WidomData standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      WidomData standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedFugacity(double beta) const
  {
    return std::exp(beta * averagedChemicalPotential(beta).total) / beta;
  }


  double averagedFugacity(std::size_t blockIndex, double beta) const
  {
    return std::exp(beta * averagedChemicalPotential(blockIndex, beta).total) / beta;
  }


  std::pair<double, double> fugacityResult(double beta) const
  {
    double average = averagedFugacity(beta);

    double sumOfSquares{};
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingFugacity[blockIndex].second / bookKeepingFugacity[0].second > 0.5)
      {
        double value = averagedFugacity(blockIndex, beta) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    double confidenceIntervalError{};
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

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyWidom &w);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyWidom &w);
};
