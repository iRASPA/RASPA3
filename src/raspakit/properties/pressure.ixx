module;

export module property_pressure;

import std;

import archive;
import double3x3;
import averages;
import pressures;
import units;
import json;

inline std::pair<Pressures, double> pair_acc_pressure(const std::pair<Pressures, double> &lhs,
                                                      const std::pair<Pressures, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyPressure
{
  PropertyPressure() {};

  PropertyPressure(std::size_t numberOfBlocks)
      : numberOfBlocks(numberOfBlocks),
        bookKeepingPressure(numberOfBlocks)
  {
  }

  bool operator==(PropertyPressure const &) const = default;

  std::uint64_t versionNumber{1};
  std::size_t numberOfBlocks;
  std::vector<std::pair<Pressures, double>> bookKeepingPressure;

  inline void addSample(std::size_t blockIndex, double idealGasPressureValue, double3x3 excessPressureTensor,
                        double weight)
  {
    double3x3 ideal_gas_pressure_tensor = double3x3(idealGasPressureValue, idealGasPressureValue, idealGasPressureValue);
    double3x3 total_pressure_tensor = excessPressureTensor + ideal_gas_pressure_tensor;
    double excess_pressure = excessPressureTensor.trace() / 3.0;
    double total_pressure = idealGasPressureValue + excess_pressure;

    bookKeepingPressure[blockIndex].first.totalPressureTensor += weight * total_pressure_tensor;
    bookKeepingPressure[blockIndex].first.excessPressureTensor += weight * excessPressureTensor;
    bookKeepingPressure[blockIndex].first.idealGasPressureTensor += weight * ideal_gas_pressure_tensor;

    bookKeepingPressure[blockIndex].first.totalPressure += weight * total_pressure;
    bookKeepingPressure[blockIndex].first.excessPressure += weight * excess_pressure;
    bookKeepingPressure[blockIndex].first.idealGasPressure += weight * idealGasPressureValue;

    bookKeepingPressure[blockIndex].second += weight;
  }

  //====================================================================================================================


  Pressures averagedPressures(std::size_t blockIndex) const
  {
    return bookKeepingPressure[blockIndex].first / bookKeepingPressure[blockIndex].second;
  }

  Pressures averagedPressures() const
  {
    std::pair<Pressures, double> summedBlocks = std::accumulate(
        bookKeepingPressure.begin(), bookKeepingPressure.end(),
        std::make_pair(Pressures(), 0.0), pair_acc_pressure);
    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<Pressures, Pressures> result() const
  {
    Pressures average = averagedPressures();

    Pressures sumOfSquares{};
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingPressure[blockIndex].second / std::max(1.0, bookKeepingPressure[0].second) > 0.5)
      {
        Pressures value = averagedPressures(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    Pressures confidenceIntervalError{};
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      Pressures standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      Pressures standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  std::string writeAveragesStatistics() const;
  nlohmann::json jsonAveragesStatistics() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPressure &e);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPressure &e);
};
