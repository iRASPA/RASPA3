export module property_energy;

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
import energy_status;
import energy_status_inter;
import energy_status_intra;
import component;

inline std::pair<EnergyStatus, double> pair_sum(const std::pair<EnergyStatus, double> &lhs, const std::pair<EnergyStatus, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyEnergy
{
  PropertyEnergy(size_t numberOfBlocks, size_t numberOfComponents) :
      numberOfBlocks(numberOfBlocks),
      numberOfComponents(numberOfComponents),
      bookKeepingEnergyStatus(std::vector<std::pair<EnergyStatus, double>>(numberOfBlocks, std::make_pair(EnergyStatus(numberOfComponents), 0.0)))
  {
  }

  size_t numberOfBlocks;
  size_t numberOfComponents;
  std::vector<std::pair<EnergyStatus, double>> bookKeepingEnergyStatus;

  void resize(size_t newNumberOfComponents)
  {
      numberOfComponents = newNumberOfComponents;
      bookKeepingEnergyStatus = std::vector<std::pair<EnergyStatus, double>>(numberOfBlocks, std::make_pair(EnergyStatus(numberOfComponents), 0.0));
  }

  inline void addSample(size_t blockIndex, const EnergyStatus &energyStatus, const double &weight)
  {
    bookKeepingEnergyStatus[blockIndex].first += weight * energyStatus;
    bookKeepingEnergyStatus[blockIndex].second += weight;
  }

  //====================================================================================================================

  EnergyStatus averagedEnergy(size_t blockIndex) const
  {
    return bookKeepingEnergyStatus[blockIndex].first / bookKeepingEnergyStatus[blockIndex].second;
  }

  EnergyStatus averagedEnergy() const
  {
    std::pair<EnergyStatus,double> summedBlocks = std::accumulate (bookKeepingEnergyStatus.begin(), bookKeepingEnergyStatus.end(), std::make_pair(EnergyStatus(numberOfComponents), 0.0), pair_sum);
    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<EnergyStatus, EnergyStatus> averageEnergy() const
  {
    EnergyStatus average = averagedEnergy();

    // Use bins that are at least 90% filled for the computation of the error
    EnergyStatus sumOfSquares(numberOfComponents);
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingEnergyStatus[blockIndex].second / bookKeepingEnergyStatus[0].second > 0.5)
      {
        EnergyStatus value = averagedEnergy(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }
    EnergyStatus confidenceIntervalError(numberOfComponents);
    if(numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      EnergyStatus standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      EnergyStatus standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }
    return std::make_pair(average, confidenceIntervalError);
  }

  std::string writeAveragesStatistics(std::vector<Component>& components) const;
};
