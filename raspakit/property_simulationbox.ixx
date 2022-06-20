export module property_simulationbox;

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
import simulationbox;

inline std::pair<SimulationBox, double> pair_sum(const std::pair<SimulationBox, double> &lhs, const std::pair<SimulationBox, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertySimulationBox
{
  PropertySimulationBox(size_t numberOfBlocks) :
      numberOfBlocks(numberOfBlocks),
      bookKeepingSimulationBox(numberOfBlocks)
  {
  }

  size_t numberOfBlocks;
  std::vector<std::pair<SimulationBox, double>> bookKeepingSimulationBox;

  inline void addSample(size_t blockIndex, const SimulationBox &box, const double &weight)
  {
    bookKeepingSimulationBox[blockIndex].first += weight * box;
    bookKeepingSimulationBox[blockIndex].second += weight;
  }

  //====================================================================================================================

  SimulationBox averagedSimulationBox(size_t blockIndex) const
  {
    return bookKeepingSimulationBox[blockIndex].first / bookKeepingSimulationBox[blockIndex].second;
  }

  SimulationBox averagedSimulationBox() const
  {
    std::pair<SimulationBox,double> summedBlocks = std::accumulate (bookKeepingSimulationBox.begin(), bookKeepingSimulationBox.end(), std::make_pair(SimulationBox(), 0.0), pair_sum);
    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<SimulationBox, SimulationBox> averageSimulationBox() const
  {
    SimulationBox average = averagedSimulationBox();

    SimulationBox sumOfSquares{};
    size_t numberOfSamples = 0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingSimulationBox[blockIndex].second / bookKeepingSimulationBox[0].second > 0.5)
      { 
        SimulationBox value = averagedSimulationBox(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }

    SimulationBox confidenceIntervalError{};
    if(numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      SimulationBox standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      SimulationBox standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }
};
