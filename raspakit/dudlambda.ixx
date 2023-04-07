export module dudlambda;

import randomnumbers;
import averages;
import property_dudlambda;
import property_lambda_probability_histogram;

import <vector>;
import <iostream>;
import <cmath>;
import <string>;
import <sstream>;
import <format>;
import <algorithm>;
import <numeric>;
import <numbers>;

export struct dUdLambda
{
  enum class Scheme : size_t
  {
    Discrete = 0,
    Random = 1,
    Continuous = 2
  };

  enum class WangLandauPhase : size_t
  {
    Initialize = 0,
    Sample = 1,
    AdjustBiasingFactors = 2,
    Finalize = 3
  };

  dUdLambda(size_t numberOfBlocks, size_t numberOfBins = 11) :
    numberOfBins(numberOfBins),
    currentBin(0),
    delta(1.0 / static_cast<double>(numberOfBins - 1)),
    histogram(numberOfBins),
    biasFactor(numberOfBins),
    dUdlambdaBookKeeping(numberOfBlocks, numberOfBins)
  {
  }

  size_t numberOfBins;
  size_t currentBin;
  double delta;

  double WangLandauScalingFactor{ 1.0 };

  std::vector<double> histogram;
  std::vector<double> biasFactor;

  PropertyDUDlambda dUdlambdaBookKeeping;

  inline double lambdaValue() const
  {
    return static_cast<double>(currentBin) * delta;
  }

  inline int selectNewBin()
  {
    return static_cast<int>(currentBin) + static_cast<int>(5 * 2.0 * (RandomNumber::Uniform() - 0.5));
  }

  inline void setCurrentBin(size_t index)
  {
    currentBin = index;
  }

  void sampledUdLambdaHistogram(size_t blockIndex, double value)
  {
    dUdlambdaBookKeeping.addSample(blockIndex, currentBin, 2.0 * value);
  }

  inline double weight() const
  {
    return std::exp(-biasFactor[currentBin]);
  }

  void WangLandauIteration(dUdLambda::WangLandauPhase phase);
  std::string writeAveragesStatistics(double beta) const;
};
