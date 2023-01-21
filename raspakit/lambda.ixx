export module lambda;

import randomnumbers;
import averages;
import property_lambda_probability_histogram;
import property_dudlambda;

import <vector>;
import <iostream>;
import <cmath>;
import <string>;
import <sstream>;
import <format>;
import <algorithm>;
import <numeric>;
import <numbers>;

export struct Lambda
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

  Lambda(size_t numberOfBlocks, size_t numberOfBins = 11) : 
      numberOfBins(numberOfBins), 
      currentBin(0), 
      delta(1.0 / static_cast<double>(numberOfBins - 1)), 
      histogram(numberOfBins),
      biasFactor(numberOfBins),
      newHistogram(numberOfBlocks, numberOfBins),
      dUdlambdaBookKeeping(numberOfBlocks,numberOfBins)
  {
  }

  size_t numberOfBins;
  size_t currentBin;
  double delta;

  double WangLandauScalingFactor{ 1.0 };

  std::vector<double> histogram;
  std::vector<double> biasFactor;

  PropertyLambdaProbabilityHistogram newHistogram;
  PropertyDUDlambda dUdlambdaBookKeeping;

  inline double lambdaValue()
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

  inline void updateHistogram()
  {
    histogram[currentBin] += 1.0;
  }

  void sampleHistogram(size_t blockIndex, double density)
  {
    newHistogram.addSample(blockIndex, currentBin, density, weight());
  }

  void sampledUdLambdaHistogram(size_t blockIndex, double value)
  {
    dUdlambdaBookKeeping.addSample(blockIndex, currentBin, 2.0 * value);
  }

  inline double weight() const
  {
    return std::exp(-biasFactor[currentBin]);
  }

  void WangLandauIteration(Lambda::WangLandauPhase phase);
};
