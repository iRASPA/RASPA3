export module property_lambda_probability_histogram;

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

import double3;
import randomnumbers;
import archive;

import averages;

inline std::pair<double, double> pair_sum(const std::pair<double, double> &lhs, const std::pair<double, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

inline std::pair<double3, double> 
pair_sum_double3(const std::pair<double3, double>& lhs, const std::pair<double3, double>& rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyLambdaProbabilityHistogram
{
  enum class WangLandauPhase : size_t
  {
    Initialize = 0,
    Sample = 1,
    AdjustBiasingFactors = 2,
    Finalize = 3
  };

  PropertyLambdaProbabilityHistogram() {};

  PropertyLambdaProbabilityHistogram(size_t numberOfBlocks, size_t numberOfBins) :
      numberOfBlocks(numberOfBlocks),
      numberOfBins(numberOfBins),
      jump_bins(numberOfBins/2),
      currentBin(0),
      delta(1.0 / static_cast<double>(numberOfBins - 1)),
      histogram(numberOfBins),
      biasFactor(numberOfBins),
      bookKeepingLambda(std::vector<std::vector<double>>(numberOfBlocks, std::vector<double>(numberOfBins))),
      bookKeepingDensity(std::vector<std::pair<double, double>>(numberOfBlocks)),
      computeDUdlambda(false),
      bookKeepingDUdlambda(std::vector<std::vector<std::pair<double3, double>>>(numberOfBlocks, 
                           std::vector<std::pair<double3, double>>(numberOfBins)))
  {
  }

  bool operator==(PropertyLambdaProbabilityHistogram const&) const = default;

  uint64_t versionNumber{ 1 };
  size_t numberOfBlocks;

  size_t numberOfBins;
  size_t jump_bins;
  size_t currentBin;
  double delta;

  double WangLandauScalingFactor{ 1.0 };

  std::vector<double> histogram;
  std::vector<double> biasFactor;

  // lambda-histogram 
  std::vector<std::vector<double>> bookKeepingLambda;

  // first: sm of weight * density, second: sum of weights
  std::vector<std::pair<double, double>> bookKeepingDensity;

  // dU/dlambda-histogram 
  bool computeDUdlambda;
  std::vector<std::vector<std::pair<double3, double>>> bookKeepingDUdlambda;

  // fractional molecule occupancy
  size_t occupancyCount{ 0 };
  size_t occupancyTotal{ 0 };  

  void clear()
  {
    std::fill(histogram.begin(), histogram.end(), 0.0);

    for (size_t i = 0; i < numberOfBlocks; ++i)
    {
      bookKeepingDensity[i] = std::make_pair<double, double>(0.0, 0.0);
      std::fill(bookKeepingLambda[i].begin(), bookKeepingLambda[i].end(), 0.0);
      std::fill(bookKeepingDUdlambda[i].begin(), bookKeepingDUdlambda[i].end(), 
                std::make_pair<double3, double>(double3(0.0, 0.0, 0.0), 0.0));
    }

    occupancyCount = 0;
    occupancyTotal = 0;
  }

  inline double lambdaValue() const
  {
    return static_cast<double>(currentBin) * delta;
  }

  inline int selectNewBin(RandomNumber &random) const
  {
    return static_cast<int>(currentBin) + 
           static_cast<int>(static_cast<double>(jump_bins) * 2.0 * (random.uniform() - 0.5));
  }

  inline int selectNewBin(RandomNumber &random, double scale) const
  {
    return static_cast<int>(currentBin) + 
           static_cast<int>(scale * static_cast<double>(numberOfBins) * 2.0 * (random.uniform() - 0.5));
  }

  inline void setCurrentBin(size_t index)
  {
    currentBin = index;
  }

  inline void updateHistogram()
  {
    histogram[currentBin] += 1.0;
  }

  void sampleOccupancy(bool containsTheFractionalMolecule)
  {
    ++occupancyTotal;
    if (containsTheFractionalMolecule)
    {
      ++occupancyCount;
    }
  }

  double occupancy() const
  {
    return static_cast<double>(occupancyCount) / static_cast<double>(std::max(size_t{1}, occupancyTotal));
  }

  void normalize(double normalizationFactor)
  {
    for (double& bias : biasFactor)
    {
      bias -= normalizationFactor;
    }
  }

  void sampleHistogram(size_t blockIndex, double density, double dUdlambda, 
                       bool containsTheFractionalMolecule, double w)
  {
    if(containsTheFractionalMolecule)
    {
      bookKeepingLambda[blockIndex][currentBin] += 1.0;
      
      bookKeepingDUdlambda[blockIndex][currentBin].first.x += dUdlambda;
      bookKeepingDUdlambda[blockIndex][currentBin].second += 1.0;
    }

    bookKeepingDensity[blockIndex].first += w * density;
    bookKeepingDensity[blockIndex].second += w;    
  }

  inline double weight() const
  {
    return std::exp(-biasFactor[currentBin]);
  }

  void WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase phase, 
                           bool containsTheFractionalMolecule, double value = 1.0);

  std::string writeAveragesStatistics(double beta, std::optional<double> imposedChemicalPotential, 
                                      std::optional<double> imposedFugacity) const;

  std::string writeDUdLambdaStatistics(double beta, std::optional<double> imposedChemicalPotential, 
                                       std::optional<double> imposedFugacity) const;

  //====================================================================================================================

  double averagedDensity(size_t blockIndex) const
  {
    return bookKeepingDensity[blockIndex].first / bookKeepingDensity[blockIndex].second;
  }

  double averagedDensity() const
  {
    std::pair<double,double> summedBlocks{0.0, 0.0};
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      summedBlocks.first += bookKeepingDensity[blockIndex].first;
      summedBlocks.second += bookKeepingDensity[blockIndex].second;
    }

    return summedBlocks.first / summedBlocks.second;
  }

  std::pair<double, double> averageDensity() const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedDensity();

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = averagedDensity(blockIndex) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }


  //====================================================================================================================

  double averagedIdealGasChemicalPotential(size_t blockIndex, double beta) const
  {
    return -std::log(1.0 / (bookKeepingDensity[blockIndex].first / bookKeepingDensity[blockIndex].second)) / beta;
  }

  double averagedIdealGasChemicalPotential(double beta) const
  {
    std::pair<double,double> summedBlocks{0.0, 0.0};
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      summedBlocks.first += bookKeepingDensity[blockIndex].first;
      summedBlocks.second += bookKeepingDensity[blockIndex].second;
    }

    return -std::log(1.0 / (summedBlocks.first / summedBlocks.second)) / beta;
  }

  std::pair<double, double> averageIdealGasChemicalPotential(double beta) const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedIdealGasChemicalPotential(beta);

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = averagedIdealGasChemicalPotential(blockIndex, beta) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }


  //====================================================================================================================

  std::vector<double> averagedProbabilityHistogram(size_t blockIndex) const
  {
    std::vector<double> averagedData(numberOfBins);
    std::transform(bookKeepingLambda[blockIndex].begin(), bookKeepingLambda[blockIndex].end(), averagedData.begin(),
                   [&](const double &sample){return sample;});
    return averagedData;
  }

  std::vector<double> averagedProbabilityHistogram() const
  {
    std::vector<double> summedBlocks(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform(summedBlocks.begin(), summedBlocks.end(), bookKeepingLambda[blockIndex].begin(), summedBlocks.begin(),
                     [](const double & a, const double & b){ return a + b; });
    }

    std::vector<double> average(numberOfBins);
    std::transform(summedBlocks.begin(), summedBlocks.end(), average.begin(),
                 [&](const double &sample){return sample / static_cast<double>(numberOfBlocks);});

    return average;
  }

  std::pair<std::vector<double>, std::vector<double>> averageProbabilityHistogram() const
  {
    size_t degreesOfFreedom = numberOfBlocks - 1;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    std::vector<double> average = averagedProbabilityHistogram();

    std::vector<double> sumOfSquares(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::vector<double> blockAverage = averagedProbabilityHistogram(blockIndex);
      for(size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
      {
        double value = blockAverage[binIndex] - average[binIndex];
        sumOfSquares[binIndex] += value * value;
      }
    }
    std::vector<double> standardDeviation(numberOfBins);
    std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(),
                   [&](const double &sumofsquares)
                   {return std::sqrt(sumofsquares / static_cast<double>(degreesOfFreedom));});

    std::vector<double> standardError(numberOfBins);
    std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
                   [&](const double &sigma){return sigma / sqrt(static_cast<double>(numberOfBlocks));});

    std::vector<double> confidenceIntervalError(numberOfBins);
    std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
                   [&](const double &error){return intermediateStandardNormalDeviate * error;});

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================


  std::vector<double> averagedLandauFreeEnergyHistogram(size_t blockIndex, double beta) const
  {
    std::vector<double> averagedData(numberOfBins);
    std::transform(bookKeepingLambda[blockIndex].cbegin(), bookKeepingLambda[blockIndex].cend(), averagedData.begin(),
                   [&](const double &sample){return -std::log(sample) / beta;});
    return averagedData;
  }

  std::vector<double> averagedLandauFreeEnergyHistogram(double beta) const
  {
    // sum all blocks into one
    std::vector<double> summedBlocks(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform(summedBlocks.begin(), summedBlocks.end(), 
                     bookKeepingLambda[blockIndex].begin(), summedBlocks.begin(),
                     [](const double & a, const double & b){ return a + b; });
    }

    std::vector<double> averagedData(numberOfBins);
    std::transform(summedBlocks.cbegin(), summedBlocks.cend(), averagedData.begin(),
                   [&](const double &sample){return -std::log(sample / static_cast<double>(numberOfBlocks)) / beta;});
    return averagedData;
  }

  std::pair<std::vector<double>, std::vector<double>> averageLandauFreeEnergyHistogram(double beta) const
  {
    size_t degreesOfFreedom = numberOfBlocks - 1;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    std::vector<double> average = averagedLandauFreeEnergyHistogram(beta);

    std::vector<double> sumOfSquares(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::vector<double> blockAverage = averagedLandauFreeEnergyHistogram(blockIndex, beta);
      for(size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
      {
        double value = blockAverage[binIndex] - average[binIndex];
        sumOfSquares[binIndex] += value * value;
      }
    }
    std::vector<double> standardDeviation(numberOfBins);
    std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(),
                   [&](const double &sumofsquares)
                   {return std::sqrt(sumofsquares / static_cast<double>(degreesOfFreedom));});

    std::vector<double> standardError(numberOfBins);
    std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
                   [&](const double &sigma){return sigma / sqrt(static_cast<double>(numberOfBlocks));});

    std::vector<double> confidenceIntervalError(numberOfBins);
    std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
                   [&](const double &error){return intermediateStandardNormalDeviate * error;});

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  double averagedExcessChemicalPotential(size_t blockIndex, double beta) const
  {
    size_t lastBin = numberOfBins - 1;
    double averagedProbabilityFirstValue = bookKeepingLambda[blockIndex][0];
    double averagedProbabilityLastValue = bookKeepingLambda[blockIndex][lastBin];
    return -std::log(averagedProbabilityLastValue / averagedProbabilityFirstValue) / beta;
  }

  double averagedExcessChemicalPotential(double beta) const
  {
    std::vector<double> summedBlocks(numberOfBins);
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform(summedBlocks.begin(), summedBlocks.end(), 
                     bookKeepingLambda[blockIndex].begin(), summedBlocks.begin(),
                     [](const double & a, const double & b){ return a + b; });
    }

    size_t lastBin = numberOfBins - 1;
    double averagedProbabilityFirstValue = summedBlocks[0];
    double averagedProbabilityLastValue = summedBlocks[lastBin];
    return -std::log(averagedProbabilityLastValue / averagedProbabilityFirstValue) / beta;
  }

  std::pair<double, double> averageExcessChemicalPotential(double beta) const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedExcessChemicalPotential(beta);

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = averagedExcessChemicalPotential(blockIndex, beta) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  std::pair<double, double> averageTotalChemicalPotential(double beta, double bias) const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedExcessChemicalPotential(beta)  + averagedIdealGasChemicalPotential(beta) + bias;

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = (averagedExcessChemicalPotential(blockIndex, beta) + 
                      averagedIdealGasChemicalPotential(blockIndex, beta) + bias) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }

  std::pair<double, double> averageFugacity(double beta, double bias) const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = std::exp(beta*(averagedExcessChemicalPotential(beta)  + 
                                    averagedIdealGasChemicalPotential(beta) + bias)) / beta;

    double sumOfSquares = 0.0;
    for(size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = std::exp(beta * (averagedExcessChemicalPotential(blockIndex, beta) + 
                                      averagedIdealGasChemicalPotential(blockIndex, beta) + bias)) / beta - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }

  //====================================================================================================================

  std::vector<double3> averagedDUdlambda(size_t blockIndex) const
  {
    std::vector<double3> averagedData(numberOfBins);
    std::transform(bookKeepingDUdlambda[blockIndex].cbegin(), bookKeepingDUdlambda[blockIndex].cend(), 
      averagedData.begin(),
      [&](const std::pair<double3, double>& sample) {return sample.first / std::max(1.0, sample.second); });
    return averagedData;
  }

  std::vector<double3> averagedDUdlambda() const
  {
    std::vector<std::pair<double3, double>> summedBlocks(numberOfBins);
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform(summedBlocks.begin(), summedBlocks.end(), bookKeepingDUdlambda[blockIndex].begin(), 
                     summedBlocks.begin(), pair_sum_double3);
    }

    std::vector<double3> averagedData(numberOfBins);
    std::transform(summedBlocks.begin(), summedBlocks.end(), averagedData.begin(),
      [&](std::pair<double3, double>& sample) {return sample.first / std::max(1.0, sample.second); });
    return averagedData;
  }

  std::pair<std::vector<double3>, std::vector<double3>> averageDuDlambda() const
  {
    size_t degreesOfFreedom = numberOfBlocks - 1;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    std::vector<double3> average = averagedDUdlambda();

    std::vector<double3> sumOfSquares(numberOfBins);
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::vector<double3> blockAverage = averagedDUdlambda(blockIndex);
      for (size_t binIndex = 0; binIndex != numberOfBins; ++binIndex)
      {
        double3 value = blockAverage[binIndex] - average[binIndex];
        sumOfSquares[binIndex] += value * value;
      }
    }
    std::vector<double3> standardDeviation(numberOfBins);
    std::transform(sumOfSquares.cbegin(), sumOfSquares.cend(), standardDeviation.begin(),
      [&](const double3& sumofsquares) {return sqrt(sumofsquares / static_cast<double>(degreesOfFreedom)); });

    std::vector<double3> standardError(numberOfBins);
    std::transform(standardDeviation.cbegin(), standardDeviation.cend(), standardError.begin(),
      [&](const double3& sigma) {return sigma / sqrt(static_cast<double>(numberOfBlocks)); });

    std::vector<double3> confidenceIntervalError(numberOfBins);
    std::transform(standardError.cbegin(), standardError.cend(), confidenceIntervalError.begin(),
      [&](const double3& error) {return intermediateStandardNormalDeviate * error; });

    return std::make_pair(average, confidenceIntervalError);
  }


  //====================================================================================================================

  double averagedExcessChemicalPotentialDUdlambda(size_t blockIndex) const
  {
    std::vector<double> averagedData(numberOfBins);
    std::transform(bookKeepingDUdlambda[blockIndex].begin(), bookKeepingDUdlambda[blockIndex].end(), 
      averagedData.begin(),
      [&](const std::pair<double3, double>& sample) 
      {return (sample.first.x + sample.first.y + sample.first.z) / std::max(1.0, sample.second); });

    // trapezoidal rule: https://en.wikipedia.org/wiki/Trapezoidal_rule
    // Calculating result
    double res = 0;
    for (size_t i = 0; i < averagedData.size(); i++)
    {
      if (i == 0 || i == averagedData.size() - 1)
        res += averagedData[i];
      else if (i % 2 != 0)
        res += 4 * averagedData[i];
      else
        res += 2 * averagedData[i];
    }
    res = res * (delta / 3.0);
    return res;

    //double sum = 0.0;
    //for(size_t i = 1; i != averagedData.size(); ++i)
    //{
    //    sum += 0.5 * (averagedData[i] + averagedData[i-1]);
    //}
    //return sum / static_cast<double>(numberOfBins - 1);
  }

  double averagedExcessChemicalPotentialDUdlambda() const
  {
    std::vector<std::pair<double3, double>> summedBlocks(numberOfBins);
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      std::transform(summedBlocks.begin(), summedBlocks.end(), bookKeepingDUdlambda[blockIndex].begin(), 
                     summedBlocks.begin(), pair_sum_double3);
    }

    std::vector<double> averagedData(numberOfBins);
    std::transform(summedBlocks.begin(), summedBlocks.end(), averagedData.begin(),
      [&](std::pair<double3, double>& sample) {return (sample.first.x + sample.first.y + sample.first.z) / std::max(1.0, sample.second); });

    double res = 0;
    for (size_t i = 0; i < averagedData.size(); i++)
    {
      if (i == 0 || i == averagedData.size() - 1)
        res += averagedData[i];
      else if (i % 2 != 0)
        res += 4 * averagedData[i];
      else
        res += 2 * averagedData[i];
    }
    res = res * (delta / 3.0);
    return res;
    //double sum = 0.0;
    //for(size_t i = 1; i != averagedData.size(); ++i)
    //{
    //    sum += 0.5 * (averagedData[i] + averagedData[i-1]);
    //}
    //return sum / static_cast<double>(numberOfBins - 1);
  }

  std::pair<double, double> averageExcessChemicalPotentialDUdlambda() const
  {
    size_t numberOfSamples = numberOfBlocks;
    size_t degreesOfFreedom = numberOfSamples - 1;
    double average = averagedExcessChemicalPotentialDUdlambda();

    double sumOfSquares = 0.0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      double value = averagedExcessChemicalPotentialDUdlambda(blockIndex) - average;
      sumOfSquares += value * value;
    }
    double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
    double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
    double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
    double confidenceIntervalError = intermediateStandardNormalDeviate * standardError;

    return std::make_pair(average, confidenceIntervalError);
  }
  
  //====================================================================================================================

  std::pair<double, double> averageTotalChemicalPotential(double beta) const
  {
    double average = averagedExcessChemicalPotentialDUdlambda() + averagedIdealGasChemicalPotential(beta);

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingDensity[blockIndex].second / bookKeepingDensity[0].second > 0.5)
      {
        double value = (averagedExcessChemicalPotentialDUdlambda(blockIndex) + 
                        averagedIdealGasChemicalPotential(blockIndex, beta)) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }
    double confidenceIntervalError = 0.0;

    if (numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  std::pair<double, double> averageFugacityDUdlambda(double beta) const
  {
    double average = std::exp(beta * (averagedExcessChemicalPotentialDUdlambda() + 
                                      averagedIdealGasChemicalPotential(beta))) / beta;

    double sumOfSquares = 0.0;
    size_t numberOfSamples = 0;
    for (size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingDensity[blockIndex].second / bookKeepingDensity[0].second > 0.5)
      {
        double value = std::exp(beta * (averagedExcessChemicalPotentialDUdlambda(blockIndex) + 
                                        averagedIdealGasChemicalPotential(blockIndex, beta))) / beta - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }
    double confidenceIntervalError = 0.0;
    if (numberOfSamples >= 3)
    {
      size_t degreesOfFreedom = numberOfSamples - 1;
      double standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      double standardError = (1.0 / sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  friend Archive<std::ofstream> 
  &operator<<(Archive<std::ofstream> &archive, const PropertyLambdaProbabilityHistogram &p);

  friend Archive<std::ifstream> 
  &operator>>(Archive<std::ifstream> &archive, PropertyLambdaProbabilityHistogram &p);
};
