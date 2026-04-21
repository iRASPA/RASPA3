module;

export module property_loading;

import std;

import archive;
import int3;
import averages;
import loading_data;
import component;

inline std::pair<LoadingData, double> pair_sum(const std::pair<LoadingData, double> &lhs,
                                            const std::pair<LoadingData, double> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

export struct PropertyLoading
{
  PropertyLoading() {};

  PropertyLoading(std::size_t numberOfBlocks, std::size_t numberOfComponents)
      : numberOfBlocks(numberOfBlocks),
        numberOfComponents(numberOfComponents),
        bookKeepingLoadingData(numberOfBlocks)
  {
    // workaround for g++
    for(std::size_t i = 0; i < numberOfBlocks; ++i)
    {
      LoadingData loadings = LoadingData(numberOfComponents);
      bookKeepingLoadingData[i] = {loadings, 0.0};
    }
  }

  bool operator==(PropertyLoading const &) const = default;

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks;
  std::size_t numberOfComponents;
  std::vector<std::pair<LoadingData, double>> bookKeepingLoadingData;

  void resize(std::size_t newNumberOfComponents)
  {
    numberOfComponents = newNumberOfComponents;

    // workaround for g++
    bookKeepingLoadingData = std::vector<std::pair<LoadingData, double>>(numberOfBlocks);
    for(std::size_t i = 0; i < numberOfBlocks; ++i)
    {
      LoadingData loadings = LoadingData(numberOfComponents);
      bookKeepingLoadingData[i] = {loadings, 0.0};
    }
  }

  inline void addSample(std::size_t blockIndex, const LoadingData &loading, const double &weight)
  {
    bookKeepingLoadingData[blockIndex].first += weight * loading;
    bookKeepingLoadingData[blockIndex].second += weight;
  }

  //====================================================================================================================

  LoadingData averagedLoading(std::size_t blockIndex) const
  {
    return bookKeepingLoadingData[blockIndex].first / std::max(1.0, bookKeepingLoadingData[blockIndex].second);
  }

  LoadingData averagedLoading() const
  {
    // workaround for g++
    LoadingData loadings = LoadingData(numberOfComponents);

    std::pair<LoadingData, double> summedBlocks = std::make_pair(loadings, 0.0);

    for(const std::pair<LoadingData, double> &bookKeepingLoading : bookKeepingLoadingData)
    {
      summedBlocks.first += bookKeepingLoading.first;
      summedBlocks.second += bookKeepingLoading.second;
    }

    return summedBlocks.first / std::max(1.0, summedBlocks.second);;
  }

  std::pair<LoadingData, LoadingData> result() const
  {
    LoadingData average = averagedLoading();

    LoadingData sumOfSquares(numberOfComponents);
    std::size_t numberOfSamples = 0;
    for (std::size_t blockIndex = 0; blockIndex != numberOfBlocks; ++blockIndex)
    {
      if (bookKeepingLoadingData[blockIndex].second / std::max(1.0, bookKeepingLoadingData[0].second) > 0.5)
      {
        LoadingData value = averagedLoading(blockIndex) - average;
        sumOfSquares += value * value;
        ++numberOfSamples;
      }
    }
    LoadingData confidenceIntervalError(numberOfComponents);
    if (numberOfSamples >= 3)
    {
      std::size_t degreesOfFreedom = numberOfSamples - 1;
      LoadingData standardDeviation = sqrt((1.0 / static_cast<double>(degreesOfFreedom)) * sumOfSquares);
      LoadingData standardError = (1.0 / std::sqrt(static_cast<double>(numberOfSamples))) * standardDeviation;
      double intermediateStandardNormalDeviate = standardNormalDeviates[degreesOfFreedom][chosenConfidenceLevel];
      confidenceIntervalError = intermediateStandardNormalDeviate * standardError;
    }

    return std::make_pair(average, confidenceIntervalError);
  }

  std::pair<double, double> averageLoadingNumberOfMolecules(std::size_t comp) const;

  std::string writeAveragesStatistics(std::vector<Component> components, 
                                      std::optional<double> frameworkMass, 
                                      std::optional<int3> numberOfUnitCells) const;

  /**
   * \brief Returns a string representation of the PropertyLoading.
   *
   * Generates a simple string indicating a test representation of the PropertyLoading.
   *
   * \return A string representing the PropertyLoading.
   */
  std::string repr() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyLoading &l);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyLoading &l);
};
