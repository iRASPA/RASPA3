module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <istream>
#include <numeric>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>
#endif

export module averages;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;

/**
 * \enum confidenceLevel
 * \brief Represents the different confidence levels for statistical calculations.
 */
export enum class confidenceLevel : int {
  percent_80 = 0,  ///< 80% confidence level.
  percent_90 = 1,  ///< 90% confidence level.
  percent_95 = 2,  ///< 95% confidence level.
  percent_98 = 3,  ///< 98% confidence level.
  percent_99 = 4   ///< 99% confidence level.
};

/**
 * \brief Table of standard normal deviates for different confidence levels and degrees of freedom.
 *
 * The table is indexed by degrees of freedom (number of bins - 1) and confidence level.
 * Values are sourced from the t-distribution critical values.
 * \see
 * https://sphweb.bumc.bu.edu/otlt/MPH-Modules/PH717-QuantCore/PH717-Module6-RandomError/PH717-Module6-RandomError11.html
 */
export constexpr inline std::array<std::array<double, 5>, 21> standardNormalDeviates{
    {{0.0, 0.0, 0.0, 0.0, 0.0},           {3.078, 6.314, 12.71, 31.82, 63.66}, {1.886, 2.920, 4.303, 6.965, 9.925},
     {1.638, 2.353, 3.182, 4.541, 5.841}, {1.533, 2.132, 2.776, 3.747, 4.604}, {1.476, 2.015, 2.571, 3.365, 4.032},
     {1.440, 1.943, 2.447, 3.143, 3.707}, {1.415, 1.895, 2.365, 2.998, 3.499}, {1.397, 1.860, 2.306, 2.896, 3.355},
     {1.383, 1.833, 2.262, 2.821, 3.250}, {1.372, 1.812, 2.228, 2.764, 3.169}, {1.362, 1.796, 2.201, 2.718, 3.106},
     {1.356, 1.782, 2.179, 2.681, 3.055}, {1.350, 1.771, 2.160, 2.650, 3.012}, {1.345, 1.761, 2.145, 2.624, 2.977},
     {1.341, 1.753, 2.131, 2.602, 2.947}, {1.337, 1.746, 2.120, 2.583, 2.921}, {1.333, 1.740, 2.110, 2.567, 2.898},
     {1.330, 1.734, 2.101, 2.552, 2.878}, {1.328, 1.729, 2.093, 2.539, 2.861}, {1.325, 1.725, 2.086, 2.528, 2.845}}};

/**
 * \brief Number of bins to use in block error estimation (must be less than 22).
 */
const int numberOfBins = 5;

/**
 * \brief Chosen confidence level for statistical calculations.
 */
export const int chosenConfidenceLevel = int(confidenceLevel::percent_95);

/**
 * \brief Standard normal deviate corresponding to the chosen confidence level and number of bins.
 */
constexpr double standardNormalDeviate = standardNormalDeviates[numberOfBins - 1][chosenConfidenceLevel];

/**
 * \brief Calculates the mean and confidence interval of a dataset.
 *
 * Computes the mean and the confidence interval based on the standard error and the standard normal deviate.
 *
 * \param data Reference to a vector containing the data samples.
 * \return A pair where the first element is the mean and the second is the confidence interval.
 */
export inline std::pair<double, double> meanConfidence(std::vector<double> &data)
{
  double size = static_cast<double>(data.size());
  double mean = std::accumulate(data.begin(), data.end(), 0.0) / size;
  double sumOfSquaresDiff = std::accumulate(data.begin(), data.end(), 0.0,
                                            [mean](double acc, double x) { return acc + (x - mean) * (x - mean); });

  double standardError = std::sqrt(sumOfSquaresDiff / (size * (size - 1.0)));

  return {mean, standardError * standardNormalDeviate};
}

/**
 * \brief Manages block error estimation for statistical analysis.
 *
 * The BlockErrorEstimation struct keeps track of the sampling process for block error analysis.
 * It divides the total number of samples into bins and updates the current bin as samples are processed.
 */
export struct BlockErrorEstimation
{
  std::size_t numberOfBins;          ///< Total number of bins.
  std::size_t currentSample{0};      ///< Current sample index.
  std::size_t numberOfSamples;       ///< Total number of samples.
  std::size_t currentBin{0};         ///< Current bin index.
  double binSize{};                  ///< Size of each bin.
  std::vector<std::size_t> nextBin;  ///< Indices of the next bin boundaries.

  /**
   * \brief Default constructor.
   *
   * Initializes a BlockErrorEstimation object with default values.
   */
  BlockErrorEstimation() {};

  /**
   * \brief Constructs a BlockErrorEstimation with specified parameters.
   *
   * Initializes the BlockErrorEstimation with the given number of bins and samples, and calculates the bin sizes.
   *
   * \param size Number of bins.
   * \param numberOfSamples Total number of samples.
   */
  BlockErrorEstimation(std::size_t size, std::size_t numberOfSamples)
      : numberOfBins(size), numberOfSamples(numberOfSamples), nextBin(size)
  {
    binSize = static_cast<double>(numberOfSamples) / static_cast<double>(numberOfBins);
    for (std::size_t i = 0; i != numberOfBins; ++i)
    {
      nextBin[i] = static_cast<std::size_t>(static_cast<double>(i + 1) * binSize);
    }
  }

  /**
   * \brief Updates the current sample and bin index.
   *
   * Sets the current sample index and updates the current bin index if the sample reaches the next bin boundary.
   *
   * \param sample The current sample index.
   */
  void setCurrentSample(std::size_t sample)
  {
    currentSample = sample;
    if (currentSample == nextBin[currentBin]) ++currentBin;
  }

  bool operator==(BlockErrorEstimation const &) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BlockErrorEstimation &blockerror);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BlockErrorEstimation &blockerror);
};
