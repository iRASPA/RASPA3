module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
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
import <array>;
import <vector>;
import <tuple>;
import <algorithm>;
import <numeric>;
import <cmath>;
import <iostream>;
import <istream>;
import <ostream>;
import <fstream>;
import <string>;
#endif

import archive;

export enum class confidenceLevel : int {
  percent_80 = 0,
  percent_90 = 1,
  percent_95 = 2,
  percent_98 = 3,
  percent_99 = 4
};

// https://sphweb.bumc.bu.edu/otlt/MPH-Modules/PH717-QuantCore/PH717-Module6-RandomError/PH717-Module6-RandomError11.html
export constexpr inline std::array<std::array<double, 5>, 21> standardNormalDeviates{
    {{0.0, 0.0, 0.0, 0.0, 0.0},           {3.078, 6.314, 12.71, 31.82, 63.66}, {1.886, 2.920, 4.303, 6.965, 9.925},
     {1.638, 2.353, 3.182, 4.541, 5.841}, {1.533, 2.132, 2.776, 3.747, 4.604}, {1.476, 2.015, 2.571, 3.365, 4.032},
     {1.440, 1.943, 2.447, 3.143, 3.707}, {1.415, 1.895, 2.365, 2.998, 3.499}, {1.397, 1.860, 2.306, 2.896, 3.355},
     {1.383, 1.833, 2.262, 2.821, 3.250}, {1.372, 1.812, 2.228, 2.764, 3.169}, {1.362, 1.796, 2.201, 2.718, 3.106},
     {1.356, 1.782, 2.179, 2.681, 3.055}, {1.350, 1.771, 2.160, 2.650, 3.012}, {1.345, 1.761, 2.145, 2.624, 2.977},
     {1.341, 1.753, 2.131, 2.602, 2.947}, {1.337, 1.746, 2.120, 2.583, 2.921}, {1.333, 1.740, 2.110, 2.567, 2.898},
     {1.330, 1.734, 2.101, 2.552, 2.878}, {1.328, 1.729, 2.093, 2.539, 2.861}, {1.325, 1.725, 2.086, 2.528, 2.845}}};

// The following two settings can be modified (keep numberOfBins smaller than 22)
const int numberOfBins = 5;
export const int chosenConfidenceLevel = int(confidenceLevel::percent_95);

constexpr double standardNormalDeviate = standardNormalDeviates[numberOfBins - 1][chosenConfidenceLevel];

export inline std::pair<double, double> meanConfidence(std::vector<double> &data)
{
  double size = static_cast<double>(data.size());
  double mean = std::accumulate(data.begin(), data.end(), 0.0) / size;
  double sumOfSquaresDiff = std::accumulate(data.begin(), data.end(), 0.0,
                                            [mean](double acc, double x) { return acc + (x - mean) * (x - mean); });

  double standardError = std::sqrt(sumOfSquaresDiff / (size * (size - 1.0)));

  return {mean, standardError * standardNormalDeviate};
}

export struct BlockErrorEstimation
{
  size_t numberOfBins;
  size_t currentSample{0};
  size_t numberOfSamples;
  size_t currentBin{0};
  double binSize{};
  std::vector<size_t> nextBin;

  BlockErrorEstimation() {};

  BlockErrorEstimation(size_t size, size_t numberOfSamples)
      : numberOfBins(size), numberOfSamples(numberOfSamples), nextBin(size)
  {
    binSize = static_cast<double>(numberOfSamples) / static_cast<double>(numberOfBins);
    for (size_t i = 0; i != numberOfBins; ++i)
    {
      nextBin[i] = static_cast<size_t>(static_cast<double>(i + 1) * binSize);
    }
  }

  bool operator==(BlockErrorEstimation const &) const = default;

  void setCurrentSample(size_t sample)
  {
    currentSample = sample;
    if (currentSample == nextBin[currentBin]) ++currentBin;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BlockErrorEstimation &blockerror);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BlockErrorEstimation &blockerror);
};
