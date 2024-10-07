module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <print>
#include <source_location>
#include <utility>
#include <vector>
#endif

module transition_matrix;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <cstddef>;
import <vector>;
import <array>;
import <map>;
import <functional>;
import <algorithm>;
import <utility>;
import <numeric>;
import <filesystem>;
import <fstream>;
import <iostream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <print>;
#endif

import archive;
import double3;

void TransitionMatrix::initialize()
{
  if (!doTMMC) return;

  cmatrix.resize(maxMacrostate - minMacrostate + 1);
  bias.resize(maxMacrostate - minMacrostate + 1);
  std::fill(bias.begin(), bias.end(), 1.0);
  lnpi.resize(maxMacrostate - minMacrostate + 1);
  forward_lnpi.resize(maxMacrostate - minMacrostate + 1);
  reverse_lnpi.resize(maxMacrostate - minMacrostate + 1);
  histogram.resize(maxMacrostate - minMacrostate + 1);
}

// C(No -> Nn) += p(o -> n)
// C(No -> No) += 1 − p(o -> n)
//
// translation: double3(0.0, 1.0, 0.0)
// insertion: double3(Pacc, 1.0 - Pacc, 0.0)
// insertion overlap-detected: double3(0.0, 1.0, 0.0)
// deletion: double3(0.0, 1.0 - Pacc, Pacc)
// deletion overlap-detected: double3(0.0, 1.0, 0.0)
void TransitionMatrix::updateMatrix(double3 Pacc, size_t oldN)
{
  if (!doTMMC) return;

  Pacc.clamp(0.0, 1.0);

  cmatrix[oldN - minMacrostate] += Pacc;
};

void TransitionMatrix::updateHistogram(size_t N)
{
  if (!doTMMC) return;

  if ((N > maxMacrostate) || (N < minMacrostate)) return;
  histogram[N - minMacrostate]++;
}

// return the biasing Factor
double TransitionMatrix::biasFactor(size_t newN, size_t oldN)
{
  if (!doTMMC || !useBias || !useTMBias) return 1.0;

  double TMMCBias = bias[newN - minMacrostate] - bias[oldN - minMacrostate];
  return std::exp(TMMCBias);
};

// From Vince Shen's pseudo code//
void TransitionMatrix::adjustBias()
{
  if (!doTMMC || !useBias || !useTMBias) return;

  if ((numberOfSteps % updateTMEvery != 0) || numberOfSteps == 0) return;

  numberOfUpdates++;

  // get the lowest and highest visited states in terms of loading
  size_t minVisitedN = static_cast<size_t>(std::distance(
      histogram.begin(), std::find_if(histogram.begin(), histogram.end(), [](const size_t &i) { return i; })));
  size_t maxVisitedN = static_cast<size_t>(
      std::distance(histogram.begin(),
                    std::find_if(histogram.rbegin(), histogram.rend(), [](const size_t &i) { return i; }).base()) -
      1);
  [[maybe_unused]] size_t nonzeroCount = maxVisitedN - minVisitedN + 1;

  lnpi[minVisitedN] = 0.0;
  double maxlnpi = lnpi[minVisitedN];
  // Update the lnpi for the sampled region//
  // x: -1; y: 0; z: +1//
  for (size_t i = minVisitedN; i < maxVisitedN; i++)
  {
    // Zhao's note: add protection to avoid numerical issues
    if (cmatrix[i].z != 0)
    {
      lnpi[i + 1] =
          lnpi[i] + std::log(cmatrix[i].z) - std::log(cmatrix[i].x + cmatrix[i].y + cmatrix[i].z);  // Forward//
    }
    forward_lnpi[i + 1] = lnpi[i + 1];
    if (cmatrix[i + 1].x != 0)
    {
      lnpi[i + 1] = lnpi[i + 1] - std::log(cmatrix[i + 1].x) +
                    std::log(cmatrix[i + 1].x + cmatrix[i + 1].y + cmatrix[i + 1].z);  // Reverse//
    }
    reverse_lnpi[i + 1] = lnpi[i + 1];
    if (lnpi[i + 1] > maxlnpi)
    {
      maxlnpi = lnpi[i + 1];
    }
  }

  // For the unsampled states, fill them with the minVisitedN/maxVisitedN stats
  for (size_t i = 0; i < minVisitedN; ++i)
  {
    lnpi[i] = lnpi[minVisitedN];
  }
  for (size_t i = maxVisitedN; i < maxMacrostate - minMacrostate + 1; ++i)
  {
    lnpi[i] = lnpi[maxVisitedN];
  }

  // Normalize
  for (size_t i = 0; i < maxMacrostate - minMacrostate + 1; ++i)
  {
    lnpi[i] -= maxlnpi;
  }
  double sumExps =
      std::accumulate(lnpi.begin(), lnpi.end(), 0.0, [](double sum, double item) { return sum + std::exp(item); });
  double normalFactor = -std::log(sumExps);

  for (size_t i = 0; i < maxMacrostate - minMacrostate + 1; ++i)
  {
    lnpi[i] += normalFactor;  // Zhao's note: mind the sign
    bias[i] = -lnpi[i];
  }
};

// Clear Collection matrix stats (used after initialization cycles)
void TransitionMatrix::clearCMatrix()
{
  if (!doTMMC || !rezeroAfterInitialization) return;

  numberOfSteps = 0;
  double3 temp = {0.0, 0.0, 0.0};
  std::fill(cmatrix.begin(), cmatrix.end(), temp);
  std::fill(histogram.begin(), histogram.end(), 0);
  std::fill(lnpi.begin(), lnpi.end(), 0.0);
  std::fill(bias.begin(), bias.end(), 1.0);
};

void TransitionMatrix::writeStatistics()
{
  std::ofstream textTMMCFile{};
  std::filesystem::path cwd = std::filesystem::current_path();

  std::string dirname = "TMMC/";
  // std::string dirname=std::print("TMMC/System_{}/", systemId);
  std::string fname = dirname + "/" + "TMMC_Statistics.txt";

  std::filesystem::path directoryName = cwd / dirname;
  std::filesystem::path fileName = cwd / fname;
  std::filesystem::create_directories(directoryName);
  textTMMCFile = std::ofstream(fileName, std::ios::out);

  if (doTMMC)
  {
    textTMMCFile << "Performed " << numberOfSteps << " Steps\n";
    textTMMCFile << "Collection Matrix Updated " << numberOfUpdates << " Times\n";
    textTMMCFile << "Min Macrostate : " << minMacrostate << '\n';
    textTMMCFile << "Max Macrostate : " << maxMacrostate << '\n';
    textTMMCFile << "N CM[-1] CM[0] CM[1] bias lnpi Forward_lnpi Reverse_lnpi histogram" << '\n';
    for (size_t j = minMacrostate; j < maxMacrostate + 1; j++)
    {
      size_t newj = j - minMacrostate;
      textTMMCFile << j << " " << cmatrix[newj].x << " " << cmatrix[newj].y << " " << cmatrix[newj].z << " "
                   << bias[newj] << " " << lnpi[newj] << " " << forward_lnpi[newj] << " " << reverse_lnpi[newj] << " "
                   << histogram[newj] << '\n';
    }
  }
};

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const TransitionMatrix &m)
{
  archive << m.versionNumber;

  archive << m.cmatrix;
  archive << m.bias;
  archive << m.lnpi;
  archive << m.forward_lnpi;
  archive << m.reverse_lnpi;
  archive << m.histogram;

  archive << m.numberOfSteps;
  archive << m.minMacrostate;
  archive << m.maxMacrostate;
  archive << m.updateTMEvery;
  archive << m.numberOfUpdates;

  archive << m.doTMMC;
  archive << m.useBias;
  archive << m.useTMBias;
  archive << m.rejectOutofBound;
  archive << m.rezeroAfterInitialization;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, TransitionMatrix &m)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > m.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'TransitionMatrix' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> m.cmatrix;
  archive >> m.bias;
  archive >> m.lnpi;
  archive >> m.forward_lnpi;
  archive >> m.reverse_lnpi;
  archive >> m.histogram;

  archive >> m.numberOfSteps;
  archive >> m.minMacrostate;
  archive >> m.maxMacrostate;
  archive >> m.updateTMEvery;
  archive >> m.numberOfUpdates;

  archive >> m.doTMMC;
  archive >> m.useBias;
  archive >> m.useTMBias;
  archive >> m.rejectOutofBound;
  archive >> m.rezeroAfterInitialization;

  return archive;
}
