module;

module transition_matrix;

import std;

import archive;
import double3;

void TransitionMatrix::initialize()
{
  if (!doTMMC) return;

  if (maxMacrostate < minMacrostate)
  {
    throw std::invalid_argument("TMMC maximum macrostate must not be smaller than its minimum macrostate");
  }

  const std::size_t numberOfMacrostates = maxMacrostate - minMacrostate + 1;
  cmatrix.assign(numberOfMacrostates, double3(0.0, 0.0, 0.0));
  bias.assign(numberOfMacrostates, 1.0);
  lnpi.assign(numberOfMacrostates, 0.0);
  forward_lnpi.assign(numberOfMacrostates, 0.0);
  reverse_lnpi.assign(numberOfMacrostates, 0.0);
  histogram.assign(numberOfMacrostates, 0uz);
}

// C(No -> Nn) += p(o -> n)
// C(No -> No) += 1 − p(o -> n)
//
// translation: double3(0.0, 1.0, 0.0)
// insertion: double3(0.0, 1.0 - Pacc, Pacc)
// insertion overlap-detected: double3(0.0, 1.0, 0.0)
// deletion: double3(Pacc, 1.0 - Pacc, 0.0)
// deletion overlap-detected: double3(0.0, 1.0, 0.0)
void TransitionMatrix::updateMatrix(double3 Pacc, std::size_t oldN)
{
  if (!doTMMC) return;

  Pacc.clamp(0.0, 1.0);

  if ((oldN < minMacrostate) || (oldN > maxMacrostate)) return;

  const std::size_t index = oldN - minMacrostate;
  if (index >= cmatrix.size()) return;

  cmatrix[index] += Pacc;
};

void TransitionMatrix::updateHistogram(std::size_t N)
{
  if (!doTMMC) return;

  if ((N > maxMacrostate) || (N < minMacrostate)) return;
  const std::size_t index = N - minMacrostate;
  if (index >= histogram.size()) return;
  histogram[index]++;
}

// return the biasing Factor
double TransitionMatrix::biasFactor(std::size_t newN, std::size_t oldN)
{
  if (!doTMMC) return 1.0;

  if ((newN < minMacrostate) || (newN > maxMacrostate) || (oldN < minMacrostate) || (oldN > maxMacrostate))
  {
    return rejectOutOfBound ? 0.0 : 1.0;
  }

  const std::size_t newIndex = newN - minMacrostate;
  const std::size_t oldIndex = oldN - minMacrostate;
  if ((newIndex >= bias.size()) || (oldIndex >= bias.size()))
  {
    return rejectOutOfBound ? 0.0 : 1.0;
  }

  if (!useBias || !useTMBias) return 1.0;

  double TMMCBias = bias[newIndex] - bias[oldIndex];
  return std::exp(TMMCBias);
};

// From Vince Shen's pseudo code//
void TransitionMatrix::adjustBias()
{
  if (!doTMMC || !useBias || !useTMBias) return;

  if ((numberOfSteps % updateTMEvery != 0) || numberOfSteps == 0) return;

  numberOfUpdates++;

  // get the lowest and highest visited states in terms of loading
  std::size_t minVisitedN = static_cast<std::size_t>(std::distance(
      histogram.begin(), std::find_if(histogram.begin(), histogram.end(), [](const std::size_t& i) { return i; })));
  std::size_t maxVisitedN = static_cast<std::size_t>(
      std::distance(histogram.begin(),
                    std::find_if(histogram.rbegin(), histogram.rend(), [](const std::size_t& i) { return i; }).base()) -
      1);
  [[maybe_unused]] std::size_t nonzeroCount = maxVisitedN - minVisitedN + 1;

  lnpi[minVisitedN] = 0.0;
  double maxlnpi = lnpi[minVisitedN];
  // Update the lnpi for the sampled region//
  // x: -1; y: 0; z: +1//
  for (std::size_t i = minVisitedN; i < maxVisitedN; i++)
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
  for (std::size_t i = 0; i < minVisitedN; ++i)
  {
    lnpi[i] = lnpi[minVisitedN];
  }
  for (std::size_t i = maxVisitedN; i < maxMacrostate - minMacrostate + 1; ++i)
  {
    lnpi[i] = lnpi[maxVisitedN];
  }

  // Normalize
  for (std::size_t i = 0; i < maxMacrostate - minMacrostate + 1; ++i)
  {
    lnpi[i] -= maxlnpi;
  }
  double sumExps =
      std::accumulate(lnpi.begin(), lnpi.end(), 0.0, [](double sum, double item) { return sum + std::exp(item); });
  double normalFactor = -std::log(sumExps);

  for (std::size_t i = 0; i < maxMacrostate - minMacrostate + 1; ++i)
  {
    lnpi[i] += normalFactor;  // Zhao's note: mind the sign
    bias[i] = -lnpi[i];
  }

  writeStatistics();
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

  std::filesystem::path fileName = cwd / statisticsFileName;
  std::filesystem::create_directories(fileName.parent_path());
  textTMMCFile = std::ofstream(fileName, std::ios::out);

  if (doTMMC)
  {
    std::print(textTMMCFile, "# performed: {} steps\n", numberOfSteps);
    std::print(textTMMCFile, "# collection matrix updated: {} times\n", numberOfUpdates);
    std::print(textTMMCFile, "# minimum microstate: {}\n", minMacrostate);
    std::print(textTMMCFile, "# maximum microstate: {}\n", maxMacrostate);
    std::print(textTMMCFile, "# column 1: N\n");
    std::print(textTMMCFile, "# column 2: CM[-1]\n");
    std::print(textTMMCFile, "# column 3: CM[ 0]0\n");
    std::print(textTMMCFile, "# column 4: CM[+1]\n");
    std::print(textTMMCFile, "# column 5: bias\n");
    std::print(textTMMCFile, "# column 6: lnpi\n");
    std::print(textTMMCFile, "# column 7: forward lnpi\n");
    std::print(textTMMCFile, "# column 8: reverse lnpi\n");
    std::print(textTMMCFile, "# column 9: histogram\n");
    std::print(textTMMCFile, "N CM[-1] CM[0] CM[1] bias lnpi Forward_lnpi Reverse_lnpi histogram\n");
    for (std::size_t j = minMacrostate; j < maxMacrostate + 1; j++)
    {
      std::size_t newj = j - minMacrostate;
      std::print(textTMMCFile, "{} {} {} {} {} {} {} {} {}\n", j, cmatrix[newj].x, cmatrix[newj].y, cmatrix[newj].z,
                 bias[newj], lnpi[newj], forward_lnpi[newj], reverse_lnpi[newj], histogram[newj]);
    }
  }
};

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const TransitionMatrix& m)
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
  archive << m.rejectOutOfBound;
  archive << m.rezeroAfterInitialization;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, TransitionMatrix& m)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > m.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
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
  archive >> m.rejectOutOfBound;
  archive >> m.rezeroAfterInitialization;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("TransitionMatrix: Error in binary restart\n"));
  }
#endif

  return archive;
}
