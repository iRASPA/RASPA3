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

  const std::size_t numberOfChainStates = this->numberOfChainStates();
  cmatrix.assign(numberOfChainStates, double3(0.0, 0.0, 0.0));
  bias.assign(numberOfChainStates, 1.0);
  lnpi.assign(numberOfChainStates, 0.0);
  forward_lnpi.assign(numberOfChainStates, 0.0);
  reverse_lnpi.assign(numberOfChainStates, 0.0);
  histogram.assign(numberOfChainStates, 0uz);
  wangLandauHistogram.assign(numberOfChainStates, 0.0);
  wangLandauVisits = 0uz;
  if (!lambdaChain()) currentLambdaBin = 0uz;
}

// C(No -> Nn) += p(o -> n)
// C(No -> No) += 1 − p(o -> n)
//
// translation: double3(0.0, 1.0, 0.0)
// insertion: double3(0.0, 1.0 - Pacc, Pacc)
// insertion overlap-detected: double3(0.0, 1.0, 0.0)
// deletion: double3(Pacc, 1.0 - Pacc, 0.0)
// deletion overlap-detected: double3(0.0, 1.0, 0.0)
// 1D-N CFCMC lambda hop: double3(0.0, 0.0, 0.0) — skip (integer N is unchanged).
// (N, λ) chain: a lambda hop is a real ±1 step and is recorded like insertion/deletion.
void TransitionMatrix::updateMatrix(double3 Pacc, std::size_t oldN)
{
  updateMatrix(Pacc, oldN, currentLambdaBin);
}

void TransitionMatrix::updateMatrix(double3 Pacc, std::size_t oldN, std::size_t oldLambdaBin)
{
  if (!doTMMC) return;

  if (Pacc.x == 0.0 && Pacc.y == 0.0 && Pacc.z == 0.0) return;

  Pacc.clamp(0.0, 1.0);

  const std::size_t index = chainIndex(oldN, oldLambdaBin);
  if (index >= cmatrix.size()) return;

  cmatrix[index] += Pacc;
};

void TransitionMatrix::updateHistogram(std::size_t N, std::size_t lambdaBin)
{
  if (!doTMMC) return;

  const std::size_t index = chainIndex(N, lambdaBin);
  if (index >= histogram.size()) return;
  histogram[index]++;
}

// Wang-Landau: penalise the state that is occupied so the walk is pushed towards the ones that are not.
void TransitionMatrix::visitWangLandau(std::size_t N, std::size_t lambdaBin)
{
  if (!doTMMC || !useBias || !useWangLandau) return;

  const std::size_t index = chainIndex(N, lambdaBin);
  if ((index >= bias.size()) || (index >= wangLandauHistogram.size())) return;

  bias[index] -= wangLandauFactor;
  wangLandauHistogram[index] += 1.0;

  if (++wangLandauVisits < wangLandauCheckEvery) return;
  wangLandauVisits = 0uz;

  // Only differences of the bias are ever used, so the running penalty is free to be shifted; holding the
  // largest entry at zero keeps it from drifting somewhere exp() of a difference cannot be represented.
  const double largest = *std::max_element(bias.begin(), bias.end());
  for (double &value : bias) value -= largest;

  if (wangLandauFactor <= wangLandauFactorFloor) return;

  double smallest = wangLandauHistogram.front();
  double total = 0.0;
  for (double visits : wangLandauHistogram)
  {
    smallest = std::min(smallest, visits);
    total += visits;
  }
  const double mean = total / static_cast<double>(wangLandauHistogram.size());
  if (!(mean > 0.0) || smallest < wangLandauFlatness * mean) return;

  // Flat over the whole window: refine the bias and start the next stage.
  wangLandauFactor = std::max(0.5 * wangLandauFactor, wangLandauFactorFloor);
  std::fill(wangLandauHistogram.begin(), wangLandauHistogram.end(), 0.0);
}

// return the biasing Factor
double TransitionMatrix::biasFactor(std::size_t newN, std::size_t oldN)
{
  if (!lambdaChain()) return biasFactor(newN, oldN, 0uz, 0uz);
  if (newN == oldN) return biasFactor(newN, oldN, currentLambdaBin, currentLambdaBin);
  if (newN > oldN) return biasFactor(newN, oldN, 0uz, lastLambdaBin());
  return biasFactor(newN, oldN, lastLambdaBin(), 0uz);
}

double TransitionMatrix::biasFactor(std::size_t newN, std::size_t oldN, std::size_t newLambdaBin,
                                   std::size_t oldLambdaBin)
{
  if (!doTMMC) return 1.0;

  const std::size_t newIndex = chainIndex(newN, newLambdaBin);
  const std::size_t oldIndex = chainIndex(oldN, oldLambdaBin);
  if ((newIndex >= bias.size()) || (oldIndex >= bias.size()))
  {
    return rejectOutOfBound ? 0.0 : 1.0;
  }

  if (!useBias || (!useTMBias && !useWangLandau)) return 1.0;

  // Clamped because a Wang-Landau bias is a running sum: early on, before the histogram has flattened, two
  // states can be far enough apart that exp() of the difference is not finite, and an infinity multiplied
  // by a zero acceptance probability is a NaN that silently rejects.
  double TMMCBias = std::clamp(bias[newIndex] - bias[oldIndex], -500.0, 500.0);
  return std::exp(TMMCBias);
};

// From Vince Shen's pseudo code//
void TransitionMatrix::adjustBias()
{
  if (!doTMMC || !useBias || !useTMBias) return;

  if ((numberOfSteps % updateTMEvery != 0) || numberOfSteps == 0) return;

  recomputeLnPiAndBias();
};

// Errington/Shen hybrid, second half: Wang-Landau explores the window, the transition-matrix bias
// flattens it. WL's per-visit penalty is what pushes a walker over ln Pi spans of several hundred
// into states the collection matrix has never seen, but its modification factor practically never
// anneals here (a window pinned at a filling fugacity keeps a lopsided histogram, so the flatness
// test keeps failing at factor 1.0) and a factor-1 bias is a sawtooth: the walker rattles around
// its equilibrium wall and the rest of the window collects tens of visits. Once the window has
// been crossed, -ln Pi from the collection matrix IS the exact flattening bias; sampling under it
// gives every macrostate comparable statistics, which is what the ln Pi increments need.
void TransitionMatrix::switchToTMBias()
{
  if (!doTMMC || !useBias || !useTMBias) return;

  useWangLandau = false;
  recomputeLnPiAndBias();
};

void TransitionMatrix::recomputeLnPiAndBias()
{
  numberOfUpdates++;

  // nothing visited yet: no statistics to derive a bias from
  if (std::find_if(histogram.begin(), histogram.end(), [](const std::size_t& i) { return i; }) == histogram.end())
  {
    return;
  }

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
    const double forwardTotal = cmatrix[i].x + cmatrix[i].y + cmatrix[i].z;
    const double reverseTotal = cmatrix[i + 1].x + cmatrix[i + 1].y + cmatrix[i + 1].z;
    const bool supportedLink =
        cmatrix[i].z > 0.0 && cmatrix[i + 1].x > 0.0 && forwardTotal > 0.0 && reverseTotal > 0.0;
    if (supportedLink)
    {
      const double forward = std::log(cmatrix[i].z) - std::log(forwardTotal);
      const double reverse = std::log(cmatrix[i + 1].x) - std::log(reverseTotal);
      forward_lnpi[i + 1] = lnpi[i] + forward;
      lnpi[i + 1] = lnpi[i] + forward - reverse;
    }
    else
    {
      // A one-sided link does not identify a detailed-balance ratio. Keep the live bias flat
      // across it until both directions have support; applying only the reverse term invents a
      // free-energy step and can create a wall before the walker has sampled the forward hop.
      lnpi[i + 1] = lnpi[i];
      forward_lnpi[i + 1] = lnpi[i];
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
  for (std::size_t i = maxVisitedN; i < lnpi.size(); ++i)
  {
    lnpi[i] = lnpi[maxVisitedN];
  }

  // Normalize
  for (std::size_t i = 0; i < lnpi.size(); ++i)
  {
    lnpi[i] -= maxlnpi;
  }
  double sumExps =
      std::accumulate(lnpi.begin(), lnpi.end(), 0.0, [](double sum, double item) { return sum + std::exp(item); });
  double normalFactor = -std::log(sumExps);

  for (std::size_t i = 0; i < lnpi.size(); ++i)
  {
    lnpi[i] += normalFactor;  // Zhao's note: mind the sign
    // Wang-Landau owns the bias when it is running. Overwriting it here would undo the exploration it has
    // paid for: outside the sampled range this ln pi is flat by construction (the loops above fill it with
    // the value at the edge), so the walker would lose the very gradient that carries it into the states
    // the collection matrix has not reached yet.
    if (!useWangLandau) bias[i] = -lnpi[i];
  }

  writeStatistics();
};

// Reset the collection matrix and the visit histogram but keep the bias (including the
// Wang-Landau state). Used at production start: the pre-production statistics were taken on
// configurations that had not yet relaxed into the deep adsorption sites, and because a row
// estimate is a mean of acceptance probabilities with a heavy upper tail, those early
// overestimates would dominate the ln Pi readout no matter how long production runs.
void TransitionMatrix::clearStatisticsKeepBias()
{
  if (!doTMMC) return;

  numberOfSteps = 0;
  std::fill(cmatrix.begin(), cmatrix.end(), double3(0.0, 0.0, 0.0));
  std::fill(histogram.begin(), histogram.end(), 0uz);
};

// Clear Collection matrix stats (used after initialization cycles)
void TransitionMatrix::clearCMatrix()
{
  if (!doTMMC || !rezeroAfterInitialization) return;

  numberOfSteps = 0;
  double3 temp = {0.0, 0.0, 0.0};
  std::fill(cmatrix.begin(), cmatrix.end(), temp);
  std::fill(histogram.begin(), histogram.end(), 0.0);
  std::fill(lnpi.begin(), lnpi.end(), 0.0);
  std::fill(bias.begin(), bias.end(), 1.0);
  std::fill(wangLandauHistogram.begin(), wangLandauHistogram.end(), 0.0);
  wangLandauVisits = 0uz;
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
    if (useWangLandau)
    {
      std::print(textTMMCFile, "# Wang-Landau modification factor: {} (floor {})\n", wangLandauFactor,
                 wangLandauFactorFloor);
    }
    std::print(textTMMCFile, "# minimum microstate: {}\n", minMacrostate);
    std::print(textTMMCFile, "# maximum microstate: {}\n", maxMacrostate);
    std::print(textTMMCFile, "# lambda bins: {}\n", lambdaBinCount());
    if (!lambdaChain())
    {
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
    else
    {
      std::print(textTMMCFile, "# column 1: N\n");
      std::print(textTMMCFile, "# column 2: lambda bin\n");
      std::print(textTMMCFile, "# column 3: CM[-1]\n");
      std::print(textTMMCFile, "# column 4: CM[ 0]\n");
      std::print(textTMMCFile, "# column 5: CM[+1]\n");
      std::print(textTMMCFile, "# column 6: bias\n");
      std::print(textTMMCFile, "# column 7: lnpi\n");
      std::print(textTMMCFile, "# column 8: forward lnpi\n");
      std::print(textTMMCFile, "# column 9: reverse lnpi\n");
      std::print(textTMMCFile, "# column 10: histogram\n");
      std::print(textTMMCFile, "N k CM[-1] CM[0] CM[1] bias lnpi Forward_lnpi Reverse_lnpi histogram\n");
      const std::size_t nLambda = lambdaBinCount();
      for (std::size_t j = minMacrostate; j < maxMacrostate + 1; j++)
      {
        for (std::size_t k = 0; k < nLambda; ++k)
        {
          const std::size_t index = chainIndex(j, k);
          if (index >= cmatrix.size()) continue;
          std::print(textTMMCFile, "{} {} {} {} {} {} {} {} {} {}\n", j, k, cmatrix[index].x, cmatrix[index].y,
                     cmatrix[index].z, bias[index], lnpi[index], forward_lnpi[index], reverse_lnpi[index],
                     histogram[index]);
        }
      }
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

  archive << m.wangLandauHistogram;
  archive << m.wangLandauFactor;
  archive << m.wangLandauFactorFloor;
  archive << m.wangLandauFlatness;
  archive << m.wangLandauCheckEvery;
  archive << m.wangLandauVisits;
  archive << m.useWangLandau;

  archive << m.statisticsFileName;
  archive << m.numberOfLambdaBins;

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

  if (versionNumber >= 2)
  {
    archive >> m.wangLandauHistogram;
    archive >> m.wangLandauFactor;
    archive >> m.wangLandauFactorFloor;
    archive >> m.wangLandauFlatness;
    archive >> m.wangLandauCheckEvery;
    archive >> m.wangLandauVisits;
    archive >> m.useWangLandau;
  }

  archive >> m.statisticsFileName;
  if (versionNumber >= 3)
  {
    archive >> m.numberOfLambdaBins;
  }
  else
  {
    m.numberOfLambdaBins = 1uz;
  }

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
