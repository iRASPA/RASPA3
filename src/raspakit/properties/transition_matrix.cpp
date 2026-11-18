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
  // the bias is a logarithm (biasFactor is exp(bias[new] - bias[old]) and the rebuild sets bias = -ln Pi): it
  // starts at zero, like ln Pi
  bias.assign(numberOfChainStates, 0.0);
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
// insertion: double3(0.0, 1.0 - Pacc, Pacc)
// insertion overlap-detected: double3(0.0, 1.0, 0.0)
// deletion: double3(Pacc, 1.0 - Pacc, 0.0)
// deletion overlap-detected: double3(0.0, 1.0, 0.0)
// moves that leave N unchanged (translation, rotation, reinsertion, volume, Widom, ...): not recorded.
//   ln Pi only uses the off-diagonal entries normalized by the row totals; a stay on the diagonal would
//   rescale every row by the same move-mix factor, which cancels in the detailed-balance ratio.
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

  // A non-finite entry (an acceptance that overflowed, or inf * 0) would poison the row for the rest of the run:
  // std::clamp passes a NaN through, and std::isnan is compiled away under -ffast-math. A finite entry is clamped
  // (an acceptance above one leaves a negative stay); an overflowed acceptance is an accepted proposal; anything
  // else (a NaN, which fails every comparison) drops the record.
  auto finite = [](double value) { return value > -1.0e300 && value < 1.0e300; };
  if (!finite(Pacc.x) || !finite(Pacc.y) || !finite(Pacc.z))
  {
    auto overflowed = [](double value) { return value >= 1.0e300; };
    if (!overflowed(Pacc.x) && !overflowed(Pacc.z)) return;
    Pacc = double3(overflowed(Pacc.x) ? 1.0 : 0.0, 0.0, overflowed(Pacc.z) ? 1.0 : 0.0);
  }

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
  if (oldIndex >= bias.size())
  {
    // A walk that starts outside the chain (an initial loading past a bound, a restart from a run whose chain
    // reached (maxN, k > 0)) is let back in: a step toward the chain is taken unbiased, any other is refused.
    // Refusing every step left such a walk where it was for the whole run.
    if (!rejectOutOfBound) return 1.0;
    const std::size_t nLambda = lambdaBinCount();
    auto position = [&](std::size_t N, std::size_t lambdaBin) -> double
    {
      return (static_cast<double>(N) - static_cast<double>(minMacrostate)) * static_cast<double>(nLambda) +
             static_cast<double>(lambdaChain() ? lambdaBin : 0uz);
    };
    const double top = static_cast<double>((maxMacrostate - minMacrostate) * nLambda);
    const double oldPosition = position(oldN, oldLambdaBin);
    const double newPosition = position(newN, newLambdaBin);
    const bool towardChain =
        (oldPosition > top && newPosition < oldPosition) || (oldPosition < 0.0 && newPosition > oldPosition);
    return towardChain ? 1.0 : 0.0;
  }
  if (newIndex >= bias.size())
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

  // (an interval of zero, which a system put together in code can carry, is no periodic rebuild: the remainder by
  // it is a division by zero)
  if (updateTMEvery == 0uz || (numberOfSteps % updateTMEvery != 0) || numberOfSteps == 0) return;

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
  if (!doTMMC) return;

  // Without a transition-matrix bias to take over, Wang-Landau is frozen rather than retired: its bias keeps
  // flattening the walk, but its un-annealed running penalty no longer changes the sampled distribution through
  // production (it did, the switch returning before it touched Wang-Landau).
  if (!useBias || !useTMBias)
  {
    wangLandauFactor = 0.0;
    return;
  }
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

  // ln Pi by the detailed-balance recursion from the lowest visited state; forward_lnpi and reverse_lnpi
  // are the running sums of the forward and reverse transition log-probabilities along the same walk
  // (diagnostic columns of the statistics file); the forward sum takes the normalization shift below, so
  // that lnpi = forward_lnpi - reverse_lnpi exactly. All three start at zero at the lowest visited state: a value left
  // there from an earlier rebuild with another visited range would otherwise be carried into every state above it.
  lnpi[minVisitedN] = 0.0;
  forward_lnpi[minVisitedN] = 0.0;
  reverse_lnpi[minVisitedN] = 0.0;
  double maxlnpi = lnpi[minVisitedN];
  // Update the lnpi for the sampled region//
  // x: -1; y: 0; z: +1//
  for (std::size_t i = minVisitedN; i < maxVisitedN; i++)
  {
    const double forwardTotal = cmatrix[i].x + cmatrix[i].y + cmatrix[i].z;
    const double reverseTotal = cmatrix[i + 1].x + cmatrix[i + 1].y + cmatrix[i + 1].z;
    if (linkSupported(i))
    {
      const double forward = std::log(cmatrix[i].z) - std::log(forwardTotal);
      const double reverse = std::log(cmatrix[i + 1].x) - std::log(reverseTotal);
      forward_lnpi[i + 1] = forward_lnpi[i] + forward;
      reverse_lnpi[i + 1] = reverse_lnpi[i] + reverse;
      lnpi[i + 1] = lnpi[i] + forward - reverse;
    }
    else
    {
      // A one-sided link does not identify a detailed-balance ratio. Keep the live bias flat
      // across it until both directions have support; applying only the reverse term invents a
      // free-energy step and can create a wall before the walker has sampled the forward hop.
      // The plateau is a bias, not an estimate: the statistics file numbers the supported runs ('segment') and
      // the reweighted loading keeps to one of them (supportedMoleculeRange).
      lnpi[i + 1] = lnpi[i];
      forward_lnpi[i + 1] = forward_lnpi[i];
      reverse_lnpi[i + 1] = reverse_lnpi[i];
    }
    if (lnpi[i + 1] > maxlnpi)
    {
      maxlnpi = lnpi[i + 1];
    }
  }

  // For the unsampled states, fill them with the minVisitedN/maxVisitedN stats
  for (std::size_t i = 0; i < minVisitedN; ++i)
  {
    lnpi[i] = lnpi[minVisitedN];
    forward_lnpi[i] = forward_lnpi[minVisitedN];
    reverse_lnpi[i] = reverse_lnpi[minVisitedN];
  }
  for (std::size_t i = maxVisitedN; i < lnpi.size(); ++i)
  {
    lnpi[i] = lnpi[maxVisitedN];
    forward_lnpi[i] = forward_lnpi[maxVisitedN];
    reverse_lnpi[i] = reverse_lnpi[maxVisitedN];
  }

  // Normalize
  for (std::size_t i = 0; i < lnpi.size(); ++i)
  {
    lnpi[i] -= maxlnpi;
  }
  double sumExps =
      std::accumulate(lnpi.begin(), lnpi.end(), 0.0, [](double sum, double item) { return sum + std::exp(item); });
  double normalFactor = -std::log(sumExps);

  // the forward sum carries the same shift, so that lnpi = forward_lnpi - reverse_lnpi holds exactly in the
  // statistics file (the reverse sum is the raw running sum, zero at the lowest visited state)
  for (std::size_t i = 0; i < lnpi.size(); ++i)
  {
    forward_lnpi[i] += normalFactor - maxlnpi;
  }
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
  numberOfUpdates = 0;
  std::fill(cmatrix.begin(), cmatrix.end(), double3(0.0, 0.0, 0.0));
  std::fill(histogram.begin(), histogram.end(), 0uz);
  // ln Pi belongs to the statistics that were dropped (the bias is kept on its own): left in place it was written
  // next to an empty collection matrix until the first rebuild, and for good when production visited nothing
  std::fill(lnpi.begin(), lnpi.end(), 0.0);
  std::fill(forward_lnpi.begin(), forward_lnpi.end(), 0.0);
  std::fill(reverse_lnpi.begin(), reverse_lnpi.end(), 0.0);
};

// Clear Collection matrix stats (used after initialization cycles)
void TransitionMatrix::clearCMatrix()
{
  if (!doTMMC || !rezeroAfterInitialization) return;

  numberOfSteps = 0;
  numberOfUpdates = 0;
  double3 temp = {0.0, 0.0, 0.0};
  std::fill(cmatrix.begin(), cmatrix.end(), temp);
  std::fill(histogram.begin(), histogram.end(), 0.0);
  std::fill(lnpi.begin(), lnpi.end(), 0.0);
  std::fill(forward_lnpi.begin(), forward_lnpi.end(), 0.0);
  std::fill(reverse_lnpi.begin(), reverse_lnpi.end(), 0.0);
  std::fill(bias.begin(), bias.end(), 0.0);
  std::fill(wangLandauHistogram.begin(), wangLandauHistogram.end(), 0.0);
  wangLandauVisits = 0uz;
  // the bias Wang-Landau built is gone, so its modification factor starts over: left annealed, a cleared bias could
  // not be rebuilt and the walk could not re-explore the macrostates
  wangLandauFactor = TransitionMatrix{}.wangLandauFactor;
};

std::optional<std::pair<std::size_t, std::size_t>> TransitionMatrix::visitedMoleculeRange() const
{
  if (!doTMMC || histogram.size() != numberOfChainStates()) return std::nullopt;
  const std::size_t nLambda = lambdaBinCount();
  std::optional<std::pair<std::size_t, std::size_t>> range{};
  for (std::size_t N = minMacrostate; N <= maxMacrostate; ++N)
  {
    std::size_t visits = 0uz;
    for (std::size_t k = 0; k < nLambda; ++k)
    {
      const std::size_t index = chainIndex(N, k);
      if (index < histogram.size()) visits += histogram[index];
    }
    if (visits == 0uz) continue;
    if (!range.has_value()) range = std::make_pair(N, N);
    range->second = N;
  }
  return range;
}

std::optional<std::pair<std::size_t, std::size_t>> TransitionMatrix::supportedMoleculeRange() const
{
  const std::optional<std::pair<std::size_t, std::size_t>> visited = visitedMoleculeRange();
  if (!visited.has_value() || cmatrix.size() != numberOfChainStates()) return visited;
  const std::size_t nLambda = lambdaBinCount();
  // the step N -> N + 1 is the walk from (N, 0) to (N + 1, 0): one link, or the nLambda links of an (N, lambda) chain
  auto stepSupported = [&](std::size_t N) -> bool
  {
    const std::size_t first = chainIndex(N, 0uz);
    for (std::size_t k = 0; k < nLambda; ++k)
    {
      if (!linkSupported(first + k)) return false;
    }
    return true;
  };
  std::pair<std::size_t, std::size_t> longest{visited->first, visited->first};
  std::size_t runStart = visited->first;
  for (std::size_t N = visited->first; N <= visited->second; ++N)
  {
    if (N - runStart > longest.second - longest.first) longest = {runStart, N};
    if (N < visited->second && !stepSupported(N)) runStart = N + 1uz;
  }
  return longest;
}

std::optional<std::pair<double, double>> TransitionMatrix::reweightedMoleculeCountMeanAndVariance() const
{
  if (!doTMMC || lnpi.size() != numberOfChainStates() || histogram.size() != numberOfChainStates()) return std::nullopt;
  // ln Pi is one estimate only inside a run of supported links: across an unsupported one it is carried flat
  const std::optional<std::pair<std::size_t, std::size_t>> range = supportedMoleculeRange();
  if (!range.has_value()) return std::nullopt;
  double largest = -1e300;
  for (std::size_t N = range->first; N <= range->second; ++N) largest = std::max(largest, lnpi[chainIndex(N, 0uz)]);
  double sum = 0.0;
  double sumN = 0.0;
  double sumNN = 0.0;
  for (std::size_t N = range->first; N <= range->second; ++N)
  {
    const double weight = std::exp(lnpi[chainIndex(N, 0uz)] - largest);
    const double n = static_cast<double>(N);
    sum += weight;
    sumN += weight * n;
    sumNN += weight * n * n;
  }
  if (!(sum > 0.0)) return std::nullopt;
  const double mean = sumN / sum;
  return std::make_pair(mean, std::max(0.0, sumNN / sum - mean * mean));
}

void TransitionMatrix::writeStatistics()
{
  // the supported runs of the chain: states joined by supported links share a number, counted from the lowest
  // visited state; -1 outside the visited range. Differences of ln Pi mean something only inside one run.
  std::vector<std::int64_t> segment(cmatrix.size(), -1);
  {
    const auto visited = [](const std::size_t& visits) { return visits != 0uz; };
    const auto first = std::find_if(histogram.begin(), histogram.end(), visited);
    if (first != histogram.end() && histogram.size() == cmatrix.size())
    {
      const std::size_t lowest = static_cast<std::size_t>(std::distance(histogram.begin(), first));
      const std::size_t highest = static_cast<std::size_t>(
          std::distance(histogram.begin(), std::find_if(histogram.rbegin(), histogram.rend(), visited).base()) - 1);
      std::int64_t current = 0;
      for (std::size_t index = lowest; index <= highest; ++index)
      {
        segment[index] = current;
        if (index < highest && !linkSupported(index)) ++current;
      }
    }
  }

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
      std::print(textTMMCFile, "# column 10: segment (states joined by links with accepted transitions in both\n");
      std::print(textTMMCFile, "#            directions; lnpi is carried flat across an unsupported link, and its\n");
      std::print(textTMMCFile, "#            differences are estimates only inside one segment; -1: not visited)\n");
      std::print(textTMMCFile, "N CM[-1] CM[0] CM[1] bias lnpi Forward_lnpi Reverse_lnpi histogram segment\n");
      for (std::size_t j = minMacrostate; j < maxMacrostate + 1; j++)
      {
        std::size_t newj = j - minMacrostate;
        std::print(textTMMCFile, "{} {} {} {} {} {} {} {} {} {}\n", j, cmatrix[newj].x, cmatrix[newj].y, cmatrix[newj].z,
                   bias[newj], lnpi[newj], forward_lnpi[newj], reverse_lnpi[newj], histogram[newj], segment[newj]);
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
      std::print(textTMMCFile, "# column 11: segment (states joined by links with accepted transitions in both\n");
      std::print(textTMMCFile, "#            directions; lnpi is carried flat across an unsupported link, and its\n");
      std::print(textTMMCFile, "#            differences are estimates only inside one segment; -1: not visited)\n");
      std::print(textTMMCFile, "# the chain ends at (maximum, 0): the window bounds are integer states\n");
      std::print(textTMMCFile, "N k CM[-1] CM[0] CM[1] bias lnpi Forward_lnpi Reverse_lnpi histogram segment\n");
      const std::size_t nLambda = lambdaBinCount();
      for (std::size_t j = minMacrostate; j < maxMacrostate + 1; j++)
      {
        for (std::size_t k = 0; k < nLambda; ++k)
        {
          const std::size_t index = chainIndex(j, k);
          if (index >= cmatrix.size()) continue;
          std::print(textTMMCFile, "{} {} {} {} {} {} {} {} {} {} {}\n", j, k, cmatrix[index].x, cmatrix[index].y,
                     cmatrix[index].z, bias[index], lnpi[index], forward_lnpi[index], reverse_lnpi[index],
                     histogram[index], segment[index]);
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
  archive << m.currentLambdaBin;

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

  archive >> m.wangLandauHistogram;
  archive >> m.wangLandauFactor;
  archive >> m.wangLandauFactorFloor;
  archive >> m.wangLandauFlatness;
  archive >> m.wangLandauCheckEvery;
  archive >> m.wangLandauVisits;
  archive >> m.useWangLandau;

  archive >> m.statisticsFileName;
  archive >> m.numberOfLambdaBins;
  archive >> m.currentLambdaBin;

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
