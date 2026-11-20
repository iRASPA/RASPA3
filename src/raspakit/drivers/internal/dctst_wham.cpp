module;

module dctst_wham;

import std;

DcTSTWham::Result DcTSTWham::solve(const std::vector<std::vector<double>>& counts,
                                   const std::vector<std::vector<double>>& reducedPotential, double tolerance,
                                   std::size_t maximumIterations, const std::vector<double>& initialFreeEnergies)
{
  Result result;
  const std::size_t numberOfStates = counts.size();
  if (numberOfStates == 0uz) return result;
  const std::size_t numberOfBins = counts.front().size();

  // A bin is usable when at least one state with samples can reach it (its reduced potential there is not the
  // forbidden ~1e30). Samples in a bin that every sampled state forbids can not be weighted: they are left out
  // of ln M_b and of ln n_k alike. Leaving a state without samples can strand further bins, hence the loop.
  // The samples a state has in a bin it forbids itself (a pilot sample just past its dividing surface, say) are
  // left out in the same way, whatever the other states do there: the state is excluded from the mixture of that
  // bin, so that counting them in n_k inflated its weight in every other bin, and counting them in M_b put
  // samples there that no state of the mixture accounts for.
  auto allows = [&](std::size_t k, std::size_t b) { return -reducedPotential[k][b] > forbiddenLogWeight; };
  std::vector<double> stateTotals(numberOfStates, 0.0);
  std::vector<bool> binUsable(numberOfBins, true);
  for (bool changed = true; changed;)
  {
    changed = false;
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      stateTotals[k] = 0.0;
      for (std::size_t b = 0; b < numberOfBins; ++b)
      {
        if (binUsable[b] && allows(k, b)) stateTotals[k] += counts[k][b];
      }
    }
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      if (!binUsable[b]) continue;
      bool reachable = false;
      for (std::size_t k = 0; k < numberOfStates && !reachable; ++k)
      {
        reachable = stateTotals[k] > 0.0 && allows(k, b);
      }
      if (!reachable)
      {
        binUsable[b] = false;
        changed = true;
      }
    }
  }

  // ln n_k and ln M_b
  std::vector<double> logStateCounts(numberOfStates, logZero);
  std::vector<double> logBinCounts(numberOfBins, logZero);
  for (std::size_t k = 0; k < numberOfStates; ++k)
  {
    if (stateTotals[k] > 0.0) logStateCounts[k] = std::log(stateTotals[k]);
  }
  for (std::size_t b = 0; b < numberOfBins; ++b)
  {
    if (!binUsable[b]) continue;
    double total{0.0};
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      if (allows(k, b)) total += counts[k][b];
    }
    if (total > 0.0) logBinCounts[b] = std::log(total);
  }

  result.freeEnergies.assign(numberOfStates, 0.0);
  if (initialFreeEnergies.size() == numberOfStates)
  {
    result.freeEnergies = initialFreeEnergies;
  }
  result.logDensity.assign(numberOfBins, logZero);
  std::vector<double> updated(numberOfStates, 0.0);
  std::vector<bool> takesPart(numberOfStates, true);

  // The gauge is fixed on the first state with samples. A state without samples does not enter Omega, and its
  // g can be anything (a forbidden bias over every sampled bin puts it at ~ +1e30): fixed on such a state the
  // gauge moved every other g by that value and the iteration alternated without converging.
  std::size_t referenceState = 0uz;
  while (referenceState < numberOfStates && logStateCounts[referenceState] <= logZeroThreshold) ++referenceState;
  if (referenceState == numberOfStates)
  {
    result.freeEnergies.assign(numberOfStates, 0.0);
    return result;
  }

  // Omega_b from the current g_k
  auto updateDensity = [&]()
  {
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      if (logBinCounts[b] <= logZeroThreshold)
      {
        result.logDensity[b] = logZero;
        continue;
      }
      double largest = logZero;
      for (std::size_t k = 0; k < numberOfStates; ++k)
      {
        if (logStateCounts[k] <= logZeroThreshold) continue;
        largest = std::max(largest, logStateCounts[k] + result.freeEnergies[k] - reducedPotential[k][b]);
      }
      // A bin no sampled state can reach: every one of them carries a forbidden bias there (a reduced potential
      // of 1e30, see the callers), so that 'largest' is about -1e30. That is far above the 'log of zero', and
      // the bin got the density ln M_b + 1e30: it then entered every g_k as the difference of two numbers of
      // 1e30, which double precision resolves to 1e14 at best. Its samples can not be weighted and are left out.
      if (largest <= forbiddenLogWeight)
      {
        result.logDensity[b] = logZero;
        continue;
      }
      double sum{0.0};
      for (std::size_t k = 0; k < numberOfStates; ++k)
      {
        if (logStateCounts[k] <= logZeroThreshold) continue;
        sum += std::exp(logStateCounts[k] + result.freeEnergies[k] - reducedPotential[k][b] - largest);
      }
      if (!(sum > 0.0))
      {
        result.logDensity[b] = logZero;
        continue;
      }
      result.logDensity[b] = logBinCounts[b] - largest - std::log(sum);
    }

  };

  double residual = 1e30;
  std::size_t iteration = 0uz;
  for (; iteration < maximumIterations && residual > tolerance; ++iteration)
  {
    updateDensity();

    // g_k
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      double largest = logZero;
      for (std::size_t b = 0; b < numberOfBins; ++b)
      {
        if (result.logDensity[b] <= logZeroThreshold) continue;
        largest = std::max(largest, result.logDensity[b] - reducedPotential[k][b]);
      }
      // A state that forbids every bin with a density (a pilot without samples in occupancy mode) has a largest
      // term of about -1e30, far above the 'log of zero': its g became +1e30 instead of the 0 of a state that
      // takes no part.
      takesPart[k] = largest > forbiddenLogWeight;
      if (!takesPart[k])
      {
        updated[k] = 0.0;
        continue;
      }
      double sum{0.0};
      for (std::size_t b = 0; b < numberOfBins; ++b)
      {
        if (result.logDensity[b] <= logZeroThreshold) continue;
        sum += std::exp(result.logDensity[b] - reducedPotential[k][b] - largest);
      }
      updated[k] = -(largest + std::log(sum));
    }

    const double gauge = updated[referenceState];
    residual = 0.0;
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      if (takesPart[k]) updated[k] -= gauge;  // the others stay at zero
      if (logStateCounts[k] <= logZeroThreshold) continue;  // an unsampled state does not take part
      residual = std::max(residual, std::abs(updated[k] - result.freeEnergies[k]));
    }
    result.freeEnergies = updated;
  }

  // The density that belongs to the returned g_k: inside the loop it is computed from the g_k of the iteration
  // before (one step behind, and by the gauge shift of that step), and with no iteration at all it was never
  // computed, every bin left at the 'log of zero'.
  updateDensity();

  result.iterations = iteration;
  result.residual = residual;
  result.connected = statesAreConnected(counts, logStateCounts, result.logDensity, reducedPotential);
  result.converged = residual <= tolerance && result.connected;
  return result;
}

bool DcTSTWham::statesAreConnected(const std::vector<std::vector<double>>& counts,
                                   const std::vector<double>& logStateCounts, const std::vector<double>& logDensity,
                                   const std::vector<std::vector<double>>& reducedPotential)
{
  // The free-energy difference of two states is measured only where both of them sampled: two states are linked
  // by a bin that holds samples of both (and that neither forbids). The share n_k exp(g_k - u_k(b)) of a state in
  // the mixture of a bin is no such evidence. WHAM knows the samples of a bin by their sum only, so that two
  // states with fully disjoint histograms and a finite bias on each other's bins have a unique fixed point: for
  // counts {{n, 0}, {0, n}} and biases {{0, u}, {u, 0}} it is g = {0, 0} for every u, each state being credited
  // with the share e^-u of the samples of the other. The iteration converged on it, the shares were far above any
  // threshold, and a free-energy difference across a barrier that no state crossed was reported as converged.
  const std::size_t numberOfStates = logStateCounts.size();
  std::vector<std::size_t> parent(numberOfStates);
  std::iota(parent.begin(), parent.end(), 0uz);
  auto find = [&](std::size_t k)
  {
    while (parent[k] != k)
    {
      parent[k] = parent[parent[k]];
      k = parent[k];
    }
    return k;
  };

  for (std::size_t b = 0; b < logDensity.size(); ++b)
  {
    if (logDensity[b] <= logZeroThreshold) continue;
    std::optional<std::size_t> first;
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      if (logStateCounts[k] <= logZeroThreshold) continue;
      if (!(counts[k][b] > 0.0) || !(-reducedPotential[k][b] > forbiddenLogWeight)) continue;
      if (!first.has_value())
      {
        first = k;
        continue;
      }
      parent[find(k)] = find(first.value());
    }
  }

  std::optional<std::size_t> root;
  for (std::size_t k = 0; k < numberOfStates; ++k)
  {
    if (logStateCounts[k] <= logZeroThreshold) continue;
    if (!root.has_value())
    {
      root = find(k);
    }
    else if (find(k) != root.value())
    {
      return false;
    }
  }
  return true;
}
