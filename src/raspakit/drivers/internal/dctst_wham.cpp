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

  // ln n_k and ln M_b
  std::vector<double> logStateCounts(numberOfStates, logZero);
  std::vector<double> logBinCounts(numberOfBins, logZero);
  for (std::size_t k = 0; k < numberOfStates; ++k)
  {
    double total{0.0};
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      total += counts[k][b];
    }
    if (total > 0.0) logStateCounts[k] = std::log(total);
  }
  for (std::size_t b = 0; b < numberOfBins; ++b)
  {
    double total{0.0};
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      total += counts[k][b];
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

  double residual = 1e30;
  std::size_t iteration = 0uz;
  for (; iteration < maximumIterations && residual > tolerance; ++iteration)
  {
    // Omega_b
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
      double sum{0.0};
      for (std::size_t k = 0; k < numberOfStates; ++k)
      {
        if (logStateCounts[k] <= logZeroThreshold) continue;
        sum += std::exp(logStateCounts[k] + result.freeEnergies[k] - reducedPotential[k][b] - largest);
      }
      result.logDensity[b] = logBinCounts[b] - largest - std::log(sum);
    }

    // g_k
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      double largest = logZero;
      for (std::size_t b = 0; b < numberOfBins; ++b)
      {
        if (result.logDensity[b] <= logZeroThreshold) continue;
        largest = std::max(largest, result.logDensity[b] - reducedPotential[k][b]);
      }
      if (largest <= logZeroThreshold)
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
      updated[k] -= gauge;
      if (logStateCounts[k] <= logZeroThreshold) continue;  // an unsampled state does not take part
      residual = std::max(residual, std::abs(updated[k] - result.freeEnergies[k]));
    }
    result.freeEnergies = updated;
  }

  result.iterations = iteration;
  result.residual = residual;
  result.connected = statesAreConnected(logStateCounts, result.freeEnergies, result.logDensity, reducedPotential);
  result.converged = residual <= tolerance && result.connected;
  return result;
}

bool DcTSTWham::statesAreConnected(const std::vector<double>& logStateCounts, const std::vector<double>& freeEnergies,
                                   const std::vector<double>& logDensity,
                                   const std::vector<std::vector<double>>& reducedPotential)
{
  // The WHAM equations fix the free-energy difference of two states only through bins both contribute to. When
  // the sampled states fall apart into groups that share no bin, every group is invariant under its own shift
  // of g and the iteration converges at once on an arbitrary offset between them. State k contributes to bin b
  // the share n_k exp(g_k - u_k(b)) / sum_l n_l exp(g_l - u_l(b)) of its mixture; two states are linked by a
  // bin where both shares exceed 'overlapShare'.
  constexpr double logOverlapShare = -27.631021115928547;  // ln 1e-12
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

  std::vector<double> logShare(numberOfStates);
  for (std::size_t b = 0; b < logDensity.size(); ++b)
  {
    if (logDensity[b] <= logZeroThreshold) continue;
    double largest = logZero;
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      logShare[k] = logZero;
      if (logStateCounts[k] <= logZeroThreshold) continue;
      logShare[k] = logStateCounts[k] + freeEnergies[k] - reducedPotential[k][b];
      largest = std::max(largest, logShare[k]);
    }
    if (largest <= logZeroThreshold) continue;
    double sum{0.0};
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      if (logShare[k] > logZeroThreshold) sum += std::exp(logShare[k] - largest);
    }
    const double logNormalization = largest + std::log(sum);
    std::optional<std::size_t> first;
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      if (logShare[k] <= logZeroThreshold || logShare[k] - logNormalization < logOverlapShare) continue;
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
