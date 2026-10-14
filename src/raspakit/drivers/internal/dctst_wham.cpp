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

    const double gauge = updated[0];
    residual = 0.0;
    for (std::size_t k = 0; k < numberOfStates; ++k)
    {
      updated[k] -= gauge;
      residual = std::max(residual, std::abs(updated[k] - result.freeEnergies[k]));
    }
    result.freeEnergies = updated;
  }

  result.iterations = iteration;
  result.residual = residual;
  result.converged = residual <= tolerance;
  return result;
}
