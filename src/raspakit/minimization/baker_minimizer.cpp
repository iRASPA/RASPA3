module;

module minimization_baker;

import std;

namespace
{
double solveLambdaMethod2(std::span<const double> eigenvalues, std::span<const double> projectedGradient)
{
  if (eigenvalues.empty())
  {
    return 0.0;
  }

  constexpr std::size_t maximumIterations = 10000;
  constexpr double convergenceTolerance = 1.0e-12;
  const double lowest = eigenvalues.front();
  double lambda = 0.0;

  auto secularSum = [&](double trial)
  {
    double sum = 0.0;
    for (std::size_t i = 0; i < eigenvalues.size(); ++i)
    {
      double denominator = trial - eigenvalues[i];
      if (std::abs(denominator) < 1.0e-12)
      {
        denominator = std::copysign(1.0e-12, denominator == 0.0 ? -1.0 : denominator);
      }
      sum += projectedGradient[i] * projectedGradient[i] / denominator;
    }
    return sum;
  };

  if (lowest < 0.0)
  {
    double upper = lowest - 1.0e-3;
    double lower = lowest - 500.0;
    for (std::size_t iteration = 0; iteration < maximumIterations; ++iteration)
    {
      lambda = 0.5 * (lower + upper);
      const double residual = secularSum(lambda) - lambda;
      if (std::abs(residual) <= convergenceTolerance * std::max(1.0, std::abs(lambda)))
      {
        return lambda;
      }
      if (residual > 0.0)
      {
        lower = lambda;
      }
      else
      {
        upper = lambda;
      }
    }
    return lambda;
  }

  for (std::size_t iteration = 0; iteration < maximumIterations; ++iteration)
  {
    const double next = secularSum(lambda);
    if (std::abs(next - lambda) <= convergenceTolerance * std::max(1.0, std::abs(lambda)))
    {
      return next;
    }
    lambda = next;
  }
  return lambda;
}
}  // namespace

BakerStep computeBakerStep(std::span<const double> hessian, std::span<const double> gradient,
                           const MinimizationOptions &options)
{
  const std::size_t size = gradient.size();
  if (hessian.size() != size * size)
  {
    throw std::invalid_argument("computeBakerStep: Hessian and gradient extents do not match");
  }
  BakerStep result{};
  result.displacement.assign(size, 0.0);
  if (size == 0)
  {
    result.converged = true;
    result.convergenceReason = "no molecular degrees of freedom";
    return result;
  }
  for (double value : gradient)
  {
    if (!std::isfinite(value))
    {
      throw std::runtime_error("computeBakerStep: non-finite gradient");
    }
  }

  double gradientNormSquared = 0.0;
  for (double value : gradient)
  {
    gradientNormSquared += value * value;
    result.maxGradient = std::max(result.maxGradient, std::abs(value));
  }
  result.rmsGradient = std::sqrt(gradientNormSquared) / static_cast<double>(size);

  const SymmetricEigenSystem eigenSystem = diagonalizeSymmetric(hessian, size);
  result.eigenvalues = eigenSystem.eigenvalues;
  result.lowestEigenvalue = eigenSystem.eigenvalues.front();
  result.highestEigenvalue = eigenSystem.eigenvalues.back();
  result.lowestMode.resize(size);
  for (std::size_t row = 0; row < size; ++row)
  {
    result.lowestMode[row] = eigenSystem.eigenvector(row, 0);
  }

  std::vector<std::size_t> activeModes;
  std::vector<double> activeEigenvalues;
  std::vector<double> projectedGradient;
  for (std::size_t mode = 0; mode < size; ++mode)
  {
    const double eigenvalue = eigenSystem.eigenvalues[mode];
    if (eigenvalue < -options.minimumEigenvalue)
    {
      ++result.negativeModes;
    }
    if (std::abs(eigenvalue) < options.minimumEigenvalue)
    {
      ++result.zeroModes;
      continue;
    }
    double projection = 0.0;
    for (std::size_t row = 0; row < size; ++row)
    {
      projection += eigenSystem.eigenvector(row, mode) * gradient[row];
    }
    activeModes.push_back(mode);
    activeEigenvalues.push_back(eigenvalue);
    projectedGradient.push_back(projection);
  }

  result.converged = result.rmsGradient < options.rmsGradientTolerance &&
                     result.maxGradient < options.maxGradientTolerance && result.negativeModes == 0;
  if (result.converged)
  {
    result.convergenceReason = "gradient tolerances met with zero negative modes";
  }
  else if (result.negativeModes != 0)
  {
    result.convergenceReason = "negative Hessian modes remain";
  }
  else
  {
    result.convergenceReason = "gradient tolerances not met";
  }
  if (result.converged || activeModes.empty())
  {
    return result;
  }

  result.shift = solveLambdaMethod2(activeEigenvalues, projectedGradient);
  for (std::size_t active = 0; active < activeModes.size(); ++active)
  {
    const double denominator = result.shift - activeEigenvalues[active];
    if (std::abs(denominator) <= options.minimumEigenvalue)
    {
      continue;
    }
    const double scale = projectedGradient[active] / denominator;
    for (std::size_t row = 0; row < size; ++row)
    {
      result.displacement[row] += scale * eigenSystem.eigenvector(row, activeModes[active]);
    }
  }

  for (double value : result.displacement)
  {
    result.stepNorm += value * value;
  }
  result.stepNorm = std::sqrt(result.stepNorm);
  const double adaptiveMaximum =
      std::min(options.maximumStepLength, std::pow(result.rmsGradient, options.convergenceFactor));
  if (result.stepNorm > adaptiveMaximum && result.stepNorm > 0.0)
  {
    const double scale = adaptiveMaximum / result.stepNorm;
    for (double &value : result.displacement)
    {
      value *= scale;
    }
    result.stepNorm = adaptiveMaximum;
  }
  return result;
}
