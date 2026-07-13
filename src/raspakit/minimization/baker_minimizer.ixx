module;

export module minimization_baker;

import std;

import minimization_options;
import symmetric_eigensolver;

export struct BakerStep
{
  std::vector<double> displacement;
  std::vector<double> eigenvalues;
  double rmsGradient{};
  double maxGradient{};
  double stepNorm{};
  double shift{};
  double lowestEigenvalue{};
  double highestEigenvalue{};
  std::size_t negativeModes{};
  std::size_t zeroModes{};
  bool converged{};
  std::string convergenceReason;
};

/**
 * Compute one RASPA2-style Baker/LambdaMethod2 local-minimum step.
 * The Hessian is row-major and is symmetrized by the eigensolver.
 */
export BakerStep computeBakerStep(std::span<const double> hessian, std::span<const double> gradient,
                                 const MinimizationOptions &options);
