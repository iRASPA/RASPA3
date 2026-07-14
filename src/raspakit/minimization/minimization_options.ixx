module;

export module minimization_options;

import std;

export struct MinimizationOptions
{
  std::size_t maximumNumberOfSteps{10000};
  std::size_t printEvery{1};
  bool computeElasticConstants{false};
  double maximumStepLength{0.3};
  double maximumCellStepLength{0.1};
  double rmsGradientTolerance{1.0e-6};
  double maxGradientTolerance{1.0e-6};
  double convergenceFactor{1.0};
  double minimumEigenvalue{1.0e-3};
  double elasticEigenvalueTolerance{1.0e-8};
};
