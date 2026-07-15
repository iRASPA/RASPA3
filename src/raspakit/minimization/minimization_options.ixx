module;

export module minimization_options;

import std;

import phonon_kpath;

export struct MinimizationOptions
{
  std::size_t maximumNumberOfSteps{10000};
  std::size_t printEvery{1};
  bool computeElasticConstants{false};
  bool computeNormalModes{false};
  bool normalModeMovies{false};
  std::size_t normalModeMoviePeriods{1};
  std::size_t normalModeMoviePointsPerPeriod{16};
  double normalModeMovieAmplitude{0.5};
  bool computePhononDispersion{false};
  std::size_t phononDispersionPointsPerSegment{20};
  /// High-symmetry nodes (fractional reciprocal coordinates) defining the band-structure path. When empty
  /// and computePhononDispersion is enabled, a default G-X, G-Y, G-Z star is used (see input_reader).
  std::vector<PhononPathNode> phononDispersionPath{};
  double maximumStepLength{0.3};
  double maximumCellStepLength{0.1};
  double rmsGradientTolerance{1.0e-6};
  double maxGradientTolerance{1.0e-6};
  double convergenceFactor{1.0};
  double minimumEigenvalue{1.0e-3};
  double elasticEigenvalueTolerance{1.0e-8};
};
