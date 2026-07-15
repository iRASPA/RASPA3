module;

export module minimization_options;

import std;

import int3;
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
  /// When enabled (and a flexible framework is present), the minimized cell is reduced to its primitive
  /// cell via SymmetryKit and the phonon dispersion is computed there, yielding the true unfolded band
  /// structure along the primitive Brillouin zone. The configured van der Waals cutoff is retained and
  /// periodic-image replicas are used to satisfy the minimum-image convention on the small primitive cell.
  /// Falls back to the simulation cell when symmetry detection fails.
  bool phononUsePrimitiveCell{true};
  /// Symmetry precision (Angstrom) used when detecting the primitive cell for the phonon computation.
  double phononPrimitiveCellSymmetryPrecision{1.0e-4};
  /// High-symmetry nodes (fractional reciprocal coordinates) defining the band-structure path. When empty
  /// and computePhononDispersion is enabled, a default G-X, G-Y, G-Z star is used (see input_reader).
  std::vector<PhononPathNode> phononDispersionPath{};
  /// When enabled, the phonon density of states is integrated over a uniform Brillouin-zone q-mesh (using the
  /// primitive cell when phononUsePrimitiveCell is set), producing a spectrum comparable to inelastic-neutron
  /// generalized phonon DOS. Uses the same primitive-cell reduction as the dispersion.
  bool computePhononDensityOfStates{false};
  /// Gamma-centered q-mesh dimensions for the density-of-states Brillouin-zone integration.
  int3 phononDOSMesh{8, 8, 8};
  /// Number of frequency bins in the density-of-states histogram.
  std::size_t phononDOSNumberOfBins{400};
  /// Gaussian broadening (cm^-1) applied to each mode frequency when accumulating the density of states.
  double phononDOSBroadening{5.0};
  double maximumStepLength{0.3};
  double maximumCellStepLength{0.1};
  double rmsGradientTolerance{1.0e-6};
  double maxGradientTolerance{1.0e-6};
  double convergenceFactor{1.0};
  double minimumEigenvalue{1.0e-3};
  double elasticEigenvalueTolerance{1.0e-8};
};
