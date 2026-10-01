module;

export module property_end_to_end_acf;

import std;

import archive;
import double3;
import atom;
import molecule;
import component;

/// One lag of the end-to-end vector autocorrelation function of a component.
export struct EndToEndAutoCorrelationFunctionData
{
  double time{0.0};             ///< Lag [ps] (or cycles when the system has no time step).
  double acf{0.0};              ///< <R(0).R(t)> [Angstrom^2].
  double normalized{0.0};       ///< <R(0).R(t)> / <R^2> [-].
  double numberOfSamples{0.0};  ///< Number of (molecule, time-origin) pairs that entered the average.
};

/// Relaxation-time estimates derived from the normalized autocorrelation function C(t) / C(0).
export struct EndToEndRelaxationTimes
{
  double meanSquaredEndToEnd{0.0};  ///< C(0) = <R^2> [Angstrom^2].
  double longestLag{0.0};           ///< The longest lag with data.

  /// Integrated correlation time, the integral of C(t)/C(0) up to the first zero crossing (or to the longest
  /// lag when the function has not crossed zero yet, in which case 'integratedIsLowerBound' is set).
  std::optional<double> integrated{};
  bool integratedIsLowerBound{false};

  /// The time at which C(t)/C(0) first falls to 1/e (linearly interpolated).
  std::optional<double> oneOverE{};

  /// Decay time of a single exponential fitted (least squares in log space) to the lags where
  /// 0.05 < C(t)/C(0) <= 0.5, the range dominated by the slowest (Rouse p = 1) mode.
  std::optional<double> exponentialFit{};
};

/**
 * \brief End-to-end vector autocorrelation function <R(0).R(t)> and the end-to-end relaxation time tau_R.
 *
 * For every component with end-to-end atoms (see Component::endToEndAtoms) the vector R between the two ends
 * (unwrapped positions) is recorded for every molecule, and the autocorrelation function
 * \f[ C(t) = \langle \mathbf{R}(0) \cdot \mathbf{R}(t) \rangle \f]
 * is accumulated with the order-N (Frenkel & Smit) blocking scheme that the mean-squared displacement uses: block
 * b holds the lags k * sampleEvery * n^b (k = 1 .. n-1), so a run covers lags from one sampling interval to
 * sampleEvery * n^(blocks) with a logarithmic density of points, at a fixed memory cost. C(0) = <R^2>.
 *
 * The number of molecules must be fixed (the molecules are tracked by index); the property is sampled in
 * System::sampleProperties and written by the molecular-dynamics drivers. In a Monte Carlo run the "time" is the
 * number of cycles, which measures the decorrelation of the chain conformations by the Monte Carlo moves.
 *
 * The relaxation time tau_R (the time the chain needs to forget its end-to-end orientation, i.e. the spacing of
 * independent samples of R) is estimated in three ways in relaxationTimes(): the integral of the normalized
 * function, its 1/e time, and an exponential fit of the tail.
 */
export struct PropertyEndToEndAutoCorrelationFunction
{
  PropertyEndToEndAutoCorrelationFunction() = default;

  PropertyEndToEndAutoCorrelationFunction(
      const std::vector<std::size_t> &numberOfMoleculesPerComponent,
      const std::vector<std::optional<std::array<std::size_t, 2>>> &endToEndAtomsPerComponent,
      std::size_t numberOfParticles, double timeStep, std::size_t numberOfBlockElements, std::size_t sampleEvery,
      std::optional<std::size_t> writeEvery);

  std::uint64_t versionNumber{1};

  std::vector<std::size_t> numberOfMoleculesPerComponent{};
  std::vector<std::optional<std::array<std::size_t, 2>>> endToEndAtomsPerComponent{};
  std::size_t numberOfComponents{0uz};
  std::size_t numberOfParticles{0uz};
  double timeStep{0.0};
  std::size_t numberOfBlockElements{25uz};
  std::size_t sampleEvery{1uz};
  std::optional<std::size_t> writeEvery{};

  std::size_t count{0uz};             ///< Number of samples taken.
  std::size_t numberOfBlocks{0uz};    ///< Blocks in use for the current count.
  std::size_t maxNumberOfBlocks{1uz}; ///< Blocks allocated.
  std::vector<std::size_t> blockLength{};

  std::vector<std::vector<std::vector<std::size_t>>> acfCount{};  ///< [block][component][k]
  std::vector<std::vector<std::vector<double3>>> blockData{};     ///< [block][molecule][k] end-to-end history
  std::vector<std::vector<std::vector<double>>> acf{};            ///< [block][component][k] sum of R(0).R(t)

  /// Samples the end-to-end vectors of all molecules (every 'sampleEvery' cycles).
  void addSample(std::size_t currentCycle, const std::vector<Molecule> &molecules, std::span<const Atom> atoms);

  /// Accumulates one sample from the end-to-end vectors given per molecule (in molecule order; the entries of
  /// components without end-to-end atoms are ignored).
  void addSampleVectors(std::span<const double3> endToEndVectors);

  /// Whether the component has end-to-end atoms and therefore an autocorrelation function.
  bool hasData(std::size_t component) const;

  /// The autocorrelation function of a component, in order of increasing lag (lag 0 first).
  std::vector<EndToEndAutoCorrelationFunctionData> result(std::size_t component) const;

  /// Relaxation-time estimates of a component from its current autocorrelation function.
  EndToEndRelaxationTimes relaxationTimes(std::size_t component) const;
  static EndToEndRelaxationTimes relaxationTimes(const std::vector<EndToEndAutoCorrelationFunctionData> &data);

  /// Lag of entry k of block 'block' in the time unit of the run (ps, or cycles without a time step).
  double lagOf(std::size_t block, std::size_t k) const;
  bool usesCycles() const { return !(timeStep > 0.0); }
  std::string timeUnit() const { return usesCycles() ? "cycles" : "ps"; }

  void writeOutput(std::size_t systemId, const std::vector<Component> &components, std::size_t currentCycle) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyEndToEndAutoCorrelationFunction &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEndToEndAutoCorrelationFunction &p);
};
