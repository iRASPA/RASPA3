module;

export module property_end_to_end_acf;

import std;

import archive;
import double3;
import atom;
import molecule;
import component;
import molecule_property_settings;

/// One lag of the end-to-end autocorrelation functions of a component.
export struct EndToEndAutoCorrelationFunctionData
{
  double time{0.0};             ///< Lag [ps] (or cycles when the system has no time step).
  double acf{0.0};              ///< <R(0).R(t)> [Angstrom^2].
  double normalized{0.0};       ///< <R(0).R(t)> / <R^2> [-].
  double numberOfSamples{0.0};  ///< Number of (molecule, time-origin) pairs that entered the average.

  /// Rotation-invariant (conformational) autocorrelation functions of the fluctuations, normalized to one at
  /// zero lag: <dR^2(0) dR^2(t)> / <(dR^2)^2> with dR^2 = R^2 - <R^2>, and the same for the squared radius of
  /// gyration. Empty when the quantity does not fluctuate (zero variance) or was not sampled.
  std::optional<double> endToEndSquaredACF{};
  std::optional<double> radiusOfGyrationSquaredACF{};
};

/// Relaxation-time estimates derived from a normalized autocorrelation function C(t) / C(0).
export struct RelaxationTimeEstimates
{
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

/// Relaxation of a rotation-invariant scalar x (R^2 or Rg^2) from the autocorrelation of its fluctuations.
export struct ScalarRelaxation
{
  bool available{false};  ///< False when the scalar does not fluctuate (or was not sampled).
  double mean{0.0};       ///< <x> [Angstrom^2].
  double variance{0.0};   ///< <x^2> - <x>^2 [Angstrom^4].
  RelaxationTimeEstimates times{};
};

/// Relaxation-time estimates of the end-to-end vector and of the conformational scalars.
export struct EndToEndRelaxationTimes
{
  double meanSquaredEndToEnd{0.0};  ///< C(0) = <R^2> [Angstrom^2].
  double longestLag{0.0};           ///< The longest lag with data.

  // the vector function <R(0).R(t)> / <R^2> (see RelaxationTimeEstimates for the meaning of the three)
  std::optional<double> integrated{};
  bool integratedIsLowerBound{false};
  std::optional<double> oneOverE{};
  std::optional<double> exponentialFit{};

  /// Conformational relaxation: the fluctuation autocorrelation of R^2 and Rg^2. The vector function above also
  /// decays by rigid-body tumbling of an unchanged conformation; these two do not, so they measure when the chain
  /// actually changes shape (e.g. when a hairpin opens).
  ScalarRelaxation endToEndSquared{};
  ScalarRelaxation radiusOfGyrationSquared{};
};

/**
 * \brief End-to-end vector autocorrelation function <R(0).R(t)> and the end-to-end relaxation time tau_R.
 *
 * For every component that asks for it (Component::endToEndACFSettings, 'ComputeEndToEndACF' in the 'Components'
 * block; the component needs end-to-end atoms, see Component::endToEndAtoms) the vector R between the two ends
 * (unwrapped positions) is recorded for every molecule, and the autocorrelation function
 * \f[ C(t) = \langle \mathbf{R}(0) \cdot \mathbf{R}(t) \rangle \f]
 * is accumulated with the order-N (Frenkel & Smit) blocking scheme that the mean-squared displacement uses: block
 * b holds the lags k * sampleEvery * n^b (k = 1 .. n-1), so a run covers lags from one sampling interval to
 * sampleEvery * n^(blocks) with a logarithmic density of points, at a fixed memory cost. C(0) = <R^2>. Every
 * component has its own sampling interval, block size and accumulators.
 *
 * The number of molecules must be fixed (the molecules are tracked by index); the property is sampled in
 * System::sampleProperties and written by the molecular-dynamics drivers. In a Monte Carlo run the "time" is the
 * number of cycles, which measures the decorrelation of the chain conformations by the Monte Carlo moves.
 *
 * The relaxation time tau_R (the time the chain needs to forget its end-to-end orientation, i.e. the spacing of
 * independent samples of R) is estimated in three ways in relaxationTimes(): the integral of the normalized
 * function, its 1/e time, and an exponential fit of the tail.
 *
 * The vector function also decays when a chain of fixed shape merely tumbles, so alongside it the same order-N
 * scheme accumulates the autocorrelation of two rotation-invariant scalars, the squared end-to-end distance R^2
 * and the squared radius of gyration Rg^2 (uniform weights over all atoms of the molecule), as fluctuation
 * functions
 * \f[ C_x(t) = \frac{\langle \delta x(0)\, \delta x(t) \rangle}{\langle \delta x^2 \rangle}, \quad
 *     \delta x = x - \langle x \rangle . \f]
 * These stay at one as long as the conformation is unchanged and measure the conformational relaxation (the
 * opening of a fold, the migration of a kink), which for a compact oligomer in a viscous liquid can be much
 * slower than the tumbling time that dominates <R(0).R(t)>.
 */
export struct PropertyEndToEndAutoCorrelationFunction
{
  /// The rotation-invariant scalars accumulated next to the end-to-end vector.
  enum Scalar : std::size_t
  {
    EndToEndSquared = 0,          ///< R^2 [Angstrom^2]
    RadiusOfGyrationSquared = 1,  ///< Rg^2 [Angstrom^2], uniform weights over all atoms
    NumberOfScalars = 2
  };
  using Scalars = std::array<double, NumberOfScalars>;

  /// Order-N accumulators of one component.
  struct ComponentData
  {
    std::size_t count{0uz};              ///< Number of samples taken.
    std::size_t numberOfBlocks{0uz};     ///< Blocks in use for the current count.
    std::size_t maxNumberOfBlocks{0uz};  ///< Blocks allocated.
    std::vector<std::size_t> blockLength{};                   ///< [block]
    std::vector<std::vector<std::size_t>> acfCount{};         ///< [block][k]
    std::vector<std::vector<std::vector<double3>>> blockData{};  ///< [block][molecule of the component][k]
    std::vector<std::vector<double>> acf{};                   ///< [block][k] sum of R(0).R(t)

    std::vector<std::vector<std::vector<Scalars>>> blockScalars{};  ///< [block][molecule][k] history of the scalars
    std::vector<std::vector<Scalars>> acfScalars{};                 ///< [block][k] sum of x(0) x(t) per scalar
    Scalars sumScalars{};                                           ///< sum of x over all (molecule, sample) pairs
    std::size_t numberOfScalarSamples{0uz};                         ///< number of (molecule, sample) pairs
    std::array<bool, NumberOfScalars> scalarSampled{};              ///< whether the scalar was supplied for every sample

    friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ComponentData &d);
    friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ComponentData &d);
  };

  PropertyEndToEndAutoCorrelationFunction() = default;

  /// Builds the accumulators for every component with 'endToEndACFSettings'. Throws when such a component has no
  /// end-to-end atoms.
  PropertyEndToEndAutoCorrelationFunction(const std::vector<Component> &components,
                                          const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                          std::size_t numberOfParticles, double timeStep);

  /// Direct form: the end-to-end atoms and settings per component (nullopt settings: not sampled).
  PropertyEndToEndAutoCorrelationFunction(
      const std::vector<std::size_t> &numberOfMoleculesPerComponent,
      const std::vector<std::optional<std::array<std::size_t, 2>>> &endToEndAtomsPerComponent,
      const std::vector<std::optional<EndToEndACFSettings>> &settingsPerComponent, std::size_t numberOfParticles,
      double timeStep);

  std::uint64_t versionNumber{3};

  std::vector<std::size_t> numberOfMoleculesPerComponent{};
  std::vector<std::size_t> moleculeOffsetPerComponent{};  ///< Index of the first molecule of each component.
  std::vector<std::optional<std::array<std::size_t, 2>>> endToEndAtomsPerComponent{};
  std::vector<std::optional<EndToEndACFSettings>> settingsPerComponent{};
  std::size_t numberOfComponents{0uz};
  std::size_t numberOfParticles{0uz};
  double timeStep{0.0};

  std::vector<ComponentData> dataPerComponent{};  ///< One per component (empty for components not sampled).

  bool isSampled(std::size_t component) const
  {
    return component < numberOfComponents && settingsPerComponent[component].has_value();
  }
  std::size_t sampleEvery(std::size_t component) const { return settingsPerComponent[component]->sampleEvery; }
  std::optional<std::size_t> writeEvery(std::size_t component) const
  {
    return settingsPerComponent[component]->writeEvery;
  }
  std::size_t numberOfBlockElements(std::size_t component) const
  {
    return settingsPerComponent[component]->numberOfBlockElements;
  }

  /// Samples the end-to-end vectors and the squared radii of gyration of the molecules of every component due
  /// this cycle.
  void addSample(std::size_t currentCycle, const std::vector<Molecule> &molecules, std::span<const Atom> atoms);

  /// Accumulates one sample for every sampled component from the end-to-end vectors given per molecule (in
  /// molecule order; the entries of components that are not sampled are ignored). 'radiusOfGyrationSquared' is
  /// optional (same indexing); a component for which it is ever omitted has no Rg^2 channel.
  void addSampleVectors(std::span<const double3> endToEndVectors, std::span<const double> radiusOfGyrationSquared = {});

  /// Accumulates one sample of one component from the end-to-end vectors (and optionally the squared radii of
  /// gyration) given per molecule (all molecules, in molecule order).
  void addSampleVectors(std::size_t component, std::span<const double3> endToEndVectors,
                        std::span<const double> radiusOfGyrationSquared = {});

  /// Whether the component is sampled and has molecules, and therefore an autocorrelation function.
  bool hasData(std::size_t component) const;

  /// The autocorrelation functions of a component, in order of increasing lag (lag 0 first).
  std::vector<EndToEndAutoCorrelationFunctionData> result(std::size_t component) const;

  /// Relaxation-time estimates of a component from its current autocorrelation functions.
  EndToEndRelaxationTimes relaxationTimes(std::size_t component) const;
  static EndToEndRelaxationTimes relaxationTimes(const std::vector<EndToEndAutoCorrelationFunctionData> &data);

  /// The three relaxation-time estimates of any normalized autocorrelation function given as (lag, C(t)/C(0))
  /// pairs in increasing lag, lag 0 first.
  static RelaxationTimeEstimates estimateRelaxationTimes(std::span<const double> lags,
                                                         std::span<const double> normalized);

  /// Squared radius of gyration of a molecule from unwrapped positions, uniform weights over all atoms.
  static double radiusOfGyrationSquared(const Molecule &molecule, std::span<const Atom> atoms);

  /// Lag of entry k of block 'block' of a component in the time unit of the run (ps, or cycles without a time
  /// step).
  double lagOf(std::size_t component, std::size_t block, std::size_t k) const;
  bool usesCycles() const { return !(timeStep > 0.0); }
  std::string timeUnit() const { return usesCycles() ? "cycles" : "ps"; }

  void writeOutput(std::size_t systemId, const std::vector<Component> &components, std::size_t currentCycle) const;

  std::string printSettings() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyEndToEndAutoCorrelationFunction &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyEndToEndAutoCorrelationFunction &p);
};
