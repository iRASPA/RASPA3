module;

export module property_polymer_chains;

import std;

import archive;
import double3;
import atom;
import forcefield;
import component;
import simulationbox;
import polymer_links;

/**
 * \brief Size and shape of the chains of a polymerization, resolved per chain length ('ComputePolymerChains').
 *
 * The chains are the connected clusters of the atom graph of the PolymerLinkTable: the bonds of the molecules that
 * no reaction broke plus the polymer bonds. The length of a chain is the number of molecules it holds (the degree
 * of polymerization for a chain of monomers); a molecule without a polymer bond is not a chain. Every chain is
 * unwrapped by walking its bonds with minimum-image bond vectors and its gyration tensor is diagonalised (uniform
 * bead weights, or mass weights with 'MassWeightedPolymerChains'), giving Rg^2, the eigenvalues, the asphericity,
 * the acylindricity, the relative shape anisotropy k^2 and the prolateness (see PropertyMoleculeShape). The
 * end-to-end vector joins the two ends of the longest shortest path of the bond graph of the chain (a double
 * breadth-first sweep: exact for a tree, a lower bound in a network); for a vinyl chain with side groups the path
 * ends in the side groups of the two terminal monomers, a constant offset that matters less as the chain grows.
 * The number of bonds along that path is reported as the contour.
 *
 * The idea is to measure the polymer as it forms rather than to equilibrate it afterwards: at 300 K a long chain
 * does not relax on the time scale of a simulation, but every polymerization run grows a new set of independent
 * chains, so the as-polymerized structure of every length is sampled by repeating the run. The averages are kept as
 * raw sums per key, so that the runs are combined by adding them (the JSON output; 'raspalib.polymer_chains'
 * combines files and computes the error from the scatter between the runs). The key of a sample is
 *
 *   length          the number of molecules of the chain,
 *   topology        linear (a path of molecules), branched (a tree in which a molecule has more than two polymer
 *                   bonds: a cross-link) or cyclic (more polymer bonds than a tree has),
 *   living          whether a group of the chain is in an active state of the scheme (a living radical),
 *   conversion bin  the conversion at the sample in 'ConversionBinsPolymerChains' equal bins over [0, 1]: an N-mer
 *                   grown in its monomer is in another environment than one grown in a polymer melt.
 *
 * A chain is measured only when its topology has been constant for the force ramp ('RelaxationSteps') plus
 * 'LagStepsPolymerChains' steps: the last change is the latest formation of one of its bonds, break of one of its
 * molecule bonds, retyping of one of its atoms or state change of one of its groups. Chains that wind through the
 * periodic boundaries (a closed bond cycle that is not closed in space: a spanning cluster, the gel) have no radius
 * of gyration and are skipped; both skips are counted.
 *
 * Sampling ('SamplingPolymerChains'):
 *
 *   'Periodic'        every chain every 'SamplePolymerChainsEvery' MD steps: an N-mer is weighted by the time it
 *                     spends at length N (for dead chains, also by how long the run continues after their death);
 *   'BeforeReaction'  a chain is sampled just before a reaction changes it (the most relaxed state of its length)
 *                     and every chain once more at the end of the run (the dead chains, which never react again):
 *                     every length of a chain counts once. A length whose ramp had not ended when the chain grew
 *                     again is not sampled (skipped while ramping).
 *
 * Output (every 'WritePolymerChainsEvery' cycles and at the end) in 'polymer_chains/': 'polymer_chains.s<i>.txt',
 * the averages per length and topology (pooled over living and conversion); 'polymer_chains_resolved.s<i>.txt',
 * the averages per key; 'polymer_chains.s<i>.json', the raw sums per key for combining runs. The standard deviation
 * in the text files is the width of the distribution over the samples, not the error of the mean (the samples of
 * one run are correlated; the error follows from the scatter between independent runs).
 */
export struct PropertyPolymerChains
{
  enum class Sampling : std::uint8_t
  {
    Periodic = 0,
    BeforeReaction = 1
  };

  enum class Topology : std::uint8_t
  {
    Linear = 0,
    Branched = 1,
    Cyclic = 2
  };

  /// The per-chain values whose sums and sums of squares are accumulated per key.
  enum Moment : std::size_t
  {
    RadiusOfGyrationSquared = 0,  ///< Rg^2 [Angstrom^2]
    RadiusOfGyration = 1,         ///< Rg [Angstrom]
    Lambda1 = 2,                  ///< largest eigenvalue of the gyration tensor [Angstrom^2]
    Lambda2 = 3,                  ///< middle eigenvalue [Angstrom^2]
    Lambda3 = 4,                  ///< smallest eigenvalue [Angstrom^2]
    Asphericity = 5,              ///< b [Angstrom^2]
    Acylindricity = 6,            ///< c [Angstrom^2]
    ShapeAnisotropy = 7,          ///< k^2 [-]
    Prolateness = 8,              ///< S [-]
    EndToEndSquared = 9,          ///< R^2 [Angstrom^2]
    EndToEnd = 10,                ///< R [Angstrom]
    ContourBonds = 11,            ///< bonds along the end-to-end path [-]
    NumberOfAtoms = 12,           ///< [-]
    Mass = 13,                    ///< [g/mol]
    NumberOfMoments = 14
  };
  using Moments = std::array<double, NumberOfMoments>;

  static constexpr std::array<std::string_view, NumberOfMoments> momentNames{
      "Rg^2", "Rg", "lambda1", "lambda2", "lambda3", "b", "c", "k^2", "S", "R^2", "R", "contour", "atoms", "mass"};
  static constexpr std::array<std::string_view, NumberOfMoments> momentUnits{
      "A^2", "A", "A^2", "A^2", "A^2", "A^2", "A^2", "-", "-", "A^2", "A", "bonds", "-", "g/mol"};

  struct Key
  {
    std::size_t length{0};
    Topology topology{Topology::Linear};
    bool living{false};
    std::size_t conversionBin{0};

    auto operator<=>(const Key &) const = default;
    bool operator==(const Key &) const = default;

    friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Key &k);
    friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Key &k);
  };

  struct Accumulator
  {
    std::uint64_t count{0};
    Moments sum{};
    Moments sumOfSquares{};

    void add(const Moments &values);
    Accumulator &operator+=(const Accumulator &other);
    double mean(Moment moment) const { return count > 0 ? sum[moment] / static_cast<double>(count) : 0.0; }
    /// The standard deviation of the samples (the width of the distribution, not the error of the mean).
    double standardDeviation(Moment moment) const;

    bool operator==(const Accumulator &) const = default;

    friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Accumulator &a);
    friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Accumulator &a);
  };

  /// One chain of the current configuration (see findChains()).
  struct Chain
  {
    std::vector<std::uint32_t> atoms{};  ///< Global atom indices, in the order of the walk.
    std::vector<double3> positions{};    ///< The unwrapped positions of 'atoms'.
    std::size_t length{0};               ///< Molecules.
    std::size_t polymerBonds{0};
    Topology topology{Topology::Linear};
    bool living{false};
    bool spanning{false};  ///< A bond cycle of the chain winds through the periodic boundaries.
    std::uint64_t lastChange{0};  ///< The MD step of the last change of its topology.
    std::size_t contourBonds{0};
    double3 endToEnd{};
  };

  PropertyPolymerChains() {}
  PropertyPolymerChains(std::size_t sampleEvery, std::size_t writeEvery, Sampling sampling, std::size_t lagSteps,
                        std::size_t numberOfConversionBins, std::size_t minimumLength, bool massWeighted);

  std::uint64_t versionNumber{1};

  std::size_t sampleEvery{1000};  ///< MD steps between the samples ('Periodic').
  std::size_t writeEvery{5000};   ///< Cycles between the writes of the output files.
  Sampling sampling{Sampling::Periodic};
  std::size_t lagSteps{0};               ///< MD steps after the force ramp before a chain is measured.
  std::size_t numberOfConversionBins{10};
  std::size_t minimumLength{2};          ///< Shorter chains (in molecules) are not sampled.
  bool massWeighted{false};

  std::map<Key, Accumulator> accumulators{};
  std::uint64_t numberOfSamplingEvents{0};
  std::uint64_t skippedRamping{0};   ///< Chains not measured because their topology changed too recently.
  std::uint64_t skippedSpanning{0};  ///< Chains not measured because they span the periodic boundaries.

  /**
   * \brief The chains of the current configuration: all of them, or only those that contain one of 'seeds' (global
   * atom indices). Needs an activated table; the adjacency is that of the table (neighboursOf()), which follows a
   * reaction at once.
   */
  static std::vector<Chain> findChains(const PolymerLinkTable &table, const std::vector<Component> &components,
                                       const SimulationBox &simulationBox, std::span<const Atom> atoms,
                                       std::span<const std::uint32_t> seeds = {});

  /// The values of one chain.
  static Moments measure(const Chain &chain, const ForceField &forceField, std::span<const Atom> atoms,
                         bool massWeighted);

  std::size_t conversionBin(double conversion) const;

  /// One sampling event: measures the eligible chains (all, or those containing one of 'seeds') at the current
  /// step of the table.
  void sample(const PolymerLinkTable &table, const ForceField &forceField, const std::vector<Component> &components,
              const SimulationBox &simulationBox, std::span<const Atom> atoms,
              std::span<const std::uint32_t> seeds = {});

  /// The accumulators pooled over living and conversion, per length and topology.
  std::map<std::pair<std::size_t, Topology>, Accumulator> pooled() const;

  /// Writes the output files when 'currentCycle' is a multiple of 'writeEvery' (or always with 'force').
  void writeOutput(std::size_t systemId, std::size_t currentCycle, const PolymerLinkTable &table, bool force) const;

  /// The line of the status report.
  std::string printStatus() const;
  std::string printSettings() const;

  static std::string_view topologyName(Topology topology);
  static std::string_view samplingName(Sampling sampling);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPolymerChains &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPolymerChains &p);
};
