module;

export module property_molecule_backbone;

import std;

import archive;
import double3;
import atom;
import component;
import molecule_property_settings;

// Samples chain statistics along the backbone of every component that asks for them
// ('Component::moleculeBackboneSettings', 'ComputeMoleculeBackbone' in the 'Components' block; see
// Component::backboneAtoms: the shortest topological path between the end-to-end atoms), as
// functions of the separation k in backbone bonds, plus the single-chain form factor.
//
//   internal distances       <r^2(k)>  = < (1/(Nb-k)) sum_i |r_{i+k} - r_i|^2 >,   k = 1 .. Nb-1
//   bond-vector correlation  C(k)      = < (1/(Nb-1-k)) sum_i b^_i . b^_{i+k} >,  k = 0 .. Nb-2
//   form factor              P(q)      = < (1/N^2) sum_ij sin(q r_ij) / (q r_ij) >   (all atoms)
//
// with Nb the number of backbone beads, b_i = r_{i+1} - r_i the bond vectors and b^_i their unit
// vectors. The scalar chain descriptors derived from these are written to the summary:
//
//   mean bond length <l>, mean squared end-to-end distance <R^2> between the backbone ends,
//   characteristic ratio  C_N = <R^2> / ((Nb-1) <l>^2),
//   Kuhn length           b_K = <R^2> / R_max  and number of Kuhn segments  N_K = R_max / b_K,
//                         with R_max the contour length of the backbone,
//   persistence length    l_p (projection) = <sum_{j>=i} b^_i . b_j>  for an end bond i, averaged over both ends,
//                         l_p (fit)        from ln C(k) = -k <l> / l_p on the initial decay,
//   Flory exponent        nu from the log-log slope of <r^2(k)> over the window Nb/8 <= k <= Nb/2.
//
// The form factor is not evaluated per molecule (that would cost a sin() per pair per wave vector
// per sample and dominates the sampling time for a dense liquid). It depends only on the
// distribution of intramolecular pair distances, so a histogram h(r) of the pair distances is
// accumulated per block (one bin increment per pair) and transformed once when averages are formed:
//
//   P(q) = (1/N^2) [ N + (2/count) sum_bins h(r_b) sin(q r_b) / (q r_b) ]
//
// with r_b the bin centres and 'count' the number of sampled molecules. The bin width (0.005
// Angstrom) makes the discretization error of order (q dr)^2 / 24, below 1e-4 at q = 5 1/Angstrom.
// The histogram grows on demand, so no upper limit on the pair distance is needed.
//
// The form factor is compared against the Debye function P_D(x) = 2 (e^-x - 1 + x) / x^2,
// x = q^2 <Rg^2>, evaluated with the sampled <Rg^2> of the molecule. All averages carry 95%
// confidence intervals from block averaging. Output is written to 'molecule_backbone/'.

export struct PropertyMoleculeBackbone
{
  // Per-molecule scalar moments accumulated per [block][component]; 'numberOfCounts' normalizes.
  enum Moment : std::size_t
  {
    BondLength = 0,               ///< mean backbone bond length of the molecule [Angstrom]
    BondLengthSquared = 1,        ///< mean squared backbone bond length [Angstrom^2]
    EndToEndSquared = 2,          ///< squared distance between the backbone ends [Angstrom^2]
    Projection = 3,               ///< sum_{j>=i} b^_i . b_j for the end bonds, averaged over both ends [Angstrom]
    RadiusOfGyrationSquared = 4,  ///< Rg^2 over all atoms, uniform weights [Angstrom^2]
    NumberOfMoments = 5
  };
  using Moments = std::array<double, NumberOfMoments>;

  // Per-molecule averages of everything sampled (per block, or combined over all blocks).
  struct Averages
  {
    std::vector<double> internalDistanceSquared;  ///< index k-1 for k = 1 .. Nb-1
    std::vector<double> bondCorrelation;          ///< index k for k = 0 .. Nb-2
    std::vector<double> formFactor;               ///< per wave vector
    Moments moments{};
  };

  PropertyMoleculeBackbone() {};

  /// Builds the accumulators for every component with 'moleculeBackboneSettings'. Throws when such a
  /// component has a backbone of fewer than three beads or an invalid wave-vector range.
  PropertyMoleculeBackbone(std::size_t numberOfBlocks, const std::vector<Component> &components);

  std::uint64_t versionNumber{3};

  std::size_t numberOfBlocks{0};
  std::size_t numberOfComponents{0};

  // Per component: the settings (nullopt when the component is not sampled) and the logarithmically
  // spaced wave vectors [1/Angstrom] of its form factor.
  std::vector<std::optional<MoleculeBackboneSettings>> settingsPerComponent{};
  std::vector<std::vector<double>> waveVectorsPerComponent{};

  // Bin width [Angstrom] of the intramolecular pair-distance histogram behind the form factor.
  double pairDistanceBinWidth{0.005};

  // Per component: the backbone atom indices (empty when not sampled: fewer than three backbone
  // beads), its contour length R_max [Angstrom], and the number of atoms of the molecule.
  std::vector<std::vector<std::size_t>> backbonePerComponent{};
  std::vector<double> contourLengthPerComponent{};
  std::vector<std::size_t> numberOfAtomsPerComponent{};

  // Accumulators indexed as [block][component][k], [block][component][bin], [block][component].
  // 'pairDistanceHistogram' counts intramolecular pairs (i < j, all atoms) per distance bin; it is
  // resized on demand when a pair falls beyond the current range.
  std::vector<std::vector<std::vector<double>>> internalDistanceSquaredSum{};
  std::vector<std::vector<std::vector<double>>> bondCorrelationSum{};
  std::vector<std::vector<std::vector<double>>> pairDistanceHistogram{};
  std::vector<std::vector<Moments>> sums{};
  std::vector<std::vector<double>> numberOfCounts{};
  std::vector<double> totalNumberOfCounts{};  ///< Sampling events per component.

  bool isSampled(std::size_t component) const { return settingsPerComponent[component].has_value(); }
  std::size_t sampleEvery(std::size_t component) const { return settingsPerComponent[component]->sampleEvery; }
  std::optional<std::size_t> writeEvery(std::size_t component) const
  {
    return settingsPerComponent[component]->writeEvery;
  }
  std::size_t numberOfBackboneBeads(std::size_t component) const { return backbonePerComponent[component].size(); }

  // Contributions of a single molecule (exposed for testing). 'backbone' indexes into 'molecule'.
  static void accumulateInternalDistances(std::span<const Atom> molecule, std::span<const std::size_t> backbone,
                                          std::span<double> internalDistanceSquared);
  static void accumulateBondCorrelation(std::span<const Atom> molecule, std::span<const std::size_t> backbone,
                                        std::span<double> bondCorrelation);
  static void accumulatePairDistances(std::span<const Atom> molecule, double binWidth, std::vector<double> &histogram);
  static Moments computeMoments(std::span<const Atom> molecule, std::span<const std::size_t> backbone);

  // Direct per-molecule evaluation of the form factor (reference implementation, used in tests).
  static void accumulateFormFactor(std::span<const Atom> molecule, std::span<const double> waveVectors,
                                   std::span<double> formFactor);

  // Form factor from a pair-distance histogram of 'count' molecules of 'numberOfAtoms' atoms each.
  static std::vector<double> formFactorFromHistogram(std::span<const double> histogram, double binWidth,
                                                     std::size_t numberOfAtoms, double count,
                                                     std::span<const double> waveVectors);

  void sample(const std::vector<Component> &components,
              const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
              std::size_t currentCycle, std::size_t block);

  // Normalized per-molecule averages of one block, or combined over all blocks.
  Averages blockAverages(std::size_t block, std::size_t component) const;
  Averages overallAverages(std::size_t component) const;

  // Block-combined mean and 95% confidence-interval error of a function of the averages: the mean is
  // evaluated on the all-block averages, the error from the scatter of the per-block evaluations.
  std::pair<double, double> statistics(std::size_t component,
                                       const std::function<double(const Averages &)> &function) const;
  // Same, on precomputed averages ('blocks' holds one entry per block; empty blocks have zero counts).
  std::pair<double, double> statistics(std::size_t component, const Averages &overall,
                                       std::span<const Averages> blocks,
                                       const std::function<double(const Averages &)> &function) const;

  // Derived chain descriptors of a set of averages (exposed for testing).
  static double persistenceLengthFromProjection(const Averages &averages);
  static double persistenceLengthFromFit(const Averages &averages);
  static double floryExponent(const Averages &averages);
  static double debyeFunction(double x);

  void writeOutput(std::size_t systemId, const std::vector<Component> &components, std::size_t currentCycle);

  std::string printSettings() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyMoleculeBackbone &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyMoleculeBackbone &p);
};
