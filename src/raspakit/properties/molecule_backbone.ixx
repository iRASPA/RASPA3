module;

export module property_molecule_backbone;

import std;

import archive;
import double3;
import atom;
import component;

// Samples chain statistics along the backbone of every component that has one (see
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

  PropertyMoleculeBackbone(std::size_t numberOfBlocks, const std::vector<Component> &components,
                          std::size_t numberOfWaveVectors, double waveVectorLowerLimit, double waveVectorUpperLimit,
                          std::size_t sampleEvery, std::optional<std::size_t> writeEvery);

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks{0};
  std::size_t numberOfComponents{0};
  std::size_t sampleEvery{10};
  std::optional<std::size_t> writeEvery{5000};

  // Logarithmically spaced wave vectors [1/Angstrom] shared by all components.
  std::size_t numberOfWaveVectors{64};
  double waveVectorLowerLimit{0.01};
  double waveVectorUpperLimit{5.0};
  std::vector<double> waveVectors{};

  // Per component: the backbone atom indices (empty when not sampled: fewer than three backbone
  // beads) and its contour length R_max [Angstrom].
  std::vector<std::vector<std::size_t>> backbonePerComponent{};
  std::vector<double> contourLengthPerComponent{};

  // Accumulators indexed as [block][component][k], [block][component][q], [block][component].
  std::vector<std::vector<std::vector<double>>> internalDistanceSquaredSum{};
  std::vector<std::vector<std::vector<double>>> bondCorrelationSum{};
  std::vector<std::vector<std::vector<double>>> formFactorSum{};
  std::vector<std::vector<Moments>> sums{};
  std::vector<std::vector<double>> numberOfCounts{};
  double totalNumberOfCounts{0.0};

  bool isSampled(std::size_t component) const { return !backbonePerComponent[component].empty(); }
  std::size_t numberOfBackboneBeads(std::size_t component) const { return backbonePerComponent[component].size(); }

  // Contributions of a single molecule (exposed for testing). 'backbone' indexes into 'molecule'.
  static void accumulateInternalDistances(std::span<const Atom> molecule, std::span<const std::size_t> backbone,
                                          std::span<double> internalDistanceSquared);
  static void accumulateBondCorrelation(std::span<const Atom> molecule, std::span<const std::size_t> backbone,
                                        std::span<double> bondCorrelation);
  static void accumulateFormFactor(std::span<const Atom> molecule, std::span<const double> waveVectors,
                                   std::span<double> formFactor);
  static Moments computeMoments(std::span<const Atom> molecule, std::span<const std::size_t> backbone);

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
