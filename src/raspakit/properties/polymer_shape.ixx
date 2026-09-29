module;

export module property_polymer_shape;

import std;

import archive;
import double3;
import double3x3;
import atom;
import forcefield;
import component;

// Samples the gyration-tensor family of single-chain shape descriptors for every component with at
// least two atoms.
//
// For each molecule the gyration tensor
//
//   S_ab = sum_i w_i (r_ia - r_cm,a) (r_ib - r_cm,b),     sum_i w_i = 1,
//
// is formed from the (unwrapped) atom positions with either uniform bead weights (the polymer-physics
// convention) or mass weights ('MassWeightedPolymerShape'), and diagonalised. With the eigenvalues
// ordered l1 >= l2 >= l3 the descriptors are
//
//   radius of gyration            Rg^2 = l1 + l2 + l3
//   asphericity                   b    = l1 - (l2 + l3) / 2                    (0 for a sphere)
//   acylindricity                 c    = l2 - l3                               (0 for a cylinder)
//   relative shape anisotropy     k^2  = (b^2 + 3 c^2 / 4) / Rg^4              (0 sphere ... 1 rod)
//   prolateness                   S    = 27 (l1 - lm)(l2 - lm)(l3 - lm) / Rg^6,  lm = Rg^2 / 3
//                                        (-1/4 oblate disk ... 0 sphere ... 2 prolate rod)
//
// Besides the per-molecule averages the ensemble-ratio forms are reported (Theodorou & Suter,
// Macromolecules 18, 1206 (1985)): <l1>:<l2>:<l3>, <b>/<Rg^2>, <c>/<Rg^2>, and
// k^2 = 1 - 3 <l1 l2 + l2 l3 + l3 l1> / <Rg^4>. The hydrodynamic radius follows in the Kirkwood
// (point-bead, pre-averaged) approximation 1/Rh = <(1/N^2) sum_{i != j} 1/r_ij>, and when the
// component has 'endToEndAtoms' the end-to-end distance is accumulated as well so the ratio
// <R^2>/<Rg^2> (6 for an ideal chain) is available in the same output.
//
// The lab-frame tensor itself is averaged too: <S_ab> (the "total" gyration tensor) is written with
// its eigenvalues and principal axes. In an isotropic bulk it tends to <Rg^2>/3 times the identity;
// its anisotropy measures how the chains are oriented by the simulation box or a framework (e.g.
// alignment along a channel direction). Note that the eigenvalues of <S> are not the <l_i> above:
// averaging the tensor before diagonalising mixes the molecular frames.
//
// Histograms of Rg, k^2 and S are accumulated per block; the Rg range defaults to 0.6 times the
// contour length of the bond-graph diameter (a bound on the Euclidean extent of the molecule) and can
// be overridden with 'RadiusOfGyrationRangePolymerShape'. Output is written to 'polymer_shape/'.
//
// For components that declare 'RepeatUnits' the same descriptors are also sampled per monomer: the
// gyration tensor of the atoms of each repeat unit (with the same weighting convention, renormalized
// within the unit), reported per unit index along the chain and pooled over all units. Units of
// fewer than two atoms are skipped.

export struct PropertyPolymerShape
{
  // Per-molecule moments accumulated per [block][component]; 'numberOfCounts' is their normalization.
  enum Moment : std::size_t
  {
    RadiusOfGyration = 0,         ///< Rg [Angstrom]
    RadiusOfGyrationSquared = 1,  ///< Rg^2 [Angstrom^2]
    RadiusOfGyrationFourth = 2,   ///< Rg^4 [Angstrom^4]
    Lambda1 = 3,                  ///< largest eigenvalue [Angstrom^2]
    Lambda2 = 4,                  ///< middle eigenvalue [Angstrom^2]
    Lambda3 = 5,                  ///< smallest eigenvalue [Angstrom^2]
    Asphericity = 6,              ///< b [Angstrom^2]
    Acylindricity = 7,            ///< c [Angstrom^2]
    ShapeAnisotropy = 8,          ///< k^2 [-]
    SecondInvariant = 9,          ///< l1 l2 + l2 l3 + l3 l1 [Angstrom^4]
    Prolateness = 10,             ///< S [-]
    InverseHydrodynamicRadius = 11,  ///< (1/N^2) sum_{i != j} 1/r_ij [1/Angstrom]
    EndToEndSquared = 12,         ///< R^2 [Angstrom^2] (zero when the component has no end-to-end atoms)
    // Lab-frame components of the gyration tensor S_ab [Angstrom^2]; their ensemble average <S_ab> is
    // the total tensor whose anisotropy reflects the orientation of the chains in the simulation box
    // (isotropic <S> = <Rg^2>/3 I in bulk; distinct diagonal entries for chains aligned by pores).
    TensorXX = 13,
    TensorYY = 14,
    TensorZZ = 15,
    TensorXY = 16,
    TensorXZ = 17,
    TensorYZ = 18,
    NumberOfMoments = 19
  };
  using Moments = std::array<double, NumberOfMoments>;

  // Instantaneous shape descriptors of a single molecule.
  struct Descriptors
  {
    double3x3 tensor{};     ///< The gyration tensor S_ab in the lab frame [Angstrom^2] (symmetric).
    double3 eigenvalues{};  ///< l1 >= l2 >= l3 of the gyration tensor [Angstrom^2].
    double radiusOfGyrationSquared{0.0};
    double asphericity{0.0};
    double acylindricity{0.0};
    double shapeAnisotropy{0.0};
    double prolateness{0.0};
  };

  PropertyPolymerShape() {};

  PropertyPolymerShape(std::size_t numberOfBlocks, const ForceField &forceField,
                       const std::vector<Component> &components, std::size_t numberOfBins, bool massWeighted,
                       std::size_t sampleEvery, std::optional<std::size_t> writeEvery,
                       std::optional<double> radiusOfGyrationRangeOverride = std::nullopt);

  std::uint64_t versionNumber{1};

  std::size_t numberOfBlocks{0};
  std::size_t numberOfBins{0};
  std::size_t numberOfComponents{0};
  bool massWeighted{false};

  std::size_t sampleEvery{10};
  std::optional<std::size_t> writeEvery{5000};

  // Per component: the normalized atom weights used for the center and the gyration tensor (empty for
  // components that are not sampled, i.e. with fewer than two atoms).
  std::vector<std::vector<double>> weightsPerComponent{};
  std::vector<std::optional<std::array<std::size_t, 2>>> endToEndAtomsPerComponent{};

  std::vector<double> radiusOfGyrationRangePerComponent{};  ///< Histogram upper limit [Angstrom].
  std::vector<double> deltaRadiusOfGyrationPerComponent{};  ///< Bin width [Angstrom].

  double shapeAnisotropyRange{1.0};      ///< k^2 in [0, 1].
  double prolatenessLowerLimit{-0.25};   ///< S in [-1/4, 2].
  double prolatenessRange{2.25};
  double deltaShapeAnisotropy{0.0};
  double deltaProlateness{0.0};

  // Histograms indexed as [block][component][bin].
  std::vector<std::vector<std::vector<double>>> radiusOfGyrationHistogram{};
  std::vector<std::vector<std::vector<double>>> shapeAnisotropyHistogram{};
  std::vector<std::vector<std::vector<double>>> prolatenessHistogram{};

  // Moment accumulators [block][component] and their normalization.
  std::vector<std::vector<Moments>> sums{};
  std::vector<std::vector<double>> numberOfCounts{};
  double totalNumberOfCounts{0.0};

  // Per-monomer sampling for components with repeat units: the atom indices of each unit (all units,
  // in chain order; units of fewer than two atoms are skipped when sampling), the normalized weights
  // per unit, and the accumulators [block][component][unit].
  std::vector<std::vector<std::vector<std::size_t>>> unitAtomsPerComponent{};
  std::vector<std::vector<std::vector<double>>> unitWeightsPerComponent{};
  std::vector<std::vector<std::vector<Moments>>> unitSums{};
  std::vector<std::vector<std::vector<double>>> unitCounts{};

  // Gyration-tensor descriptors of one molecule for the given normalized weights (one per atom).
  static Descriptors computeDescriptors(std::span<const Atom> molecule, std::span<const double> weights);

  // The same for the subset 'atoms' (indices into 'molecule'), with one normalized weight per
  // listed atom.
  static Descriptors computeDescriptors(std::span<const Atom> molecule, std::span<const std::size_t> atoms,
                                        std::span<const double> weights);

  // Accumulates the shape moments of a set of descriptors (the Kirkwood and end-to-end entries are
  // left to the caller).
  static void accumulateShapeMoments(Moments &moments, const Descriptors &descriptors);

  // Kirkwood sum (1/N^2) sum_{i != j} 1/r_ij of one molecule [1/Angstrom].
  static double computeInverseHydrodynamicRadius(std::span<const Atom> molecule);

  bool isSampled(std::size_t component) const { return !weightsPerComponent[component].empty(); }
  std::size_t numberOfUnits(std::size_t component) const { return unitAtomsPerComponent[component].size(); }

  void sample(const std::vector<Component> &components,
              const std::vector<std::size_t> &numberOfMoleculesPerComponent, std::span<const Atom> moleculeAtoms,
              std::size_t currentCycle, std::size_t block);

  // Block-combined mean and 95% confidence-interval error of a function of the per-molecule moment
  // averages: the mean is evaluated on the moments averaged over all blocks, the error from the scatter
  // of the function evaluated on each block's averages.
  std::pair<double, double> statistics(std::size_t component,
                                       const std::function<double(const Moments &)> &function) const;

  // Mean and error of one accumulated moment.
  std::pair<double, double> momentStatistics(std::size_t component, Moment moment) const;

  // The same for one repeat unit, or pooled over all units of the component (unit = nullopt).
  std::pair<double, double> unitStatistics(std::size_t component, std::optional<std::size_t> unit,
                                           const std::function<double(const Moments &)> &function) const;

  // Shared block-statistics kernel: 'sumOf(block)' and 'countOf(block)' provide the per-block
  // accumulators.
  std::pair<double, double> blockStatistics(const std::function<Moments(std::size_t)> &sumOf,
                                            const std::function<double(std::size_t)> &countOf,
                                            const std::function<double(const Moments &)> &function) const;

  // Returns {binCenters, averageProbabilityDensity, confidenceIntervalError} for one histogram.
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> result(
      const std::vector<std::vector<std::vector<double>>> &histogram, std::size_t component, double delta,
      double rangeStart) const;

  void writeOutput(std::size_t systemId, const std::vector<Component> &components, std::size_t currentCycle);

  std::string printSettings() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPolymerShape &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPolymerShape &p);
};
