module;

export module forcefield_settings;

import std;

import archive;
import uint3;
import pseudo_atom;
import json;

/**
 * \brief The sampling and numerical settings read from the force-field file, as opposed to the force-field
 * parameters themselves.
 *
 * A 'ForceField' is split in two: the parameters that define the Hamiltonian (pseudo-atoms, pair potentials,
 * cut-offs, truncation, charge method and Ewald parameters, external field) live directly in 'ForceField' and
 * are what every energy, gradient and Hessian routine reads; everything that only steers HOW a simulation
 * samples or evaluates that Hamiltonian lives here and is reached through 'ForceField::settings':
 *
 *  - the configurational-bias / recoil-growth sampling parameters (trial counts, Rosenbluth threshold, ring
 *    closure probabilities, the dual cut-off scheme), which the CBMC code copies into its 'GrowthSettings';
 *  - the framework interpolation grids (which pseudo-atoms, spacing or number of points, interpolation
 *    scheme, test points, output);
 *  - the interpolation grid of the external field (whether to use it, its size, output).
 *
 * None of these changes a single pair energy, so the pair kernels never read this struct. The keys are read
 * from the same JSON file as the parameters ('readFromJSON'), serialized with the force field and printed as
 * part of its status.
 */
export struct ForceFieldSettings
{
  enum class InterpolationGridType : std::size_t
  {
    LennardJones = 0,
    LennardJonesRepulsion = 1,
    LennardJonesAttraction = 2,
    EwaldReal = 3
  };

  enum class InterpolationScheme : std::size_t
  {
    Polynomial = 1,
    Tricubic = 8,
    Triquintic = 27
  };

  std::uint64_t versionNumber{1};  ///< Version number of the settings format.

  // Configurational-bias Monte Carlo sampling parameters. None of these changes the sampled distribution,
  // only the efficiency (and cost) of the sampling.
  std::size_t numberOfTrialDirections{10};         ///< Trial directions 'k' per configurational-bias step.
  std::size_t numberOfTorsionTrialDirections{100}; ///< Trial spins per trial direction for the torsion selection.
  std::size_t numberOfFirstBeadPositions{10};      ///< Trial positions of the first bead.
  std::size_t numberOfTrialMovesPerOpenBead{150};  ///< Internal Metropolis moves per placed bead (tilt, ring closure).
  double minimumRosenbluthFactor{1e-150};          ///< Minimum allowed Rosenbluth factor.

  // Internal ring-closure Monte-Carlo tuning (CBMC growth of cyclic clusters). Per internal-MC trial,
  // the conformer-hopping crankshaft is attempted with probability 'cbmcRingCrankshaftProbability';
  // otherwise a whole-ring junction tilt is attempted with probability 'cbmcRingTiltProbability' and a
  // local displacement/rotation with the remainder. These affect sampling efficiency only (the moves
  // carry no Rosenbluth weight), so any value in [0, 1] is valid.
  double cbmcRingCrankshaftProbability{0.2};  ///< Attempt probability of the large-angle ring crankshaft.
  double cbmcRingTiltProbability{0.25};       ///< Attempt probability of the whole-ring junction tilt.

  // Recoil-growth (RG) options for flexible molecules (Consta et al., Mol. Phys. 97, 1243 (1999)).
  // When 'useRecoilGrowth' is true, the flexible-molecule chain is grown/retraced with the recoil
  // growth algorithm instead of configurational-bias Monte Carlo (CBMC).
  bool useRecoilGrowth{false};                         ///< Use recoil growth instead of CBMC for flexible molecules.
  std::size_t recoilGrowthMaximumRecoilLength{2};      ///< Feeler / recoil length 'l' (look-ahead depth).
  std::size_t recoilGrowthNumberOfTrialDirections{5};  ///< Number of trial directions 'k' per segment in RG.

  // Dual cut-off scheme: the CBMC grow and retrace evaluate the external energies at the inner cut-off
  // 'dualCutOff' and the result is corrected to the full cut-offs afterwards (see CBMC::GrowContext).
  bool useDualCutOff{false};  ///< Indicates if the dual cut-off scheme is used.
  double dualCutOff{6.0};     ///< Inner cut-off distance when using the dual cut-off scheme.

  // Framework interpolation grids.
  std::vector<std::size_t> gridPseudoAtomIndices{};  ///< Pseudo-atoms whose framework interaction is tabulated.
  double spacingVDWGrid{0.15};
  double spacingCoulombGrid{0.15};
  std::optional<uint3> numberOfVDWGridPoints{};
  std::optional<uint3> numberOfCoulombGridPoints{};
  std::size_t numberOfGridTestPoints{100000};
  bool interpolationSchemeAuto{true};
  InterpolationScheme interpolationScheme{InterpolationScheme::Polynomial};
  bool writeFrameworkInterpolationGrids{false};

  // Interpolation grid of the external field.
  bool useExternalFieldGrid{true};
  uint3 numberOfExternalFieldGridPoints{8, 8, 8};
  bool writeExternalFieldInterpolationGrid{false};

  bool operator==(const ForceFieldSettings &other) const = default;

  /**
   * \brief Reads the settings from the parsed force-field file.
   *
   * \param parsed_data The parsed JSON of the force-field file.
   * \param pseudoAtoms The pseudo-atoms of the force field (to resolve the names in 'UseInterpolationGrids').
   * \param smallestExplicitCutOff The smallest of the explicitly given (non-automatic) full cut-offs, used to
   *        validate 'DualCutOff'; std::numeric_limits<double>::max() when every full cut-off is automatic.
   */
  void readFromJSON(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data,
                    const std::vector<PseudoAtom> &pseudoAtoms, double smallestExplicitCutOff);

  /// The CBMC sampling part of the force-field status.
  std::string printSamplingStatus() const;
  /// The interpolation-grid part of the force-field status (empty when no grids are used).
  std::string printInterpolationGridStatus() const;

  /// The keys of the force-field file read by this struct.
  static const std::set<std::string> options;
  static bool isOption(const std::string &key);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ForceFieldSettings &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ForceFieldSettings &s);
};
