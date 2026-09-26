module;

export module cbmc_grow_context;

import atom;
import forcefield;
import framework;
import simulationbox;
import interpolation_energy_grid;

import std;

export namespace CBMC
{
/**
 * \brief Which cut-offs a growth context evaluates the external energies with.
 *
 *  - Growth: the cut-offs the CBMC grow/retrace itself uses -- the inner 'dualCutOff' for all three
 *    when the force field enables the dual cut-off scheme, the full cut-offs otherwise. The dual
 *    cut-off decision is made here, once; a caller that grows with this mode must apply
 *    'computeDualCutOffCorrection' afterwards when 'forceField.useDualCutOff' is set.
 *  - Full: the force field's full cut-offs regardless of the dual cut-off setting (e.g. the initial
 *    configuration and the reference grows, which never apply a correction).
 *  - Inner: the inner 'dualCutOff' for all three (the other side of the dual cut-off correction).
 */
enum class CutOffMode : std::size_t
{
  Growth = 0,
  Full = 1,
  Inner = 2,
};

/**
 * \brief Which scheme grows and retraces the chain beyond the first bead.
 *
 *  - ConfigurationalBias: the Rosenbluth scheme; the per-step weight is the sum of the Boltzmann
 *    factors of the trial directions. Its weight is the classic Rosenbluth weight, so its average over
 *    fresh insertions is the Widom estimator of the excess chemical potential.
 *  - RecoilGrowth: the recoil-growth scheme (Consta et al. 1999); the per-step weight counts the
 *    available (open, feeler-viable) directions and divides by the openness probability. Its weight
 *    is only established as a valid factor in a Metropolis acceptance RATIO; it is not the Rosenbluth
 *    weight whose average is the Widom estimator. Widom sampling therefore always grows with
 *    configurational bias, whatever the production moves use.
 *
 * The force field's 'useRecoilGrowth' sets the default of a context; a caller that needs a specific
 * scheme derives it with 'withChainScheme'.
 */
enum class ChainScheme : std::size_t
{
  ConfigurationalBias = 0,
  RecoilGrowth = 1,
};

/**
 * \brief The sampling parameters of a CBMC grow or retrace.
 *
 * None of these changes the sampled distribution, only the efficiency (and cost) of the sampling.
 * They are read from the force field file and copied into every 'GrowContext' at construction
 * ('fromForceField'); a caller that needs different parameters for one grow (a test, the Widom
 * estimator, the ideal-gas reference grows) derives a context with 'GrowContext::withSettings' or
 * 'withChainScheme' instead of editing the force field.
 */
struct GrowthSettings
{
  /// Chain scheme beyond the first bead.
  ChainScheme chainScheme{ChainScheme::ConfigurationalBias};
  /// Trial positions of the first bead ('FirstBeadScheme::MultipleFirstBead' and 'Reinsertion').
  std::size_t numberOfFirstBeadPositions{10};
  /// Trial directions 'k' per configurational-bias step.
  std::size_t numberOfTrialDirections{10};
  /// Trial spins per trial direction among which the torsion orientation is Rosenbluth-selected.
  std::size_t numberOfTorsionTrialDirections{100};
  /// Internal Metropolis moves per placed bead of the rigid-body tilt and ring-closure samplers.
  std::size_t numberOfTrialMovesPerOpenBead{150};
  /// Attempt probability of the large-angle ring crankshaft per internal ring move.
  double ringCrankshaftProbability{0.2};
  /// Attempt probability of the whole-ring junction tilt per internal ring move.
  double ringTiltProbability{0.25};
  /// A grow whose weight (of any step) falls below this is reported as failed (std::nullopt).
  double minimumRosenbluthFactor{1e-150};
  /// Trial directions 'k' per recoil-growth step.
  std::size_t recoilGrowthNumberOfTrialDirections{5};
  /// Recoil (feeler) length 'l' of recoil growth.
  std::size_t recoilGrowthMaximumRecoilLength{2};

  [[nodiscard]] static GrowthSettings fromForceField(const ForceField &forceField)
  {
    return GrowthSettings{
        .chainScheme = forceField.useRecoilGrowth ? ChainScheme::RecoilGrowth : ChainScheme::ConfigurationalBias,
        .numberOfFirstBeadPositions = forceField.numberOfFirstBeadPositions,
        .numberOfTrialDirections = forceField.numberOfTrialDirections,
        .numberOfTorsionTrialDirections = forceField.numberOfTorsionTrialDirections,
        .numberOfTrialMovesPerOpenBead = forceField.numberOfTrialMovesPerOpenBead,
        .ringCrankshaftProbability = forceField.cbmcRingCrankshaftProbability,
        .ringTiltProbability = forceField.cbmcRingTiltProbability,
        .minimumRosenbluthFactor = forceField.minimumRosenbluthFactor,
        .recoilGrowthNumberOfTrialDirections = std::max<std::size_t>(1, forceField.recoilGrowthNumberOfTrialDirections),
        .recoilGrowthMaximumRecoilLength = std::max<std::size_t>(1, forceField.recoilGrowthMaximumRecoilLength)};
  }
};

/**
 * \brief Everything a CBMC grow or retrace needs to know about its environment: the force field,
 * the box, the framework and its atoms, the background molecule atoms, the inverse temperature, the
 * cut-offs to evaluate external energies with, and the sampling parameters ('GrowthSettings').
 *
 * Holds references and spans into the owning system, so it is a cheap value and must not outlive it.
 * Build one with 'System::makeGrowContext' (or the constructor for an environment that is not a
 * system, e.g. the ideal-gas grows), and derive variants with the 'with...' members: a different
 * background ('withMoleculeAtoms', used by moves that grow against an edited copy of the molecule
 * atoms; 'withSkippedMolecule', used by moves whose trial molecule carries a different id than the
 * molecule it replaces), different cut-offs ('withCutOffs', 'withFullCutOffs', 'withInnerCutOffs',
 * used by the dual cut-off correction), or different sampling parameters ('withSettings',
 * 'withChainScheme'; Widom sampling and the ideal-gas reference grows must use configurational
 * bias). Those are the only fields a caller ever varies; every other field is fixed by the system,
 * which is why there is no aggregate initialization to keep in sync.
 *
 * Background exclusion: the inter-molecular energy never pairs a trial atom with a background atom
 * of the same molecule id, so a molecule regrown under its own id (reinsertion, partial
 * reinsertion, reptation) is automatically excluded from its own background. Only when the trial
 * atoms carry a NEW id while the old molecule is still in the background (identity change) must
 * the old molecule be excluded explicitly, with 'withSkippedMolecule'. The exclusion is part of the
 * environment, so every evaluation through the context (grow, retrace, dual cut-off correction)
 * applies it consistently.
 */
struct GrowContext
{
  GrowContext(bool hasExternalField, const ForceField &forceField, const SimulationBox &simulationBox,
              const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
              const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid,
              const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
              std::span<const Atom> moleculeAtoms, double beta, CutOffMode cutOffMode = CutOffMode::Growth)
      : hasExternalField(hasExternalField),
        forceField(forceField),
        simulationBox(simulationBox),
        interpolationGrids(interpolationGrids),
        externalFieldInterpolationGrid(externalFieldInterpolationGrid),
        framework(framework),
        frameworkAtoms(frameworkAtoms),
        moleculeAtoms(moleculeAtoms),
        beta(beta),
        cutOffFrameworkVDW(frameworkVDWCutOff(forceField, cutOffMode)),
        cutOffMoleculeVDW(moleculeVDWCutOff(forceField, cutOffMode)),
        cutOffCoulomb(coulombCutOff(forceField, cutOffMode)),
        settings(GrowthSettings::fromForceField(forceField))
  {
  }

  bool hasExternalField;
  const ForceField &forceField;
  const SimulationBox &simulationBox;
  const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids;
  const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid;
  const std::optional<Framework> &framework;
  std::span<const Atom> frameworkAtoms;
  std::span<const Atom> moleculeAtoms;
  /// Molecule id whose atoms in 'moleculeAtoms' are ignored by every energy evaluation (see the
  /// class comment); std::nullopt skips nothing beyond the same-id rule.
  std::optional<std::size_t> skipBackgroundMolecule{};
  double beta;
  double cutOffFrameworkVDW;
  double cutOffMoleculeVDW;
  double cutOffCoulomb;
  GrowthSettings settings;

  /// The same environment with the molecule background replaced (e.g. the system's molecule atoms
  /// with a pair removed, or with already grown group members appended).
  [[nodiscard]] GrowContext withMoleculeAtoms(std::span<const Atom> background) const
  {
    GrowContext copy(*this);
    copy.moleculeAtoms = background;
    return copy;
  }

  /// The same environment with the atoms of background molecule 'moleculeId' ignored (identity
  /// change: the trial molecule has a new id, the old molecule is still in the background).
  [[nodiscard]] GrowContext withSkippedMolecule(std::size_t moleculeId) const
  {
    GrowContext copy(*this);
    copy.skipBackgroundMolecule = moleculeId;
    return copy;
  }

  /// The same environment grown with the given sampling parameters.
  [[nodiscard]] GrowContext withSettings(const GrowthSettings &growthSettings) const
  {
    GrowContext copy(*this);
    copy.settings = growthSettings;
    return copy;
  }

  /// The same environment grown with the given chain scheme (see 'ChainScheme').
  [[nodiscard]] GrowContext withChainScheme(ChainScheme scheme) const
  {
    GrowContext copy(*this);
    copy.settings.chainScheme = scheme;
    return copy;
  }

  /// The same environment evaluated with the given cut-offs.
  [[nodiscard]] GrowContext withCutOffs(double frameworkVDW, double moleculeVDW, double coulomb) const
  {
    GrowContext copy(*this);
    copy.cutOffFrameworkVDW = frameworkVDW;
    copy.cutOffMoleculeVDW = moleculeVDW;
    copy.cutOffCoulomb = coulomb;
    return copy;
  }

  /// The same environment evaluated with the cut-offs of 'mode'.
  [[nodiscard]] GrowContext withCutOffMode(CutOffMode mode) const
  {
    return withCutOffs(frameworkVDWCutOff(forceField, mode), moleculeVDWCutOff(forceField, mode),
                       coulombCutOff(forceField, mode));
  }
  [[nodiscard]] GrowContext withFullCutOffs() const { return withCutOffMode(CutOffMode::Full); }
  [[nodiscard]] GrowContext withInnerCutOffs() const { return withCutOffMode(CutOffMode::Inner); }

 private:
  static bool usesInner(const ForceField &forceField, CutOffMode mode)
  {
    return mode == CutOffMode::Inner || (mode == CutOffMode::Growth && forceField.useDualCutOff);
  }
  static double frameworkVDWCutOff(const ForceField &forceField, CutOffMode mode)
  {
    return usesInner(forceField, mode) ? forceField.dualCutOff : forceField.cutOffFrameworkVDW;
  }
  static double moleculeVDWCutOff(const ForceField &forceField, CutOffMode mode)
  {
    return usesInner(forceField, mode) ? forceField.dualCutOff : forceField.cutOffMoleculeVDW;
  }
  static double coulombCutOff(const ForceField &forceField, CutOffMode mode)
  {
    return usesInner(forceField, mode) ? forceField.dualCutOff : forceField.cutOffCoulomb;
  }
};
}  // namespace CBMC
