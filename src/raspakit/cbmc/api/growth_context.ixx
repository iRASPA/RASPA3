module;

export module cbmc_growth_context;

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
 * \brief Everything a CBMC grow or retrace needs to know about its environment: the force field,
 * the box, the framework and its atoms, the background molecule atoms, the inverse temperature, and
 * the cut-offs to evaluate external energies with.
 *
 * Holds references and spans into the owning system, so it is a cheap value and must not outlive it.
 * Build one with 'System::makeGrowContext' (or the constructor for an environment that is not a
 * system, e.g. the ideal-gas grows), and derive variants with the 'with...' members: a different
 * background ('withMoleculeAtoms', used by moves that grow against an edited copy of the molecule
 * atoms) or different cut-offs ('withCutOffs', 'withFullCutOffs', 'withInnerCutOffs', used by the
 * dual cut-off correction). The cut-offs are the only fields a caller ever varies; every other field
 * is fixed by the system, which is why there is no aggregate initialization to keep in sync.
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
        cutOffCoulomb(coulombCutOff(forceField, cutOffMode))
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
  double beta;
  double cutOffFrameworkVDW;
  double cutOffMoleculeVDW;
  double cutOffCoulomb;

  /// The same environment with the molecule background replaced (e.g. the system's molecule atoms
  /// with a pair removed, or with already grown group members appended).
  [[nodiscard]] GrowContext withMoleculeAtoms(std::span<const Atom> background) const
  {
    GrowContext copy(*this);
    copy.moleculeAtoms = background;
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
