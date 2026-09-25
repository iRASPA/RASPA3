module;

export module cbmc_growth_plan;

import std;

import atom;
import connectivity_table;
import fragment_graph;
import intra_molecular_potentials;
import chiral_center;
import bond_potential;
import bend_potential;

export namespace CBMC
{
/**
 * \brief The frozen per-step-signature constants of the base coupling of a flexible attach step.
 *
 * Estimated once per step signature and temperature by a fixed-seed Monte Carlo that mirrors the base
 * sampler's own independent per-bead draws (see 'estimateBaseCouplingConstants'):
 *  - referenceEnergy u_ref: the offset of the rejection acceptance min(1, e^{-beta (u - u_ref)}). It
 *    does NOT need to be a rigorous bound: whenever u < u_ref the acceptance clamps at one and the
 *    excess e^{-beta (u - u_ref)} > 1 rides the trial's Rosenbluth weight instead, so the sampled
 *    density x weight product is exactly Boltzmann for any u_ref -- the constant only tunes efficiency.
 *  - logMeanClampedBoltzmann log<a>: the log mean CLAMPED acceptance over the independent base. The
 *    step's base normalization is Z_indep x <a> x e^{-beta u_ref}; the frozen relative error (< 0.5%)
 *    enters the reptation acceptance as a constant at the same sub-percent level.
 * Congruent steps (same signature) must share one estimate so their factors cancel exactly in the
 * reptation acceptance; 'Component::prepareGrowthPlans' guarantees this by memoising per signature.
 */
struct BaseCouplingConstants
{
  double referenceEnergy{0.0};
  double logMeanClampedBoltzmann{0.0};
};

/**
 * \brief One typed operator of the deterministic CBMC growth plan.
 *
 * The growth plan is a sequence of operators over the component's fragment graph. Each operator
 * places the beads 'nextBeads', growing from the placed anchor 'currentBead', optionally with a
 * 'previousBead' defining the bend/torsion reference:
 *  - PlaceSeedFragment: no previous bead exists (the growth seed). A single-atom fragment is placed
 *    with a Boltzmann bond length in a uniformly random direction; a rigid-body fragment
 *    ('rigidBody') is hinged on the anchor with a uniformly random orientation. A fully rigid
 *    molecule is one PlaceSeedFragment covering all atoms but the first bead.
 *  - AttachFragment: the unified attach operator. Single-atom fragments (a branch point can place
 *    several at once) are sampled with the bond/bend Monte-Carlo plus the coupled-decoupled torsion
 *    selection; a rigid-body fragment ('rigidBody') is hinged on its already-placed connecting atom,
 *    with the junction bends sampled by a rigid-body tilt Monte-Carlo and the spin about the
 *    junction bond biased by the crossing torsions in the same torsion selection.
 *  - CloseRing: 'nextBeads' are all remaining atoms of a cyclic cluster (the 2-edge-connected
 *    component of the fragment graph), grown together with ring-closure CBMC: the internal cluster
 *    conformation is sampled from its Boltzmann distribution by an internal Monte-Carlo (the closure
 *    bonds keep every ring closed -- simple, fused, and bridged rings alike), and the placement is
 *    biased by the junction-crossing terms through the torsion selection.
 *
 * The plan is deterministic (it only depends on the fragment graph and the set of already-placed
 * beads), so grow and retrace generate exactly the same sequence, as required for detailed balance.
 */
struct GrowStep
{
  enum class Kind : std::size_t
  {
    PlaceSeedFragment = 0,  ///< No previous bead: the seed of the growth.
    AttachFragment = 1,     ///< Grown from a placed anchor with a bend/torsion reference.
    CloseRing = 2,          ///< A cyclic cluster grown with ring-closure CBMC.
  };

  Kind kind{Kind::AttachFragment};
  std::optional<std::size_t> previousBead{};  ///< Bend/torsion reference (absent for a seed).
  std::size_t currentBead{};                  ///< The placed anchor the step grows from.
  std::vector<std::size_t> nextBeads{};       ///< The beads placed by this step.
  bool rigidBody{false};                      ///< Whether 'nextBeads' are hinged as one rigid body.
  Potentials::IntraMolecularPotentials intra{};  ///< Interactions affecting the placement of 'nextBeads'.

  // -------------------------------------------------------------------------------------------------
  // Derived per-step data, filled once by 'buildGrowthPlan' so the operator engine only does lookups.
  // Everything below is a pure function of the fields above (and the component's chiral centres),
  // except 'baseCouplingConstants', which depends on the temperature and is filled by
  // 'Component::prepareGrowthPlans'.
  // -------------------------------------------------------------------------------------------------

  /// A flexible attach step (AttachFragment, not rigid, with a previous bead): its base sampler and
  /// torsion selection handle every internal term of the step themselves (see
  /// 'stepHandlesUnsampledInternalTerms').
  bool flexibleAttach{false};

  /// Per next bead (same order as 'nextBeads'): its bond to the current bead, and its anchor bend
  /// previous-current-next (centred on the current bead). Absent when the topology declares none.
  /// Filled for every step with a previous bead; used by the flexible base sampler.
  std::vector<std::optional<BondPotential>> nextBeadBonds{};
  std::vector<std::optional<BendPotential>> nextBeadAnchorBends{};

  /// Sibling bends next_i - current - next_j of the step: spin-invariant, imposed on the base
  /// conformation by rejection.
  std::vector<BendPotential> siblingBends{};

  /// The unsampled internal terms of a flexible attach step that are fully determined by the step's
  /// own sampled coordinates (base coupling, imposed by rejection). 'hasBaseCouplingTerms' is its
  /// non-emptiness; 'hasBaseCoupling' also counts the sibling bends.
  Potentials::IntraMolecularPotentials baseCouplingTerms{};
  bool hasBaseCouplingTerms{false};
  bool hasBaseCoupling{false};

  /// Bends of the step that change under the spin about the previous-current axis (they involve a
  /// placed atom other than the previous bead); weighted in the torsion-spin selection.
  std::vector<BendPotential> spinVariantBends{};

  /// The potentials evaluated by the torsion-spin selection: the step's torsions (for a ring-closure
  /// step only the junction-crossing ones) plus, for a flexible attach step, the spin-routed share of
  /// the unsampled terms. 'torsionSelectionHasUnsampledTerms' tells whether that share is non-empty.
  Potentials::IntraMolecularPotentials torsionSelectionPotentials{};
  bool torsionSelectionHasUnsampledTerms{false};

  /// Ring-closure steps only: the bonded terms sampled by the internal conformational Monte-Carlo of
  /// the cyclic cluster -- its bonds, and the spin-invariant bends and torsions (the spin-variant
  /// ones are weighted in the torsion selection instead, never double counted).
  Potentials::IntraMolecularPotentials ringInternalPotentials{};

  /// Declared chiral centres fully determined by this step (centred on the current bead, every
  /// neighbour the previous bead or grown here); enforced by parity rejection in the base sampler.
  std::vector<ChiralCenter> determinedChiralCenters{};

  /// Temperature-independent memo key of the base coupling: the step's per-bead samplers and every
  /// coupling term, with atom identifiers mapped to step-local roles so congruent steps share it.
  /// Empty when the step has no base coupling.
  std::string baseCouplingSignature{};

  /// The base-coupling constants for the component's temperature ('Component::prepareGrowthPlans');
  /// absent until prepared, and for steps without base coupling.
  std::optional<BaseCouplingConstants> baseCouplingConstants{};
};

/**
 * \brief Builds the deterministic growth plan starting from 'beadsAlreadyPlaced'.
 *
 * The plan follows the fragment graph: fully flexible molecules reproduce the bead-by-bead order of
 * 'ConnectivityTable::nextBeads'; rigid-body fragments are emitted as a single hinged step (their
 * connecting atom is grown first as an ordinary flexible bead, so the junction bond and bend/torsion
 * are sampled); cyclic clusters are emitted as a single CloseRing step.
 *
 * Building a plan filters the intramolecular potentials per step and derives the per-step sampler
 * data (see the derived fields of 'GrowStep'), which is not cheap: simulation code should use the
 * cached plans through 'Component::growthPlan' instead of calling this directly. The returned steps
 * carry no 'baseCouplingConstants' yet (those need the temperature).
 */
std::vector<GrowStep> buildGrowthPlan(const ConnectivityTable &connectivity, const FragmentGraph &fragmentGraph,
                                      const Potentials::IntraMolecularPotentials &intraMolecularPotentials,
                                      const std::vector<std::size_t> &beadsAlreadyPlaced);
}  // namespace CBMC
