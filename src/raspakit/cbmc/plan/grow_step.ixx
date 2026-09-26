module;

export module cbmc_grow_step;

import std;

import intra_molecular_potentials;
import chiral_center;
import bond_potential;
import bend_potential;
import cbmc_closure_guide;

// One operator of the deterministic growth plan: its topology (what is grown from where) and the
// derived, step-constant data every sampler of the operator engine reads. Building the topology is
// 'cbmc_growth_plan', filling the derived data 'cbmc_step_terms'; both run once per plan, at setup.
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
 *    several at once) are sampled with the exact bonded base sampler plus the coupled-decoupled
 *    torsion selection; a rigid-body fragment ('rigidBody') is hinged on its already-placed connecting
 *    atom, with the junction bends sampled by a rigid-body tilt Monte-Carlo and the spin about the
 *    junction bond biased by the crossing torsions in the same torsion selection.
 *  - CloseRing: 'nextBeads' are all remaining atoms of a cyclic cluster (the 2-edge-connected
 *    component of the fragment graph), grown together with ring-closure CBMC: the internal cluster
 *    conformation is sampled from its Boltzmann distribution by an internal Monte-Carlo (the closure
 *    bonds keep every ring closed -- simple, fused, and bridged rings alike), and the placement is
 *    biased by the junction-crossing terms through the torsion selection.
 *  - CloseBridge: the last bead of a regrown interior segment (fixed-endpoint regrowth). The single
 *    next bead is bonded to the anchor 'currentBead' AND to a second placed bead 'closureBead'; it is
 *    sampled exactly in bipolar coordinates about the anchor-closure axis (both bond lengths from
 *    their Boltzmann densities, the position on the resulting circle spun by the torsion selection),
 *    see cbmc_bridge_closure. The attach steps that lead up to it carry closure guides
 *    ('spin.guides') that steer their spin towards closable geometries.
 *
 * The plan is deterministic (it only depends on the fragment graph and the set of already-placed
 * beads), so grow and retrace generate exactly the same sequence, as required for detailed balance.
 *
 * The struct has two layers. The TOPOLOGY ('kind' .. 'intra') is set by 'buildGrowthPlan' and is what
 * a step IS. Everything after it is derived SAMPLER DATA: pure functions of the topology (and the
 * component's chiral centres), filled once by 'prepareStep' and grouped per consumer -- 'base' for the
 * exact flexible base sampler, 'spin' for the torsion-spin selection, 'rigidTilt' and 'ring' for the
 * two internal Monte-Carlo samplers -- so the operator engine only does lookups. The one exception is
 * 'base.couplingConstants', which depends on the temperature and is filled by
 * 'Component::prepareGrowthPlans'.
 */
struct GrowStep
{
  enum class Kind : std::size_t
  {
    PlaceSeedFragment = 0,  ///< No previous bead: the seed of the growth.
    AttachFragment = 1,     ///< Grown from a placed anchor with a bend/torsion reference.
    CloseRing = 2,          ///< A cyclic cluster grown with ring-closure CBMC.
    CloseBridge = 3,        ///< The closing bead of an interior segment, bonded to two placed beads.
  };

  // ----------------------------------------------------------------------------------------------
  // Topology.
  // ----------------------------------------------------------------------------------------------
  Kind kind{Kind::AttachFragment};
  std::optional<std::size_t> previousBead{};  ///< Bend/torsion reference (absent for a seed).
  std::size_t currentBead{};                  ///< The placed anchor the step grows from.
  std::vector<std::size_t> nextBeads{};       ///< The beads placed by this step.
  bool rigidBody{false};                      ///< Whether 'nextBeads' are hinged as one rigid body.
  Potentials::IntraMolecularPotentials intra{};  ///< Interactions affecting the placement of 'nextBeads'.
  /// CloseBridge only: the second placed bead the single next bead is bonded to.
  std::optional<std::size_t> closureBead{};

  // ----------------------------------------------------------------------------------------------
  // Derived sampler data.
  // ----------------------------------------------------------------------------------------------

  /// A flexible attach step (AttachFragment, not rigid, with a previous bead): its base sampler and
  /// torsion selection handle every internal term of the step themselves (see
  /// 'stepHandlesUnsampledInternalTerms').
  bool flexibleAttach{false};

  /// What the exact flexible base sampler reads (filled for every step with a previous bead; the
  /// coupling members only for a flexible attach step).
  struct BaseSamplerData
  {
    /// Per next bead (same order as 'nextBeads'): its bond to the current bead, and its anchor bend
    /// previous-current-next (centred on the current bead). Absent when the topology declares none.
    std::vector<std::optional<BondPotential>> bonds{};
    std::vector<std::optional<BendPotential>> anchorBends{};

    /// Sibling bends next_i - current - next_j of the step: spin-invariant, imposed on the base
    /// conformation by rejection.
    std::vector<BendPotential> siblingBends{};

    /// The unsampled internal terms of a flexible attach step that are fully determined by the step's
    /// own sampled coordinates (base coupling, imposed by rejection). 'hasCouplingTerms' is its
    /// non-emptiness; 'hasCoupling' also counts the sibling bends.
    Potentials::IntraMolecularPotentials couplingTerms{};
    bool hasCouplingTerms{false};
    bool hasCoupling{false};

    /// Declared chiral centres fully determined by this step (centred on the current bead, every
    /// neighbour the previous bead or grown here); enforced by parity rejection.
    std::vector<ChiralCenter> determinedChiralCenters{};

    /// Temperature-independent memo key of the base coupling: the step's per-bead samplers and every
    /// coupling term, with atom identifiers mapped to step-local roles so congruent steps share it.
    /// Empty when the step has no base coupling.
    std::string couplingSignature{};

    /// The base-coupling constants for the component's temperature ('Component::prepareGrowthPlans');
    /// absent until prepared, and for steps without base coupling.
    std::optional<BaseCouplingConstants> couplingConstants{};
  };
  BaseSamplerData base{};

  /// What the torsion-spin selection about the previous-current axis reads.
  struct SpinSelectionData
  {
    /// The potentials evaluated per spin trial: the step's torsions (for a ring-closure step only the
    /// junction-crossing ones) plus, for a flexible attach step, the spin-routed share of the
    /// unsampled terms. 'hasUnsampledTerms' tells whether that share is non-empty.
    Potentials::IntraMolecularPotentials potentials{};
    bool hasUnsampledTerms{false};

    /// Bends of the step that change under the spin (they involve a placed atom other than the
    /// previous bead); weighted alongside the torsions.
    std::vector<BendPotential> variantBends{};

    /// Closure guides of an attach step of a fixed-endpoint regrowth (see cbmc_closure_guide): per
    /// guided next bead, the placed bead it must eventually reach and the tabulated bias g(D) of
    /// their distance. The bias steers the spin selection and is divided out of its weight again.
    /// The table depends on the temperature and is filled by 'Component::prepareGrowthPlans'.
    struct ClosureGuide
    {
      std::size_t nextBeadIndex{};  ///< Index into 'nextBeads' of the guided bead.
      std::size_t targetBead{};     ///< The placed bead the guided bead's segment closes onto.
      ClosureGuidePath path{};      ///< The bonded terms along the shortest path bead -> target.
      std::string signature{};      ///< Memo key of 'path' (congruent paths share one table).
      std::shared_ptr<const ClosureGuideTable> table{};  ///< Absent until prepared.
    };
    std::vector<ClosureGuide> guides{};
  };
  SpinSelectionData spin{};

  /// CloseBridge steps: the two bonds of the closing bead (to the anchor and to the closure bead) and
  /// the bend centred on it (anchor - next - closure), which is invariant under the closure spin and
  /// enters the base weight; every other bonded term of the step is spin-variant and weighted in the
  /// torsion selection.
  struct BridgeClosureData
  {
    std::optional<BondPotential> anchorBond{};
    std::optional<BondPotential> closureBond{};
    std::optional<BendPotential> midBend{};
  };
  BridgeClosureData bridge{};

  /// Rigid-body steps with a junction: the body atom bonded to the anchor (the 'inner' atom of the
  /// junction bend previous-current-inner) and that bend when the topology declares one.
  struct RigidTiltData
  {
    std::size_t innerBead{};
    std::optional<BendPotential> junctionBend{};
  };
  RigidTiltData rigidTilt{};

  /// A conformer-hopping move of the ring-closure Monte-Carlo: 'atom' rotated about the line through
  /// its positioned neighbours 'axisA' and 'axisB'.
  struct RingCrankshaft
  {
    std::size_t atom{};
    std::size_t axisA{};
    std::size_t axisB{};
  };

  /// Ring-closure steps: the step-constant data of the internal conformational Monte-Carlo.
  struct RingSamplerData
  {
    /// The bonded terms the internal Monte-Carlo samples -- the cluster's bonds, and its spin-invariant
    /// bends and torsions (the spin-variant ones are weighted in the torsion selection instead, never
    /// double counted).
    Potentials::IntraMolecularPotentials internalPotentials{};
    /// Move units: a flexible ring atom alone, a rigid sub-fragment's atoms together.
    std::vector<std::vector<std::size_t>> moveUnits{};
    /// Per atom (indexed by atom id, sized to the molecule): its neighbours through Fixed bonds of the
    /// step, the pivots any move of that atom must rotate about.
    std::vector<std::vector<std::size_t>> fixedNeighbors{};
    /// Crankshaft candidates (flexible atoms with two positioned bonded neighbours, fixed ones first).
    std::vector<RingCrankshaft> crankshafts{};
    /// Declared chiral centres whose four atoms are all positioned during this step (ring body, anchor,
    /// junction neighbour); their parity is kept during the internal MC.
    std::vector<ChiralCenter> monitoredChiralCenters{};
  };
  RingSamplerData ring{};
};
}  // namespace CBMC
