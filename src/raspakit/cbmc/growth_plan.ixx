module;

export module cbmc_growth_plan;

import std;

import atom;
import connectivity_table;
import fragment_graph;
import intra_molecular_potentials;

export namespace CBMC
{
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
};

/**
 * \brief Builds the deterministic growth plan starting from 'beadsAlreadyPlaced'.
 *
 * The plan follows the fragment graph: fully flexible molecules reproduce the bead-by-bead order of
 * 'ConnectivityTable::nextBeads'; rigid-body fragments are emitted as a single hinged step (their
 * connecting atom is grown first as an ordinary flexible bead, so the junction bond and bend/torsion
 * are sampled); cyclic clusters are emitted as a single CloseRing step.
 *
 * Building a plan filters the intramolecular potentials per step, which is not cheap: simulation
 * code should use the cached plans through 'Component::growthPlan' instead of calling this directly.
 */
std::vector<GrowStep> buildGrowthPlan(const ConnectivityTable &connectivity, const FragmentGraph &fragmentGraph,
                                      const Potentials::IntraMolecularPotentials &intraMolecularPotentials,
                                      const std::vector<std::size_t> &beadsAlreadyPlaced);
}  // namespace CBMC
