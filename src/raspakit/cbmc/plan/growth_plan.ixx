module;

export module cbmc_growth_plan;

import std;

import connectivity_table;
import fragment_graph;
import intra_molecular_potentials;
export import cbmc_grow_step;

export namespace CBMC
{
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
 * carry no 'base.couplingConstants' yet (those need the temperature).
 */
std::vector<GrowStep> buildGrowthPlan(const ConnectivityTable &connectivity, const FragmentGraph &fragmentGraph,
                                      const Potentials::IntraMolecularPotentials &intraMolecularPotentials,
                                      const std::vector<std::size_t> &beadsAlreadyPlaced);

/**
 * \brief The steps of a plan that the operator engine can only grow with a built-in default geometry:
 * a flexible bead whose bond to its anchor has no bond potential (placed at
 * 'Constants::defaultBondLength'), and a hinged rigid body whose junction bend previous-current-inner
 * has no bend potential while the step carries other bends (tilted to
 * 'Constants::defaultRigidJunctionBendAngle').
 *
 * Such a topology is usually an omission in the component file (a 'Connectivity' entry without a
 * matching 'Bonds' term); the sampling stays exact for the potentials that ARE declared, but the
 * default geometry is not something the user chose. One human-readable message per affected step,
 * for the component reader to report as a warning.
 */
std::vector<std::string> growthPlanDefaultGeometryWarnings(const std::vector<GrowStep> &plan);
}  // namespace CBMC
