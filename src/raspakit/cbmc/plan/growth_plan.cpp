module;

module cbmc_growth_plan;

import std;

import atom;
import connectivity_table;
import fragment;
import fragment_graph;
import intra_molecular_potentials;
import units;
import chiral_center;
import bond_potential;
import bend_potential;
import torsion_potential;
import cbmc_constants;
import cbmc_step_terms;
import cbmc_closure_guide;

// First placed neighbor (ascending atom index) of 'anchor' that satisfies 'acceptable'; used as the
// deterministic bend/torsion reference of a step.
template <typename Predicate>
static std::optional<std::size_t> firstPlacedNeighbor(const ConnectivityTable &connectivity, std::size_t anchor,
                                                      const std::vector<bool> &placed, Predicate acceptable)
{
  for (std::size_t i = 0; i != connectivity.numberOfBeads; ++i)
  {
    if (connectivity[i, anchor] && placed[i] && acceptable(i))
    {
      return i;
    }
  }
  return std::nullopt;
}

// The placed neighbours of 'bead' other than 'except' (ascending atom index).
static std::vector<std::size_t> placedNeighborsExcept(const ConnectivityTable &connectivity, std::size_t bead,
                                                      const std::vector<bool> &placed, std::size_t except)
{
  std::vector<std::size_t> result{};
  for (std::size_t i = 0; i != connectivity.numberOfBeads; ++i)
  {
    if (i != except && connectivity[i, bead] && placed[i]) result.push_back(i);
  }
  return result;
}

// Fixed-endpoint regrowth: a placed set with more than one connected piece leaves interior segments
// whose growth must eventually reach a second placed bead. The plan supports this for flexible
// single-atom beads (the closing bead gets a CloseBridge step, the beads before it closure guides). A
// rigid body or a cyclic cluster that is grown INTO a second placed bead has no exact sampler here and
// is rejected at plan-building time (i.e. when the component file is read), with the offending bond.
static void rejectExternalClosure(const ConnectivityTable &connectivity, const std::vector<bool> &placed,
                                  std::size_t currentBead, std::span<const std::size_t> unitAtoms,
                                  std::string_view unitName)
{
  std::vector<bool> inUnit(connectivity.numberOfBeads, false);
  for (std::size_t atom : unitAtoms) inUnit[atom] = true;
  for (std::size_t atom : unitAtoms)
  {
    if (placed[atom]) continue;
    for (std::size_t p : placedNeighborsExcept(connectivity, atom, placed, currentBead))
    {
      if (inUnit[p]) continue;
      throw std::runtime_error(std::format(
          "CBMC growth plan: atom {} of a {} is bonded to the placed atom {} while the {} is grown from atom {}. "
          "Fixed-endpoint regrowth (a partial reinsertion whose fixed atoms enclose the regrown part) is only "
          "supported for flexible single-atom beads; regrow the whole {} or keep it fixed.\n",
          atom, unitName, p, unitName, currentBead, unitName));
    }
  }
}

// Shortest path (through unplaced beads) from 'start' to a placed bead other than 'currentBead': the
// closure target of a guided bead and the atoms along the way (start first, target last). Breadth
// first, so the path is the shortest; ties are broken deterministically (lowest target index at the
// shortest distance). Empty when no such placed bead is reachable (ordinary tail growth).
static std::vector<std::size_t> shortestPathToPlacedBead(const ConnectivityTable &connectivity,
                                                         const std::vector<bool> &placed, std::size_t currentBead,
                                                         std::size_t start)
{
  const std::size_t numberOfBeads = connectivity.numberOfBeads;
  std::vector<std::optional<std::size_t>> parent(numberOfBeads);
  std::vector<bool> visited(numberOfBeads, false);
  std::vector<std::size_t> frontier{start};
  visited[start] = true;

  while (!frontier.empty())
  {
    // Every target adjacent to this level is at the same distance; take the lowest index.
    std::optional<std::size_t> target{};
    std::optional<std::size_t> targetParent{};
    std::vector<std::size_t> next{};
    for (std::size_t u : frontier)
    {
      for (std::size_t v = 0; v != numberOfBeads; ++v)
      {
        if (!connectivity[u, v]) continue;
        if (placed[v])
        {
          if (v != currentBead && (!target.has_value() || v < target.value()))
          {
            target = v;
            targetParent = u;
          }
          continue;
        }
        if (visited[v]) continue;
        visited[v] = true;
        parent[v] = u;
        next.push_back(v);
      }
    }
    if (target.has_value())
    {
      std::vector<std::size_t> path{target.value()};
      for (std::optional<std::size_t> u = targetParent; u.has_value(); u = parent[u.value()]) path.push_back(u.value());
      std::reverse(path.begin(), path.end());
      return path;
    }
    frontier = std::move(next);
  }
  return {};
}

// The bonded terms of the full potential along a path (an absent term is the engine's default).
static CBMC::ClosureGuidePath guidePathPotentials(const Potentials::IntraMolecularPotentials &intra,
                                                  std::span<const std::size_t> path)
{
  CBMC::ClosureGuidePath result{};
  for (std::size_t i = 0; i + 1 < path.size(); ++i)
  {
    result.bonds.push_back(intra.findBondPotential(path[i], path[i + 1]));
  }
  for (std::size_t i = 0; i + 2 < path.size(); ++i)
  {
    std::optional<BendPotential> found{};
    for (const BendPotential &bend : intra.bends)
    {
      const auto &ids = bend.identifiers;
      if (ids[1] != path[i + 1]) continue;
      if ((ids[0] == path[i] && ids[2] == path[i + 2]) || (ids[2] == path[i] && ids[0] == path[i + 2]))
      {
        found = bend;
        break;
      }
    }
    result.bends.push_back(found);
  }
  for (std::size_t i = 0; i + 3 < path.size(); ++i)
  {
    std::optional<TorsionPotential> found{};
    for (const TorsionPotential &torsion : intra.torsions)
    {
      const auto &ids = torsion.identifiers;
      const bool forward = ids[0] == path[i] && ids[1] == path[i + 1] && ids[2] == path[i + 2] && ids[3] == path[i + 3];
      const bool backward = ids[3] == path[i] && ids[2] == path[i + 1] && ids[1] == path[i + 2] && ids[0] == path[i + 3];
      if (forward || backward)
      {
        found = torsion;
        break;
      }
    }
    result.torsions.push_back(found);
  }
  return result;
}

// Attaches the closure guides of a flexible attach step: for every next bead from which a placed
// bead other than the anchor is reachable through unplaced beads (a bead of an interior segment), the
// shortest path to it and the bonded terms along that path. The tables themselves need the
// temperature and are filled later ('prepareClosureGuides').
static void attachClosureGuides(CBMC::GrowStep &step, const ConnectivityTable &connectivity,
                                const std::vector<bool> &placed,
                                const Potentials::IntraMolecularPotentials &intraMolecularPotentials)
{
  for (std::size_t k = 0; k != step.nextBeads.size(); ++k)
  {
    const std::vector<std::size_t> path =
        shortestPathToPlacedBead(connectivity, placed, step.currentBead, step.nextBeads[k]);
    if (path.size() < 3) continue;  // no target, or (size 2) a closure bead, which is never an attach bead

    CBMC::GrowStep::SpinSelectionData::ClosureGuide guide{};
    guide.nextBeadIndex = k;
    guide.targetBead = path.back();
    guide.path = guidePathPotentials(intraMolecularPotentials, path);
    guide.signature = guide.path.signature();
    step.spin.guides.push_back(std::move(guide));
  }
}

// Determine the next growth step given the set of already-placed beads.
//
// Mirrors 'ConnectivityTable::nextBeads' for the fully flexible case, but rigid-body fragments are
// hinged as one rigid body once their connecting atom is placed, and cyclic clusters (all rings:
// simple, fused, bridged) are grown as one ring-closure step.
static CBMC::GrowStep nextGrowthStep(const ConnectivityTable &connectivity, const FragmentGraph &graph,
                                     const Potentials::IntraMolecularPotentials &intraMolecularPotentials,
                                     const std::vector<std::size_t> &placedBeads)
{
  std::size_t numberOfBeads = connectivity.numberOfBeads;

  std::vector<bool> placed(numberOfBeads, false);
  for (std::size_t bead : placedBeads) placed[bead] = true;

  // 'filteredInteractions' takes mutable spans, so keep a mutable copy of the placed beads.
  std::vector<std::size_t> placedVec(placedBeads.begin(), placedBeads.end());

  auto makeStep = [&](CBMC::GrowStep::Kind kind, std::optional<std::size_t> previousBead, std::size_t currentBead,
                      std::vector<std::size_t> nextBeads, bool rigidBody) -> CBMC::GrowStep
  {
    Potentials::IntraMolecularPotentials intra =
        intraMolecularPotentials.filteredInteractions(numberOfBeads, placedVec, nextBeads);
    return {kind, previousBead, currentBead, std::move(nextBeads), rigidBody, std::move(intra)};
  };

  // A partially placed rigid-body fragment is completed first: the anchor is its placed atom (the
  // connecting atom grown as an ordinary flexible bead in the previous step, or a bead of
  // 'beadsAlreadyPlaced'), and the remaining atoms are hinged on it as one rigid body. This rule
  // also covers a fully rigid molecule (no frontier bond exists when there is no connectivity).
  for (std::size_t anchor : placedBeads)
  {
    std::size_t fragmentIndex = graph.atomFragmentIds[anchor];
    const Fragment &fragment = graph.fragments[fragmentIndex];
    if (!fragment.isRigidBody()) continue;

    std::vector<std::size_t> nextBeads{};
    nextBeads.reserve(fragment.atoms.size());
    for (std::size_t atom : fragment.atoms)
    {
      if (!placed[atom]) nextBeads.push_back(atom);
    }
    if (nextBeads.empty()) continue;

    rejectExternalClosure(connectivity, placed, anchor, fragment.atoms, "rigid body");

    // Reference bead for the junction bend/torsion bias: a placed neighbor of the anchor outside the
    // fragment. Absent only when the fragment is the growth seed.
    std::optional<std::size_t> previousBead = firstPlacedNeighbor(
        connectivity, anchor, placed, [&](std::size_t i) { return graph.atomFragmentIds[i] != fragmentIndex; });

    CBMC::GrowStep::Kind kind =
        previousBead.has_value() ? CBMC::GrowStep::Kind::AttachFragment : CBMC::GrowStep::Kind::PlaceSeedFragment;
    return makeStep(kind, previousBead, anchor, std::move(nextBeads), true);
  }

  // Search for the first frontier bond (a placed bead 'k' connected to an unplaced bead 'j'),
  // scanning placed beads in order and neighbors by ascending index (identical ordering to
  // 'ConnectivityTable::nextBeads' so the fully flexible case is unchanged).
  std::optional<std::size_t> currentBeadOpt{};
  std::optional<std::size_t> firstNextBead{};
  for (std::size_t k : placedBeads)
  {
    for (std::size_t j = 0; j != numberOfBeads; ++j)
    {
      if (connectivity[j, k] && !placed[j])
      {
        currentBeadOpt = k;
        firstNextBead = j;
        break;
      }
    }
    if (currentBeadOpt.has_value()) break;
  }

  if (!currentBeadOpt.has_value())
  {
    throw std::runtime_error(std::format("Error in CBMC: No bead can be grown\n"));
  }

  std::size_t currentBead = currentBeadOpt.value();
  std::size_t nextBead = firstNextBead.value();

  auto clusterOf = [&](std::size_t bead) { return graph.fragmentCyclicClusterIds[graph.atomFragmentIds[bead]]; };

  // The frontier bond leads into a cyclic cluster.
  if (std::optional<std::size_t> cluster = clusterOf(nextBead); cluster.has_value())
  {
    rejectExternalClosure(connectivity, placed, currentBead, graph.cyclicClusters[*cluster], "flexible ring");

    if (clusterOf(currentBead) == cluster)
    {
      // The anchor is a placed atom of the cluster: grow the remaining cluster atoms as one
      // ring-closure step (the closure bonds in 'intra' keep every ring of the cluster closed).
      const std::vector<std::size_t> &clusterAtoms = graph.cyclicClusters[*cluster];
      std::vector<std::size_t> nextBeads{};
      nextBeads.reserve(clusterAtoms.size());
      for (std::size_t atom : clusterAtoms)
      {
        if (!placed[atom]) nextBeads.push_back(atom);
      }

      // Reference bead for the junction bend/torsion bias: a placed neighbor of the anchor outside
      // the cluster. Absent when the cluster is the growth seed.
      std::optional<std::size_t> previousBead = firstPlacedNeighbor(connectivity, currentBead, placed,
                                                                    [&](std::size_t i) { return clusterOf(i) != cluster; });

      return makeStep(CBMC::GrowStep::Kind::CloseRing, previousBead, currentBead, std::move(nextBeads), false);
    }

    // The frontier bond leads from outside into the cluster. Grow the single connecting atom first
    // as an ordinary flexible bead (its junction bond, and its bend/torsion when a previous bead
    // exists, are sampled); once placed it becomes the anchor of the ring-closure step.
    std::optional<std::size_t> previousBead =
        firstPlacedNeighbor(connectivity, currentBead, placed, [](std::size_t) { return true; });
    CBMC::GrowStep::Kind kind =
        previousBead.has_value() ? CBMC::GrowStep::Kind::AttachFragment : CBMC::GrowStep::Kind::PlaceSeedFragment;
    return makeStep(kind, previousBead, currentBead, {nextBead}, false);
  }

  // The frontier bond leads into a rigid-body fragment (the anchor is outside it, otherwise the
  // completion rule above would have fired). Grow the single connecting atom first as an ordinary
  // flexible bead; the remainder of the fragment is hinged on it in the next step.
  if (graph.fragments[graph.atomFragmentIds[nextBead]].isRigidBody())
  {
    rejectExternalClosure(connectivity, placed, currentBead, graph.fragments[graph.atomFragmentIds[nextBead]].atoms,
                          "rigid body");

    std::optional<std::size_t> previousBead =
        firstPlacedNeighbor(connectivity, currentBead, placed, [](std::size_t) { return true; });
    CBMC::GrowStep::Kind kind =
        previousBead.has_value() ? CBMC::GrowStep::Kind::AttachFragment : CBMC::GrowStep::Kind::PlaceSeedFragment;
    return makeStep(kind, previousBead, currentBead, {nextBead}, false);
  }

  // Fixed-endpoint regrowth: a flexible frontier bead that is ALSO bonded to a second placed bead is the
  // closing bead of an interior segment (a placed set with several connected pieces; in a tree
  // topology with a connected placed set this never happens). It is grown alone by a CloseBridge step
  // that samples it exactly on the circle spanned by its two bonds (cbmc_bridge_closure). A bead bonded
  // to three or more placed beads is over-determined (no single-bead sampler exists) and is rejected.
  {
    const std::vector<std::size_t> closurePartners = placedNeighborsExcept(connectivity, nextBead, placed, currentBead);
    if (closurePartners.size() >= 2)
    {
      throw std::runtime_error(std::format(
          "CBMC growth plan: atom {} is bonded to {} placed atoms ({}, {}, ...) while it is grown from atom {}. "
          "Fixed-endpoint regrowth closes a segment through one bond; a bead bonded to three or more fixed atoms "
          "cannot be regrown on its own. Regrow at least one of its neighbours as well.\n",
          nextBead, closurePartners.size() + 1, currentBead, closurePartners[0], currentBead));
    }
    if (closurePartners.size() == 1)
    {
      const std::size_t closureBead = closurePartners.front();
      for (const ChiralCenter &center : intraMolecularPotentials.chiralCenters)
      {
        if (std::ranges::find(center.ids, nextBead) != center.ids.end())
        {
          throw std::runtime_error(std::format(
              "CBMC growth plan: atom {} closes a regrown segment (bonded to the fixed atoms {} and {}) but takes part "
              "in the declared chiral centre [{}, {}, {}, {}]; the closure sampler does not enforce chirality. Choose "
              "a fixed set whose closing bead is not part of a chiral centre.\n",
              nextBead, currentBead, closureBead, center.ids[0], center.ids[1], center.ids[2], center.ids[3]));
        }
      }
      std::optional<std::size_t> previousBead =
          firstPlacedNeighbor(connectivity, currentBead, placed, [](std::size_t) { return true; });
      CBMC::GrowStep step = makeStep(CBMC::GrowStep::Kind::CloseBridge, previousBead, currentBead, {nextBead}, false);
      step.closureBead = closureBead;
      return step;
    }
  }

  // Flexible step. Reproduce 'ConnectivityTable::nextBeads': determine the previous bead and, when
  // there is a previous bead, grow all unplaced flexible neighbors of 'currentBead' together (branch
  // point); when there is no previous bead grow only the single frontier bead.
  std::vector<std::size_t> nextBeads{};
  std::size_t numberOfPreviousBeads{};
  std::optional<std::size_t> previousBead{};
  for (std::size_t i = 0; i != numberOfBeads; ++i)
  {
    if (!connectivity[i, currentBead]) continue;

    if (!placed[i])
    {
      // Rigid-body or cyclic-cluster neighbors are grown as their own step, never mixed into a
      // flexible branch; nor is a closing bead (it gets its own CloseBridge step once the branch is
      // placed).
      if (!graph.fragments[graph.atomFragmentIds[i]].isRigidBody() && !clusterOf(i).has_value() &&
          placedNeighborsExcept(connectivity, i, placed, currentBead).empty())
      {
        nextBeads.push_back(i);
      }
    }
    else
    {
      ++numberOfPreviousBeads;
      // Keep the first (lowest-index) placed neighbor as the deterministic bend/torsion reference.
      if (!previousBead.has_value()) previousBead = i;
    }
  }

  if (numberOfPreviousBeads == 0)
  {
    // A seed grown from a fixed bead without placed neighbours (e.g. a fixed chain end of a
    // fixed-endpoint regrowth) still gets its closure guides: the seed operator applies them to the
    // direction (there is no junction to spin about).
    CBMC::GrowStep step = makeStep(CBMC::GrowStep::Kind::PlaceSeedFragment, std::nullopt, currentBead, {nextBead}, false);
    attachClosureGuides(step, connectivity, placed, intraMolecularPotentials);
    return step;
  }

  // 'currentBead' with more than one placed neighbor is expected when it belongs to an already-placed
  // rigid fragment or cyclic cluster (e.g. a flexible tail growing off a ring atom). It can also
  // happen at a purely flexible acyclic branch point. Both cases have a well-defined growth order via
  // the spanning tree of 'FragmentGraph', and grow and retrace build from the same deterministic
  // plan, so the lowest-index placed neighbour is a consistent bend/torsion reference that preserves
  // detailed balance. The former code rejected the flexible acyclic case with "Multiple previous
  // beads"; the spanning tree removes the need for that restriction, so it is simply allowed.
  CBMC::GrowStep step =
      makeStep(CBMC::GrowStep::Kind::AttachFragment, previousBead, currentBead, std::move(nextBeads), false);
  attachClosureGuides(step, connectivity, placed, intraMolecularPotentials);
  return step;
}

std::vector<CBMC::GrowStep> CBMC::buildGrowthPlan(
    const ConnectivityTable &connectivity, const FragmentGraph &fragmentGraph,
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials,
    const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  std::size_t numberOfBeads = connectivity.numberOfBeads;

  std::vector<GrowStep> plan{};
  std::vector<std::size_t> placed(beadsAlreadyPlaced.begin(), beadsAlreadyPlaced.end());

  while (placed.size() < numberOfBeads)
  {
    GrowStep step = nextGrowthStep(connectivity, fragmentGraph, intraMolecularPotentials, placed);
    placed.insert(placed.end(), step.nextBeads.begin(), step.nextBeads.end());
    CBMC::prepareStep(step, connectivity, fragmentGraph, intraMolecularPotentials.chiralCenters);
    plan.push_back(std::move(step));
  }

  return plan;
}

std::vector<std::string> CBMC::growthPlanDefaultGeometryWarnings(const std::vector<GrowStep> &plan)
{
  std::vector<std::string> warnings{};

  for (const GrowStep &step : plan)
  {
    if (step.kind == GrowStep::Kind::CloseRing) continue;

    // Bridge closure: both bonds of the closing bead are sampled.
    if (step.kind == GrowStep::Kind::CloseBridge)
    {
      if (!step.bridge.anchorBond.has_value())
      {
        warnings.push_back(std::format(
            "bond {}-{} has no bond potential; CBMC closes atom {} at the default length of {:g} Angstrom from atom {}. "
            "Declare the bond in 'Bonds'.",
            step.currentBead, step.nextBeads[0], step.nextBeads[0], Constants::defaultBondLength, step.currentBead));
      }
      if (!step.bridge.closureBond.has_value())
      {
        warnings.push_back(std::format(
            "bond {}-{} has no bond potential; CBMC closes atom {} at the default length of {:g} Angstrom from atom {}. "
            "Declare the bond in 'Bonds'.",
            step.nextBeads[0], step.closureBead.value(), step.nextBeads[0], Constants::defaultBondLength,
            step.closureBead.value()));
      }
      continue;
    }

    if (step.rigidBody)
    {
      // Mirrors the rigid tilt: a junction with bends but without the previous-current-inner bend
      // draws the tilt at the default angle.
      if (step.previousBead.has_value() && !step.intra.bends.empty() && !step.rigidTilt.junctionBend.has_value())
      {
        warnings.push_back(std::format(
            "rigid body hinged on atom {} has no bend potential {}-{}-{} for its junction; CBMC tilts it to the "
            "default angle of {:g} degrees. Declare the junction bend in 'Bends'.",
            step.currentBead, step.previousBead.value(), step.currentBead, step.rigidTilt.innerBead,
            Constants::defaultRigidJunctionBendAngle / Units::DegreesToRadians));
      }
      continue;
    }

    // Flexible seed: the single bond current-next is the front of the step's bonds when declared.
    if (!step.previousBead.has_value())
    {
      if (step.intra.bonds.empty())
      {
        warnings.push_back(std::format(
            "bond {}-{} has no bond potential; CBMC places atom {} at the default length of {:g} Angstrom. "
            "Declare the bond in 'Bonds'.",
            step.currentBead, step.nextBeads[0], step.nextBeads[0], Constants::defaultBondLength));
      }
      continue;
    }

    // Flexible attach: one bond per grown bead.
    for (std::size_t i = 0; i != step.nextBeads.size(); ++i)
    {
      if (step.base.bonds[i].has_value()) continue;
      warnings.push_back(std::format(
          "bond {}-{} has no bond potential; CBMC places atom {} at the default length of {:g} Angstrom. Declare "
          "the bond in 'Bonds'.",
          step.currentBead, step.nextBeads[i], step.nextBeads[i], Constants::defaultBondLength));
    }
  }

  return warnings;
}
