module;

module cbmc_growth_plan;

import std;

import atom;
import component;
import connectivity_table;
import fragment;
import fragment_graph;
import intra_molecular_potentials;

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

// Determine the next growth step given the set of already-placed beads.
//
// Mirrors 'ConnectivityTable::nextBeads' for the fully flexible case, but rigid-body fragments are
// hinged as one rigid body once their connecting atom is placed, and cyclic clusters (all rings:
// simple, fused, bridged) are grown as one ring-closure step.
static CBMC::GrowStep nextGrowthStep(const Component &component, const std::vector<std::size_t> &placedBeads)
{
  const ConnectivityTable &connectivity = component.connectivityTable;
  const FragmentGraph &graph = component.fragmentGraph;
  std::size_t numberOfBeads = connectivity.numberOfBeads;

  std::vector<bool> placed(numberOfBeads, false);
  for (std::size_t bead : placedBeads) placed[bead] = true;

  // 'filteredInteractions' takes mutable spans, so keep a mutable copy of the placed beads.
  std::vector<std::size_t> placedVec(placedBeads.begin(), placedBeads.end());

  auto makeStep = [&](CBMC::GrowStep::Kind kind, std::optional<std::size_t> previousBead, std::size_t currentBead,
                      std::vector<std::size_t> nextBeads, bool rigidBody) -> CBMC::GrowStep
  {
    Potentials::IntraMolecularPotentials intra =
        component.intraMolecularPotentials.filteredInteractions(numberOfBeads, placedVec, nextBeads);
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
  std::make_signed_t<std::size_t> cluster = clusterOf(nextBead);
  if (cluster >= 0)
  {
    if (clusterOf(currentBead) == cluster)
    {
      // The anchor is a placed atom of the cluster: grow the remaining cluster atoms as one
      // ring-closure step (the closure bonds in 'intra' keep every ring of the cluster closed).
      const std::vector<std::size_t> &clusterAtoms = graph.cyclicClusters[static_cast<std::size_t>(cluster)];
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
    std::optional<std::size_t> previousBead =
        firstPlacedNeighbor(connectivity, currentBead, placed, [](std::size_t) { return true; });
    CBMC::GrowStep::Kind kind =
        previousBead.has_value() ? CBMC::GrowStep::Kind::AttachFragment : CBMC::GrowStep::Kind::PlaceSeedFragment;
    return makeStep(kind, previousBead, currentBead, {nextBead}, false);
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
      // flexible branch.
      if (!graph.fragments[graph.atomFragmentIds[i]].isRigidBody() && clusterOf(i) < 0)
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
    return makeStep(CBMC::GrowStep::Kind::PlaceSeedFragment, std::nullopt, currentBead, {nextBead}, false);
  }

  // More than one placed neighbor is only allowed when 'currentBead' is part of an already-placed
  // rigid fragment or cyclic cluster (e.g. a flexible tail growing off a ring atom): pick the first
  // placed neighbor as the reference. For a purely flexible acyclic anchor this signals an
  // inconsistent fragment graph.
  if (numberOfPreviousBeads > 1 && !graph.fragments[graph.atomFragmentIds[currentBead]].isRigidBody() &&
      clusterOf(currentBead) < 0)
  {
    throw std::runtime_error(std::format("Error in CBMC: Multiple previous beads\n"));
  }

  return makeStep(CBMC::GrowStep::Kind::AttachFragment, previousBead, currentBead, std::move(nextBeads), false);
}

std::vector<CBMC::GrowStep> CBMC::buildGrowthPlan(const Component &component,
                                                  const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  std::size_t numberOfBeads = component.connectivityTable.numberOfBeads;

  std::vector<GrowStep> plan{};
  std::vector<std::size_t> placed(beadsAlreadyPlaced.begin(), beadsAlreadyPlaced.end());

  while (placed.size() < numberOfBeads)
  {
    GrowStep step = nextGrowthStep(component, placed);
    placed.insert(placed.end(), step.nextBeads.begin(), step.nextBeads.end());
    plan.push_back(std::move(step));
  }

  return plan;
}
