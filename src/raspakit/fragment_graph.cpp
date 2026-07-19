module;

module fragment_graph;

import std;

import archive;
import double3;
import atom;
import molecule;
import connectivity_table;
import fragment;

std::vector<std::vector<std::size_t>> FragmentGraph::partitionAtoms(
    std::size_t numberOfBeads, const std::vector<std::vector<std::size_t>> &rigidBodies)
{
  std::vector<std::make_signed_t<std::size_t>> assigned(numberOfBeads, -1);

  // Rigid bodies get a stable temporary id; every listed atom must be in range and listed once.
  std::vector<std::vector<std::size_t>> rigidBodyAtoms(rigidBodies.size());
  for (std::size_t b = 0; b != rigidBodies.size(); ++b)
  {
    for (std::size_t atom : rigidBodies[b])
    {
      if (atom >= numberOfBeads)
      {
        throw std::runtime_error(std::format(
            "[FragmentGraph]: rigid-body atom index {} out of range (molecule has {} atoms)\n", atom, numberOfBeads));
      }
      if (assigned[atom] != -1)
      {
        throw std::runtime_error(
            std::format("[FragmentGraph]: atom {} is listed in more than one rigid body\n", atom));
      }
      assigned[atom] = static_cast<std::make_signed_t<std::size_t>>(b);
    }
    rigidBodyAtoms[b] = rigidBodies[b];
  }

  // Emit fragments ordered by their lowest atom index: when scanning atoms in ascending order, the
  // first time we reach an atom of a rigid body we emit that whole body; an unlisted atom becomes a
  // singleton fragment.
  std::vector<std::vector<std::size_t>> partition{};
  std::vector<bool> bodyEmitted(rigidBodies.size(), false);
  for (std::size_t atom = 0; atom != numberOfBeads; ++atom)
  {
    if (assigned[atom] == -1)
    {
      partition.push_back({atom});
      continue;
    }
    std::size_t body = static_cast<std::size_t>(assigned[atom]);
    if (bodyEmitted[body]) continue;
    bodyEmitted[body] = true;

    std::vector<std::size_t> sortedBody = rigidBodyAtoms[body];
    std::sort(sortedBody.begin(), sortedBody.end());
    partition.push_back(std::move(sortedBody));
  }

  return partition;
}

void FragmentGraph::build(const ConnectivityTable &connectivity,
                          const std::vector<std::vector<std::size_t>> &partition, std::size_t startingBead,
                          std::span<const double3> referencePositions, std::span<const double> masses)
{
  std::size_t numberOfBeads = connectivity.numberOfBeads;

  fragments.clear();
  fragments.reserve(partition.size());
  atomFragmentIds.assign(numberOfBeads, 0);

  for (std::size_t f = 0; f != partition.size(); ++f)
  {
    Fragment fragment(partition[f]);
    fragment.computeRigidProperties(referencePositions, masses);
    for (std::size_t atom : partition[f]) atomFragmentIds[atom] = f;
    fragments.push_back(std::move(fragment));
  }

  std::size_t numberOfFragments = fragments.size();
  rootFragment = (startingBead < numberOfBeads) ? atomFragmentIds[startingBead] : 0;

  // Inter-fragment bonds, each stored once with atomA < atomB, in ascending order.
  std::vector<std::array<std::size_t, 2>> interFragmentBonds{};
  for (std::size_t i = 0; i != numberOfBeads; ++i)
  {
    for (std::size_t j = i + 1; j != numberOfBeads; ++j)
    {
      if (connectivity[i, j] && atomFragmentIds[i] != atomFragmentIds[j])
      {
        interFragmentBonds.push_back({i, j});
      }
    }
  }

  // Adjacency: per fragment, the incident inter-fragment bonds with the local/remote atom identified.
  struct Incidence
  {
    std::size_t bondIndex;
    std::size_t neighborFragment;
    std::size_t localAtom;
    std::size_t remoteAtom;
  };
  std::vector<std::vector<Incidence>> adjacency(numberOfFragments);
  for (std::size_t b = 0; b != interFragmentBonds.size(); ++b)
  {
    std::size_t a = interFragmentBonds[b][0];
    std::size_t c = interFragmentBonds[b][1];
    std::size_t fa = atomFragmentIds[a];
    std::size_t fc = atomFragmentIds[c];
    adjacency[fa].push_back({b, fc, a, c});
    adjacency[fc].push_back({b, fa, c, a});
  }

  parentFragment.assign(numberOfFragments, std::nullopt);
  parentBond.assign(numberOfFragments, std::nullopt);
  growthOrder.clear();
  growthOrder.reserve(numberOfFragments);
  closureBonds.clear();

  std::vector<bool> visited(numberOfFragments, false);
  std::vector<bool> bondUsed(interFragmentBonds.size(), false);

  // Deterministic BFS starting from the root fragment, then any remaining components (ascending
  // fragment index) so a disconnected molecule still yields a full forest.
  auto bfs = [&](std::size_t start)
  {
    std::queue<std::size_t> queue;
    visited[start] = true;
    queue.push(start);
    growthOrder.push_back(start);
    while (!queue.empty())
    {
      std::size_t f = queue.front();
      queue.pop();
      for (const Incidence &inc : adjacency[f])
      {
        if (!visited[inc.neighborFragment])
        {
          visited[inc.neighborFragment] = true;
          bondUsed[inc.bondIndex] = true;
          parentFragment[inc.neighborFragment] = f;
          parentBond[inc.neighborFragment] = std::array<std::size_t, 2>{inc.localAtom, inc.remoteAtom};
          growthOrder.push_back(inc.neighborFragment);
          queue.push(inc.neighborFragment);
        }
      }
    }
  };

  bfs(rootFragment);
  for (std::size_t f = 0; f != numberOfFragments; ++f)
  {
    if (!visited[f]) bfs(f);
  }

  for (std::size_t b = 0; b != interFragmentBonds.size(); ++b)
  {
    if (!bondUsed[b]) closureBonds.push_back(interFragmentBonds[b]);
  }

  // Cyclic clusters: the 2-edge-connected components of the fragment graph. Every closure bond,
  // together with the spanning-tree path between its endpoint fragments, forms a cycle; overlapping
  // cycles merge into one cluster (fused/bridged polycyclic systems).
  std::vector<std::size_t> depth(numberOfFragments, 0);
  for (std::size_t f : growthOrder)
  {
    depth[f] = parentFragment[f].has_value() ? depth[parentFragment[f].value()] + 1 : 0;
  }

  std::vector<std::size_t> unionParent(numberOfFragments);
  std::iota(unionParent.begin(), unionParent.end(), std::size_t{0});
  std::function<std::size_t(std::size_t)> findRoot = [&](std::size_t x) -> std::size_t
  {
    while (unionParent[x] != x)
    {
      unionParent[x] = unionParent[unionParent[x]];
      x = unionParent[x];
    }
    return x;
  };
  auto unite = [&](std::size_t a, std::size_t b) { unionParent[findRoot(a)] = findRoot(b); };

  std::vector<bool> onCycle(numberOfFragments, false);
  for (const std::array<std::size_t, 2> &closure : closureBonds)
  {
    std::size_t fa = atomFragmentIds[closure[0]];
    std::size_t fb = atomFragmentIds[closure[1]];
    onCycle[fa] = true;
    onCycle[fb] = true;
    while (fa != fb)
    {
      if (depth[fa] < depth[fb]) std::swap(fa, fb);
      std::size_t up = parentFragment[fa].value();
      onCycle[up] = onCycle[fa] = true;
      unite(fa, up);
      fa = up;
    }
  }

  cyclicClusters.clear();
  fragmentCyclicClusterIds.assign(numberOfFragments, -1);
  std::map<std::size_t, std::size_t> rootToCluster{};
  for (std::size_t f = 0; f != numberOfFragments; ++f)
  {
    if (!onCycle[f]) continue;
    std::size_t root = findRoot(f);
    auto [it, inserted] = rootToCluster.insert({root, cyclicClusters.size()});
    if (inserted) cyclicClusters.push_back({});
    fragmentCyclicClusterIds[f] = static_cast<std::make_signed_t<std::size_t>>(it->second);
    for (std::size_t atom : fragments[f].atoms) cyclicClusters[it->second].push_back(atom);
  }
  for (std::vector<std::size_t> &cluster : cyclicClusters) std::sort(cluster.begin(), cluster.end());

  finalizeCache();
}

void FragmentGraph::finalizeCache()
{
  rigidFragmentCount = 0;
  for (const Fragment &fragment : fragments)
  {
    if (fragment.isRigidBody()) ++rigidFragmentCount;
  }
  semiFlexible = (rigidFragmentCount > 0) && (fragments.size() > 1);
}

bool FragmentGraph::isInsideRigidFragment(std::span<const std::size_t> ids) const
{
  if (ids.empty() || atomFragmentIds.empty()) return false;
  std::size_t f = atomFragmentIds[ids[0]];
  if (!fragments[f].isRigidBody()) return false;
  for (std::size_t id : ids)
  {
    if (atomFragmentIds[id] != f) return false;
  }
  return true;
}

std::optional<std::vector<std::size_t>> FragmentGraph::nextGrowBeads(
    const std::vector<std::size_t> &placedBeads, const ConnectivityTable &connectivity) const
{
  std::size_t numberOfBeads = connectivity.numberOfBeads;

  std::vector<bool> placed(numberOfBeads, false);
  for (std::size_t bead : placedBeads) placed[bead] = true;

  // First frontier bond: a placed bead 'k' bonded to an unplaced bead 'j' (placed beads in order,
  // neighbors by ascending index).
  std::optional<std::size_t> nextBead{};
  for (std::size_t k : placedBeads)
  {
    for (std::size_t j = 0; j != numberOfBeads; ++j)
    {
      if (connectivity[k, j] && !placed[j])
      {
        nextBead = j;
        break;
      }
    }
    if (nextBead.has_value()) break;
  }
  if (!nextBead.has_value()) return std::nullopt;

  // Grow the whole fragment that contains the frontier bead. Every other atom of that fragment must
  // still be unplaced: a fragment is placed as a single unit.
  std::size_t fragmentIndex = atomFragmentIds[nextBead.value()];
  const Fragment &fragment = fragments[fragmentIndex];

  std::vector<std::size_t> nextBeads{};
  nextBeads.reserve(fragment.atoms.size());
  for (std::size_t atom : fragment.atoms)
  {
    if (placed[atom]) return std::nullopt;
    nextBeads.push_back(atom);
  }
  return nextBeads;
}

std::vector<std::size_t> FragmentGraph::rigidFragmentIndices() const
{
  std::vector<std::size_t> result{};
  for (std::size_t f = 0; f != fragments.size(); ++f)
  {
    if (fragments[f].isRigidBody()) result.push_back(f);
  }
  return result;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const FragmentGraph &g)
{
  archive << g.versionNumber;
  archive << g.fragments;
  archive << g.atomFragmentIds;
  archive << g.rootFragment;
  archive << g.parentFragment;
  archive << g.parentBond;
  archive << g.growthOrder;
  archive << g.closureBonds;
  archive << g.cyclicClusters;
  archive << g.fragmentCyclicClusterIds;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, FragmentGraph &g)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > g.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'FragmentGraph' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> g.fragments;
  archive >> g.atomFragmentIds;
  archive >> g.rootFragment;
  archive >> g.parentFragment;
  archive >> g.parentBond;
  archive >> g.growthOrder;
  archive >> g.closureBonds;
  archive >> g.cyclicClusters;
  archive >> g.fragmentCyclicClusterIds;
  g.finalizeCache();
  return archive;
}
