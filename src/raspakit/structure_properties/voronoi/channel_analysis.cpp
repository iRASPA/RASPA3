module;

module voronoi_channels;

import std;

import int3;
import double3;
import atom;
import framework;
import forcefield;
import voronoi_network;

int latticeVectorRank(const std::vector<int3>& vectors)
{
  // In integers throughout. These are lattice translations, so how many directions they span is a property of
  // the integers themselves and there is no tolerance to choose: a vector is off a line when its cross product
  // with that line does not vanish, and off a plane when its product with the plane's normal does not, and
  // both tests are decided rather than estimated. Rounding here would be a threshold on how nearly parallel
  // two directions of a channel system may be before they are called one, and there is no such thing.
  using vector = std::array<std::int64_t, 3>;
  auto cross = [](const vector& u, const vector& v) {
    return vector{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
  };
  auto nonZero = [](const vector& u) { return u[0] != 0 || u[1] != 0 || u[2] != 0; };

  vector along{}, normal{};
  int rank = 0;
  for (const int3& v : vectors)
  {
    vector row{v.x, v.y, v.z};
    if (!nonZero(row)) continue;

    if (rank == 0)
    {
      along = row;
      rank = 1;
    }
    else if (rank == 1)
    {
      vector off = cross(along, row);
      if (nonZero(off))
      {
        normal = off;
        rank = 2;
      }
    }
    else if (normal[0] * row[0] + normal[1] * row[1] + normal[2] * row[2] != 0)
    {
      return 3;
    }
  }
  return rank;
}

ChannelAnalysis ChannelAnalysis::compute(const VoronoiNetwork& network, double probeRadius)
{
  ChannelAnalysis analysis;
  const std::size_t numberOfNodes = network.nodes.size();
  analysis.nodePoreId.assign(numberOfNodes, -1);
  analysis.nodeLatticeOffset.assign(numberOfNodes, int3(0, 0, 0));

  // Adjacency list restricted to nodes/edges wide enough for the probe.
  std::vector<bool> nodeActive(numberOfNodes, false);
  for (std::size_t i = 0; i < numberOfNodes; ++i) nodeActive[i] = network.nodes[i].radius > probeRadius;

  std::vector<std::vector<std::pair<std::size_t, int3>>> adjacency(numberOfNodes);
  for (const VoronoiEdge& edge : network.edges)
  {
    if (edge.radius <= probeRadius) continue;
    if (!nodeActive[edge.from] || !nodeActive[edge.to]) continue;
    adjacency[edge.from].push_back({edge.to, edge.delta});
  }

  std::vector<bool> visited(numberOfNodes, false);

  for (std::size_t start = 0; start < numberOfNodes; ++start)
  {
    if (!nodeActive[start] || visited[start]) continue;

    // Depth-first walk tracking the accumulated periodic displacement of each node.
    std::unordered_map<std::size_t, int3> displacement;
    std::vector<std::size_t> componentNodes;
    std::vector<int3> loopVectors;

    std::vector<std::pair<std::size_t, int3>> stack;
    stack.push_back({start, int3(0, 0, 0)});
    displacement[start] = int3(0, 0, 0);
    visited[start] = true;

    while (!stack.empty())
    {
      auto [node, disp] = stack.back();
      stack.pop_back();
      componentNodes.push_back(node);

      for (const auto& [neighbor, delta] : adjacency[node])
      {
        int3 newDisp = disp + delta;
        auto it = displacement.find(neighbor);
        if (it == displacement.end())
        {
          displacement[neighbor] = newDisp;
          visited[neighbor] = true;
          stack.push_back({neighbor, newDisp});
        }
        else
        {
          int3 loop = newDisp - it->second;
          if (loop.x != 0 || loop.y != 0 || loop.z != 0) loopVectors.push_back(loop);
        }
      }
    }

    int dimensionality = latticeVectorRank(loopVectors);

    VoronoiPore pore;
    pore.isChannel = dimensionality > 0;
    pore.dimensionality = dimensionality;
    pore.nodeIndices = componentNodes;

    std::int32_t poreId = static_cast<std::int32_t>(analysis.pores.size());
    for (std::size_t node : componentNodes)
    {
      analysis.nodePoreId[node] = poreId;
      analysis.nodeLatticeOffset[node] = displacement[node];
    }
    if (pore.isChannel)
      ++analysis.numberOfChannels;
    else
      ++analysis.numberOfPockets;
    analysis.pores.push_back(std::move(pore));
  }

  return analysis;
}

void VoronoiChannels::run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom)
{
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (const Atom& atom : framework.unitCellAtoms)
  {
    fractionalPositions.push_back(framework.simulationBox.inverseCell * atom.position);
    std::size_t type = static_cast<std::size_t>(atom.type);
    radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }

  run(forceField, framework, probePseudoAtom,
      VoronoiNetwork::create(framework.simulationBox, fractionalPositions, radii));
}

void VoronoiChannels::run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
                          const VoronoiNetwork& network)
{
  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("VoronoiChannels: Unknown probe-atom type\n");
  }
  double probeRadius = 0.5 * forceField[probeType.value()].sizeParameter();

  result = ChannelAnalysis::compute(network, probeRadius);

  std::ofstream myfile;
  myfile.open(framework.name + ".voronoi.chan.txt");
  std::print(myfile, "# Channel and pocket analysis from the Voronoi network\n");
  std::print(myfile, "# Framework: {}\n", framework.name);
  std::print(myfile, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, probeRadius);
  std::print(myfile, "# Number of Voronoi nodes: {}\n", network.nodes.size());
  std::print(myfile, "Number of channels: {}\n", result.numberOfChannels);
  std::print(myfile, "Number of pockets:  {}\n", result.numberOfPockets);
  for (std::size_t i = 0; i < result.pores.size(); ++i)
  {
    const VoronoiPore& pore = result.pores[i];
    std::print(myfile, "  pore {}: {} dimensionality={} nodes={}\n", i,
               pore.isChannel ? "channel" : "pocket", pore.dimensionality, pore.nodeIndices.size());
  }
  myfile.close();
}
