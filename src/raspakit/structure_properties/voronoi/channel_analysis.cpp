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
  // Gaussian elimination on the vectors as rows; count non-zero pivot rows.
  std::array<std::array<double, 3>, 3> basis{};
  int rank = 0;
  for (const int3& v : vectors)
  {
    std::array<double, 3> row{static_cast<double>(v.x), static_cast<double>(v.y), static_cast<double>(v.z)};
    for (int r = 0; r < rank; ++r)
    {
      double dot = row[0] * basis[r][0] + row[1] * basis[r][1] + row[2] * basis[r][2];
      double norm = basis[r][0] * basis[r][0] + basis[r][1] * basis[r][1] + basis[r][2] * basis[r][2];
      if (norm > 0.0)
      {
        double factor = dot / norm;
        for (int d = 0; d < 3; ++d) row[d] -= factor * basis[r][d];
      }
    }
    double residual = std::sqrt(row[0] * row[0] + row[1] * row[1] + row[2] * row[2]);
    if (residual > 1.0e-6)
    {
      basis[rank] = row;
      ++rank;
      if (rank == 3) break;
    }
  }
  return rank;
}

ChannelAnalysis ChannelAnalysis::compute(const VoronoiNetwork& network, double probeRadius)
{
  ChannelAnalysis analysis;
  const std::size_t numberOfNodes = network.nodes.size();
  analysis.nodePoreId.assign(numberOfNodes, -1);

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
    for (std::size_t node : componentNodes) analysis.nodePoreId[node] = poreId;
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
  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("VoronoiChannels: Unknown probe-atom type\n");
  }
  double probeRadius = 0.5 * forceField[probeType.value()].sizeParameter();

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (const Atom& atom : framework.unitCellAtoms)
  {
    fractionalPositions.push_back(framework.simulationBox.inverseCell * atom.position);
    std::size_t type = static_cast<std::size_t>(atom.type);
    radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }

  VoronoiNetwork network = VoronoiNetwork::create(framework.simulationBox, fractionalPositions, radii);
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
