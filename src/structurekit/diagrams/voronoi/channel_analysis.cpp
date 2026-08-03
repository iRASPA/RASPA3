module;

module voronoi_channels;

import std;

import double3;
import crystal;
import pair_interactions;
import voronoi_network;
import channel_analysis;

void VoronoiChannels::run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom)
{
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (const CrystalAtom& atom : framework.atoms)
  {
    fractionalPositions.push_back(framework.unitCell.inverseCell * atom.position);
    std::size_t type = atom.type;
    radii.push_back(0.5 * interactions(type, type).sizeParameter);
  }

  run(interactions, framework, probePseudoAtom,
      VoronoiNetwork::create(framework.unitCell, fractionalPositions, radii));
}

void VoronoiChannels::run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom,
                          const VoronoiNetwork& network)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("VoronoiChannels: Unknown probe-atom type\n");
  }
  double probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  result = ChannelAnalysis::compute(network, probeRadius);

  std::ofstream myfile;
  myfile.open(framework.name + ".voronoi.chan.txt");
  std::print(myfile, "# Channel and pocket analysis from the Voronoi network\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
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
