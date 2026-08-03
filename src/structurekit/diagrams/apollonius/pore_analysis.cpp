module;

module apollonius_pore_analysis;

import std;

import double3;
import crystal;
import pair_interactions;
import apollonius_network;
import voronoi_pore_diameters;
import voronoi_channels;
import pore_window;

void ApolloniusPoreAnalysis::run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("ApolloniusPoreAnalysis: Unknown probe-atom type\n");
  }
  double probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  fractionalPositions.reserve(framework.atoms.size());
  radii.reserve(framework.atoms.size());
  for (const CrystalAtom& atom : framework.atoms)
  {
    fractionalPositions.push_back(framework.unitCell.inverseCell * atom.position);
    std::size_t type = atom.type;
    radii.push_back(0.5 * interactions(type, type).sizeParameter);
  }

  ApolloniusPoreNetwork pores = ApolloniusPoreNetwork::create(framework.unitCell, fractionalPositions, radii);

  diameters = PoreDiameters::compute(pores.network);
  channels = ChannelAnalysis::compute(pores.network, probeRadius);
  freeSphere = freeSphereWindow(pores.network);
  windows = channelWindows(pores.network, channels);

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  std::ofstream diameterFile;
  diameterFile.open(framework.name + ".apollonius.res.txt");
  std::print(diameterFile, "# Pore diameters from the Apollonius diagram (Di, Df, Dif)\n");
  std::print(diameterFile, "# Crystal: {}\n", framework.name);
  std::print(diameterFile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(diameterFile, "# Number of framework atoms: {}\n", framework.atoms.size());
  pores.writeHeader(diameterFile);
  std::print(diameterFile, "# CPU Timing: {} [s]\n", timing.count());
  std::print(diameterFile, "Di (largest included sphere):            {} [Å]\n", diameters.includedSphereDiameter);
  std::print(diameterFile, "Df (largest free sphere):                {} [Å]\n", diameters.freeSphereDiameter);
  std::print(diameterFile, "Dif (included sphere along free path):   {} [Å]\n",
             diameters.includedAlongFreePathDiameter);

  std::print(diameterFile, "\n");
  std::print(diameterFile, "# The shape of the narrowest window, measured in the plane across the passage where the\n");
  std::print(diameterFile, "# bottleneck sits. Df is the inscribed circle of that window and says nothing about how far\n");
  std::print(diameterFile, "# it reaches the other ways, which these two pairs of numbers do: the free width is the\n");
  std::print(diameterFile, "# shortest and the longest free chord through the bottleneck, and the ellipse is the\n");
  std::print(diameterFile, "# largest-area one that fits there, which may be narrower than Df one way in return for\n");
  std::print(diameterFile, "# being longer the other. Both are lower bounds, each direction being followed only as far\n");
  std::print(diameterFile, "# as the first atom it meets.\n");
  std::print(diameterFile, "#\n");
  std::print(diameterFile, "# Neither is a criterion for a non-spherical molecule, whose passage turns on its own shape\n");
  std::print(diameterFile, "# and on how it may turn along the way. They describe the window.\n");

  auto writeWindow = [&](const PoreWindow& measured)
  {
    std::print(diameterFile, " free width {:.5f} - {:.5f} [Å] ellipse {:.5f} x {:.5f} [Å] ring {} atoms{}\n",
               measured.smallestFreeWidth, measured.largestFreeWidth, measured.minorAxis, measured.majorAxis,
               measured.boundingAtoms,
               measured.clipped ? " (open in some direction, cut off at the reach of the cell)" : "");
    std::print(diameterFile,
               "  window at ({:.5f}, {:.5f}, {:.5f}) [Å], across ({:.5f}, {:.5f}, {:.5f}), major axis along "
               "({:.5f}, {:.5f}, {:.5f})\n",
               measured.position.x, measured.position.y, measured.position.z, measured.normal.x, measured.normal.y,
               measured.normal.z, measured.majorAxisDirection.x, measured.majorAxisDirection.y,
               measured.majorAxisDirection.z);
  };

  std::print(diameterFile, "Df window:");
  if (freeSphere.measured)
    writeWindow(freeSphere);
  else
    std::print(diameterFile, " none, the structure does not percolate\n");

  // The channels are those a probe of this size can travel, so a channel's own bottleneck is its own and
  // is at most the one above.
  std::print(diameterFile, "\n# Per channel, at the bottleneck of that channel's own widest percolating path for probe\n");
  std::print(diameterFile, "# {} of radius {} [Å]. Channels are numbered as in the .chan file.\n", probePseudoAtom,
             probeRadius);
  if (windows.empty())
  {
    std::print(diameterFile, "no channels for this probe\n");
  }
  for (const ChannelWindow& channel : windows)
  {
    std::print(diameterFile, "channel {} (dimensionality {}): Df {:.5f} [Å]", channel.poreIndex,
               channel.dimensionality, channel.freeSphereDiameter);
    if (!channel.window.measured)
    {
      std::print(diameterFile, " window not measured (bottleneck clearance {:.5f} [Å])\n", channel.window.freeRadius);
      continue;
    }
    writeWindow(channel.window);
  }
  diameterFile.close();

  std::ofstream channelFile;
  channelFile.open(framework.name + ".apollonius.chan.txt");
  std::print(channelFile, "# Channel and pocket analysis from the Apollonius diagram\n");
  std::print(channelFile, "# Crystal: {}\n", framework.name);
  std::print(channelFile, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, probeRadius);
  pores.writeHeader(channelFile);
  std::print(channelFile, "Number of channels: {}\n", channels.numberOfChannels);
  std::print(channelFile, "Number of pockets:  {}\n", channels.numberOfPockets);
  for (std::size_t i = 0; i < channels.pores.size(); ++i)
  {
    const VoronoiPore& pore = channels.pores[i];
    std::print(channelFile, "  pore {}: {} dimensionality={} nodes={}\n", i, pore.isChannel ? "channel" : "pocket",
               pore.dimensionality, pore.nodeIndices.size());
  }
  channelFile.close();
}
