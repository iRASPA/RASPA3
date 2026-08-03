module;

module mc_channels;

import std;

import channel_analysis;
import sampled_structure;
import sampled_roadmap;
import mc_backend;

void MC_Channels::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                      std::optional<std::size_t> numberOfInnerSteps)
{
  run(structure, SampledRoadmap::build(structure, samplingBackendCPU(), numberOfIterations, numberOfInnerSteps));
}

void MC_Channels::run(const SampledStructure &structure, const SampledRoadmap &roadmap)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  this->roadmap = roadmap;
  this->result = roadmap.components;

  this->dimensionality = this->result.dimensionality();
  this->numberOfChannels = roadmap.numberOfChannels;
  this->numberOfPockets = roadmap.numberOfPockets;
  this->numberOfUnresolvedPieces = roadmap.numberOfUnresolvedPieces;

  std::size_t inAChannel = 0;
  std::size_t volumeNodes = roadmap.numberOfVolumeNodes();
  for (std::size_t node = 0; node < volumeNodes; ++node)
  {
    if (roadmap.isReachable(node)) ++inAChannel;
  }
  this->channelShareOfVoid =
      volumeNodes > 0 ? static_cast<double>(inAChannel) / static_cast<double>(volumeNodes) : 0.0;

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;
  this->seconds = timing.count();

  std::ofstream myfile;
  myfile.open(std::format("{}.mc.chan.{}.txt", structure.name, roadmap.backend));
  std::print(myfile, "# Channel and pocket analysis by sampling the void\n");
  structure.writeHeader(myfile);
  roadmap.writeHeader(myfile);
  std::print(myfile, "# Timing (channels, on the processor either way): {} [s]\n", this->seconds);
  std::print(myfile, "# A dimensionality here is a lower bound and a pocket count an upper one: hops the\n");
  std::print(myfile, "# sample missed break a pore into pieces, and no piece runs in more directions than\n");
  std::print(myfile, "# the whole did. Every direction reported is one a sphere was shown to travel in.\n");
  std::print(myfile, "Number of channels: {}\n", this->numberOfChannels);
  std::print(myfile, "Number of pockets:  {}\n", this->numberOfPockets);
  std::print(myfile, "Number of pieces too thinly sampled to say: {}\n", this->numberOfUnresolvedPieces);
  std::print(myfile, "Pore system dimensionality: {}\n", this->dimensionality);
  std::print(myfile, "Share of the void in a channel: {}\n", this->channelShareOfVoid);
  for (std::size_t i = 0; i < this->result.pores.size(); ++i)
  {
    const VoronoiPore &pore = this->result.pores[i];
    std::print(myfile, "  pore {}: {} dimensionality={} nodes={}\n", i,
               pore.isChannel               ? "channel"
               : this->roadmap.poreIsResolved[i] != 0 ? "pocket"
                                                      : "unresolved",
               pore.dimensionality, pore.nodeIndices.size());
  }
  myfile.close();
}
