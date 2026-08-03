module;

export module mc_channels;

import std;

import channel_analysis;
import sampled_structure;
import sampled_roadmap;

// Which way the pores of a crystal run, and how many of them go nowhere, by sampling.
//
// The dimensionality of a pore system is an integer and it is not a measurement: a channel either closes on
// a translate of itself along some lattice direction or it does not, and how many independent directions it
// does that in is 1, 2 or 3. So there is nothing here for a sample to converge to. What the sample decides
// is which points belong to the same pore, and once that is settled the rank is exact -- `ChannelAnalysis`
// is the same walk over the same integers whether the graph beneath it came from a Voronoi tessellation, an
// Apollonius diagram or a roadmap of thrown points.
//
// Because it is a rank, it can only be got wrong in one direction. A hop the sample failed to find breaks a
// pore into pieces, and pieces run in no more directions than the whole did and often fewer; a hop it did
// find is a hop a sphere was shown to make. So a sampled dimensionality is a lower bound on the true one
// and a sampled pocket count an upper bound, and a run that finds a 3D system has found three directions a
// probe can really travel in.
//
// The pruning the diagram routes do for a probe radius has already happened here: the contact radii the
// structure carries have the probe in them, so a point kept at all is a point the probe fits at.
export struct MC_Channels
{
  ChannelAnalysis result;

  int dimensionality{0};            // of the pore system as a whole: the widest-running of its channels
  std::size_t numberOfChannels{0};  // components that run off across the crystal
  std::size_t numberOfPockets{0};   // those that close on themselves and were resolved
  std::size_t numberOfUnresolvedPieces{0};  // and those too thinly sampled to call either way

  // Of the points that landed in the void, how many are in a channel. A pore is one component however many
  // points fell in it, so counting components says nothing about how much of the void runs anywhere.
  double channelShareOfVoid{0.0};

  SampledRoadmap roadmap;
  double seconds{0.0};

  void run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
           std::optional<std::size_t> numberOfInnerSteps);

  void run(const SampledStructure &structure, const SampledRoadmap &roadmap);
};
