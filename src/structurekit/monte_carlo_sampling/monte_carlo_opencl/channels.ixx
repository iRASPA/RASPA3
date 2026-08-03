module;

export module mc_opencl_channels;

import std;

import sampled_structure;
import sampled_roadmap;
import mc_channels;

// The channels, the pockets and the dimensionality from a roadmap of the void built on the GPU.
//
// The reading itself is `MC_Channels`, which this is, and it is not repeated here. What a roadmap costs is
// the arithmetic that goes into building it -- a clearance at every point thrown, a measured way along
// every hop between them, a walk uphill from every peak and a search for a bent way between every pair of
// pockets -- and that is what moves to the device. What is left over is a graph to walk about on, which is
// branchy, serial and cheap, and which is the same code either way.
//
// So this is not a second implementation to be kept in step with the first. It is the first one, reading a
// roadmap whose spheres were measured somewhere else, and the fields it leaves behind are the same fields
// in the same units. The files it writes are named `.gpu.` where the processor's are named `.cpu.`, so that
// a run of each leaves both to be compared rather than one on top of the other.
export struct MC_OpenCL_Channels : MC_Channels
{
  // The reading of a roadmap that was built already is inherited unchanged: a roadmap does not care what
  // measured it, and one built on the device can be read here or handed to the processor's own reader.
  using MC_Channels::run;

  void run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
           std::optional<std::size_t> numberOfInnerSteps);
};
