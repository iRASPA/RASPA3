module;

export module mc_blocking_pockets;

import std;

import blocking_spheres;
import sampled_structure;
import sampled_roadmap;

// Spheres covering the pockets a probe cannot reach, found by sampling, for a simulation to reject
// insertions in.
//
// The two things a blocking sphere has to do pull against each other: between them the spheres of a pocket
// must hold all of it, or a molecule still gets into the part they miss, and none of them may reach a
// channel, or the simulation loses part of its own pore. Sampling answers both from the points themselves.
// Every point here is a position for the probe's *centre*, the probe already being in the contact radii, so
// distances between points are in those terms and a radius read off them needs no further allowance for the
// probe's size.
//
// A sphere reaches as far as the furthest point of its pocket still uncovered and no further than the
// nearest point a probe can occupy from a channel; and below both of those sits a floor that holds
// regardless, the clearance at its own centre, since a sphere of that radius holds no atom and so lies in
// the free space, and being connected it cannot cross the neck that made the pocket a pocket. Both of the
// first two carry an allowance for how finely the void was sampled, and both allowances fall away as the
// sample grows, which is what keeps the radii a property of the structure rather than of the sample size.
//
// Where this stands relative to the exact routes is worth being plain about, because the direction of the
// error is not the same as it is for a diameter. A pocket no point landed in is a pocket left open: the
// diagram routes cannot miss one, having a vertex in every cavity however small, and this one can. Against
// that, a sphere that is written was placed among points measured to be shut in, and the roadmap has
// already been made to try every way out of every piece it found before calling it shut, so a sphere over
// part of a working pore takes both a missed way and a sample too thin to notice. The count of pockets and
// the fraction of the void they hold are reported for that reason: they are what says whether the sample
// was equal to the structure.
export struct MC_BlockingPockets
{
  std::vector<BlockingSphere> spheres;

  std::size_t numberOfPockets{0};           // pieces of the void that run nowhere
  std::size_t numberOfCoveredPockets{0};    // of those, the ones a sphere was written for
  std::size_t numberOfUnresolvedPieces{0};  // pieces too thinly sampled to say, which are left open
  double blockedFractionOfVoid{0.0};        // the share of the sampled void shut in a pocket

  double resolution{0.0};  // Å, the mean spacing of the samples, which is what the sample can resolve

  SampledRoadmap roadmap;
  double seconds{0.0};

  void run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
           std::optional<std::size_t> numberOfInnerSteps);

  void run(const SampledStructure &structure, const SampledRoadmap &roadmap);
};
