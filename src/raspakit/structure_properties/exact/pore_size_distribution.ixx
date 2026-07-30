module;

export module exact_pore_size_distribution;

import std;

import double3;
import simulationbox;
import pore_accessibility;
import exact_solvent_excluded;

// The pore-size distribution of Gelb and Gubbins, evaluated rather than sampled.
//
// One probe radius at a time, the volume the probe opens up and the derivative of that volume are both closed
// forms over the boundary of the union of the inflated atoms; see exact_solvent_excluded for what they are and
// why. What is left here is the sweep over probe radii and the bookkeeping that turns the two into the curve
// that is usually plotted: the cumulative pore volume against pore diameter, and its derivative.
//
// Nothing is binned. Each row is the distribution at that diameter, evaluated there, so the spacing of the
// rows is the resolution of the plot and not the resolution of the answer. That is the practical difference
// from the sampled route: halving the spacing costs twice the work and buys twice the detail, where halving
// the bin width of a histogram costs nothing and buys noise.
export struct PoreSizeDistributionPoint
{
  double diameter{0.0};  // Å, twice the probe radius

  double poreVolume{0.0};      // Å³, the volume of the union of the balls of this radius that fit in the void
  double cumulative{0.0};      // that volume over the void volume, so one at zero diameter and zero beyond D_i
  double distribution{0.0};    // 1/Å, minus the derivative of the cumulative in the diameter

  // The share of the distribution on surfaces a probe of this row's own diameter can reach, on surfaces sealed
  // off from it, and on surfaces the network could not place. Every reentrant patch belongs to one connected
  // surface and that surface faces one pore, so this is a division of the boundary and not a classification of
  // samples. The probe here moves with the row: what these three divide is the volume leaving at diameter d
  // between the pores a probe of that same diameter d can get into and the ones it cannot.
  double accessible{0.0};
  double inaccessible{0.0};
  double undecided{0.0};

  // The same row for one fixed probe instead: the pore volume, at this diameter, of the part of the void that
  // probe can reach, and the derivative of it. The probe is the one the analysis was given, so these are the
  // distribution of the pore sizes a molecule of that size actually meets, where the three above are a
  // property of the geometry at every diameter separately.
  //
  // Below the probe's own diameter they are constant and zero respectively: the region the probe can occupy is
  // a union of balls of its radius, so every point of it has a pore size of at least that diameter and none of
  // it has left yet.
  double probeAccessiblePoreVolume{0.0};    // Å³
  double probeAccessibleCumulative{0.0};    // over the pore volume the probe can reach, so one below 2 r_p
  double probeAccessibleDistribution{0.0};  // 1/Å, normalised by the same volume

  // What the boundary was made of at this diameter, kept to be looked at rather than plotted.
  std::size_t numberOfArcs{0};
  std::size_t cuspedArcs{0};
  std::size_t numberOfVertices{0};
  std::size_t clippedVertices{0};
  std::size_t degenerateVertices{0};
};

// A pore size that holds volume of its own.
//
// The cumulative volume does not only slide down, it also falls off cliffs, and it is worth being clear that
// this is the definition speaking and not the arithmetic. Take a cavity that is a single ball of radius a. The
// probe of radius r fits in it while r <= a, and the union of the positions it then occupies is the whole
// cavity whatever r is, so the volume of pores at least 2r across is the whole cavity up to 2a and nothing
// above it. Every point of that cavity has the same pore size. The same happens in a cylinder, and it happens
// in a real pore too: a cage holds all of its volume at the diameter of the largest sphere that fits in it,
// and it is only the corrugation of the walls, the corners a maximal sphere cannot reach into, that puts any
// volume at smaller sizes.
//
// So the pore-size distribution of a crystal is not a function. It is a measure, with a continuous part from
// the corrugation and an atom at the largest sphere of every family of pores the void falls into. What a
// histogram of trial spheres shows there is a tall narrow peak whose height is the weight divided by the bin
// width, which is why the peak heights of a sampled distribution are not comparable between two runs that
// binned differently, and why they are the part of a published distribution that is hardest to read.
export struct PoreSizeSpike
{
  double diameter{0.0};  // Å
  double weight{0.0};    // the fraction of the void that has exactly this pore size
  double bracket{0.0};   // Å, the width of the interval it was finally cornered in
};

export struct PoreSizeDistributionCurve
{
  std::vector<PoreSizeDistributionPoint> points;

  double cellVolume{0.0};  // Å³
  double voidVolume{0.0};  // Å³, the pore volume at zero diameter, which normalises the distribution

  // The probe the accessible curve belongs to, and the volume that probe can reach: the pore volume at its own
  // diameter, less the pores it cannot get into. It is what normalises the accessible curve, so that curve
  // integrates to one over the region it describes rather than to the accessible share of the void.
  double probeRadius{0.0};            // Å
  double probeAccessibleVolume{0.0};  // Å³

  // The integral of the continuous part of the distribution, and the spikes, whose weights make up the rest.
  //
  // These come by two routes with no term in common: the first is the derivative, integrated over the
  // diameters; the second is what the volume lost between one diameter and the next that the derivative does
  // not account for. Together they have to come to one, since the void has some pore size everywhere, and how
  // near they come is the check of the volume against its own derivative over the whole range.
  double integral{0.0};
  double singularWeight{0.0};
  std::vector<PoreSizeSpike> spikes;

  // The same two, and the same spikes, for the part of the void the fixed probe can reach. They are found in
  // the same refinement of the same intervals, a volume that goes over at one size in a reachable pore doing so
  // in the total as well, and they are weights of the reachable volume rather than of the void.
  double probeAccessibleIntegral{0.0};
  double probeAccessibleSingularWeight{0.0};
  std::vector<PoreSizeSpike> probeAccessibleSpikes;

  // What is left at the end of the range, which ought to be nothing: a range that stops short of the largest
  // sphere in the void leaves volume unaccounted, and this is how much.
  double truncatedWeight{0.0};

  // The largest diameter at which any volume is left, which is the diameter of the largest sphere that fits
  // in the void. It is the position of the last spike, so it is resolved to that spike's bracket and not to
  // the spacing of the rows.
  double largestDiameter{0.0};

  // The same two for the accessible curve: what is left of the reachable volume at the end of the range, and the
  // largest sphere that fits anywhere the probe can reach, which is the largest pore of the framework a molecule
  // of that size ever meets.
  double probeAccessibleTruncatedWeight{0.0};
  double probeAccessibleLargestDiameter{0.0};

  double seconds{0.0};
  std::size_t numberOfEvaluations{0};
};

// The curve, over `numberOfBins` diameters evenly spaced up to `maximumDiameter`, taken at the midpoint of
// each step as the sampled routes report theirs.
//
// The spikes are found rather than assumed, and they are found by the one thing that tells an atom of the
// measure from a steep stretch of the continuous part: what happens under refinement. Wherever the volume lost
// across a step exceeds what the derivative accounts for, the step is bisected and the excess is looked for
// again in the halves. A corner of the continuous part leaves an excess that is the truncation error of the
// trapezium rule and falls off as the cube of the step, so it is gone within a few bisections. A cliff leaves
// the same excess however far the step is narrowed, and what survives `refinements` of them is a spike, its
// weight known to the trapezium error over an interval that narrow.
//
// Two curves come out of the one sweep. The first is the whole void: every pore of the framework, whether or
// not anything can get into it, which is the distribution the volume of the cell is made of. The second is the
// void one fixed probe can reach, `probeRadius` being that probe, which is the distribution a molecule of that
// size meets. They differ by the pores that are sealed off from that probe, and the difference is not a
// rescaling: a sealed cage can be the largest pore in a framework.
//
// The second is not the accessible share of the first row by row. At a diameter d the first divides the volume
// leaving there by what a probe of that same diameter can reach, which moves along the curve; the second holds
// the probe still. So a cage that the fixed probe cannot enter is outside the accessible curve at every
// diameter, and a bulge in a channel that only a larger sphere fits into is inside it, though at that diameter
// it is a pore of its own that nothing can get into.
//
// `build` makes the accessibility for a given probe radius: it is a callback because which diagram the pores
// are taken from is the caller's business and the geometry here is the same either way. It is called at
// `probeRadius` as well, the pores of that one network being what the accessible curve is divided by at every
// diameter above it, and once at vanishing probe for the void volume. The diameter of the largest sphere in
// that vanishing-probe network is where the cumulative hits zero, so rows of the report past it plus one bin
// are left at zero without being evaluated.
export PoreSizeDistributionCurve exactPoreSizeDistribution(
    const std::function<PoreAccessibility(double)>& build, double cellVolume, double maximumDiameter,
    std::size_t numberOfBins, std::size_t subdivisions, double probeRadius = 0.0, std::size_t refinements = 12);
