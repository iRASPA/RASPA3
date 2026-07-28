module;

export module exact_pore_size_distribution;

import std;

import double3;
import simulationbox;
import voronoi_accessibility;
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

  // The share of the distribution carried on surfaces the probe can reach, on surfaces sealed off from it, and
  // on surfaces the network could not place. Every reentrant patch belongs to one connected surface and that
  // surface faces one pore, so this is a division of the boundary and not a classification of samples.
  double accessible{0.0};
  double inaccessible{0.0};
  double undecided{0.0};

  // What the boundary was made of at this diameter, kept to be looked at rather than plotted.
  std::size_t numberOfArcs{0};
  std::size_t cuspedArcs{0};
  std::size_t numberOfVertices{0};
  std::size_t clippedVertices{0};
  std::size_t degenerateVertices{0};
};

// A pore size that carries volume of its own.
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

  // The integral of the continuous part of the distribution, and the spikes, whose weights make up the rest.
  //
  // These come by two routes with no term in common: the first is the derivative, integrated over the
  // diameters; the second is what the volume lost between one diameter and the next that the derivative does
  // not account for. Together they have to come to one, since the void has some pore size everywhere, and how
  // near they come is the check of the volume against its own derivative over the whole range.
  double integral{0.0};
  double singularWeight{0.0};
  std::vector<PoreSizeSpike> spikes;

  // What is left at the end of the range, which ought to be nothing: a range that stops short of the largest
  // sphere in the void leaves volume unaccounted, and this is how much.
  double truncatedWeight{0.0};

  // The largest diameter at which any volume is left, which is the diameter of the largest sphere that fits
  // in the void. It is the position of the last spike, so it is resolved to that spike's bracket and not to
  // the spacing of the rows.
  double largestDiameter{0.0};

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
// `build` makes the accessibility for a given probe radius: it is a callback because which diagram the pores
// are taken from is the caller's business and the geometry here is the same either way.
export PoreSizeDistributionCurve exactPoreSizeDistribution(
    const std::function<VoronoiAccessibility(double)>& build, double cellVolume, double maximumDiameter,
    std::size_t numberOfBins, std::size_t subdivisions, std::size_t refinements = 12);
